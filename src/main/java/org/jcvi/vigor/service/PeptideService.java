package org.jcvi.vigor.service;

import com.google.common.collect.Sets;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.align.pairwise.PairwiseAlignmentBuilder;
import org.jcvi.jillion.align.pairwise.ProteinPairwiseSequenceAlignment;
import org.jcvi.jillion.core.DirectedRange;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStore;
;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.service.exception.ServiceException;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Service
public class PeptideService implements PeptideMatchingService {

    private static final long PROXIMITY_MIN = 30L;
    private static Logger LOGGER = LogManager.getLogger(PeptideService.class);

    private static class Scores {
        final double minidentity;
        final double mincoverage;
        final double minsimilarity;

        public Scores(double minidentity, double mincoverage, double minsimilarity) {
            this.minidentity = minidentity;
            this.mincoverage = mincoverage;
            this.minsimilarity = minsimilarity;
        }

        public static Scores of(double identity, double coverage, double similarity) {
            return new Scores(identity, coverage, similarity);
        }
    }


    private static class PeptideMatch {
        public final ProteinFastaRecord peptide;
        public final ViralProtein protein;
        public final ProteinPairwiseSequenceAlignment alignment;

        public PeptideMatch(ProteinFastaRecord peptide, ViralProtein protein, ProteinPairwiseSequenceAlignment alignment) {
            this.peptide = peptide;
            this.protein = protein;
            this.alignment = alignment;
        }

        public static PeptideMatch of(ProteinFastaRecord peptide, ViralProtein protein, ProteinPairwiseSequenceAlignment alignment) {
            return new PeptideMatch(peptide, protein, alignment);
        }

    }

    private static class PeptideProfile {
        final private Map<String, Long> profile;
        final private long peptideLength;

        public PeptideProfile(Map<String, Long> profile, long length) {
            this.profile = profile;
            this.peptideLength = length;
        }

        public static PeptideProfile profileFromSequence(ProteinSequence sequence) {
            return new PeptideProfile(getKmerProfile(sequence), sequence.getLength());
        }
        //TODO using KMers for this is overkill but convenient
        private static Map<String,Long> getKmerProfile(ProteinSequence sequence){
            return Stream.of(1,2,3).flatMap(i -> sequence.kmers(i)).collect(Collectors.groupingBy(String::valueOf, Collectors.counting()));
        }

        public Long get(String kmer) {
            return profile.get(kmer);
        }

        public Long getOrDefault(String kmer, Long defaultValue) {
            return profile.getOrDefault(kmer, defaultValue);
        }

        public Set<String> keySet() {
            return profile.keySet();
        }

        public long getPeptideLength() {
            return peptideLength;
        }
    }

    static final Comparator<PeptideMatch> byQueryAlignment = new Comparator<PeptideMatch>() {
        @Override
        public int compare(PeptideMatch a, PeptideMatch b) {
            DirectedRange aDirectedRange = a.alignment.getSubjectRange();
            DirectedRange bDirectedRange = b.alignment.getSubjectRange();

            int result = aDirectedRange.getDirection().compareTo(bDirectedRange.getDirection());
            if (result == 0) {
                result = Range.Comparators.ARRIVAL.compare(aDirectedRange.getRange(),
                        bDirectedRange.getRange());
            }
            return result;
        }
    };

    @Autowired
    public PeptideService() {

    }

    public List<ProteinSequence> findPeptides(ViralProtein protein, File peptideDatabase) throws ServiceException {
        return findPeptides(protein, peptideDatabase, Scores.of(0.25d, .40d, .50d));
    }

    public List<ProteinSequence> findPeptides(ViralProtein protein, File peptideDatabase, Scores minscores) throws ServiceException {
        // outline from VIGOR3. Leading - means not implementing in Vigor4
        // - 1) fill in gaps and truncations from referenceSequence (
        // 2) find alignments (blast or jillion)
        // 3) filter hits based on
        //     a) overlap
        //     b) coverage
        //     c) similarity
        //     d) identity
        //     e) other?
        // 4) weight scores based on edge adjustments made
        // 5) remove redundant peptides
        // 6) create tree where nodes are connected if end and start are within a epsilon gap
        // 7) find best path for each start
        // 8) remove redundant paths
        // 9) adjust peptide edges
        //

        // filter
        Predicate<PeptideMatch> filterByScore = match -> {
            ProteinPairwiseSequenceAlignment alignment = match.alignment;
            LOGGER.debug("alignment for subject {} to query {} %identity {} min %identity {} %similarity {} min %similarity {}",
                    alignment.getGappedSubjectAlignment(),
                    alignment.getGappedQueryAlignment(),
                    alignment.getPercentIdentity(), minscores.minidentity,
                    computeSimilarity(match), minscores.minsimilarity);

            return alignment.getPercentIdentity() >= minscores.minidentity &&
                    computeSimilarity(match) >= minscores.minsimilarity;
        };

        try (Stream<PeptideMatch> alignments = getAlignments(protein, peptideDatabase);) {


            List<Range> ranges = new ArrayList<>();

            // assumes sorted
            Function<PeptideMatch, Range> binByRange = match -> {
                Range subjectRange = match.alignment.getSubjectRange().asRange();

                if (ranges.size() != 0) {
                    Range rangeKey = ranges.get(ranges.size() - 1);

                    if (Math.abs(rangeKey.getBegin() - subjectRange.getBegin()) < PROXIMITY_MIN &&
                            Math.abs(rangeKey.getEnd() - subjectRange.getEnd()) < PROXIMITY_MIN) {
                        return rangeKey;
                    }
                }
                ranges.add(subjectRange);
                return subjectRange;
            };

            Map<Range, List<PeptideMatch>> alignmentsByRange = alignments.filter(filterByScore)
                                                                         .sorted(byQueryAlignment::compare)
                                                                         .collect(Collectors.groupingBy(binByRange));

            LOGGER.debug(() -> String.format("%s alignments after binning into %s ranges: %s",
                    alignmentsByRange.values().stream()
                                     .flatMap(l -> l.stream())
                                     .count(),
                    alignmentsByRange.keySet()
                                     .size(),
                    alignmentsByRange.keySet().stream()
                                     .sorted(Range.Comparators.ARRIVAL)
                                     .map(r -> String.format("%d-%d", r.getBegin(), r.getEnd()))
                                     .collect(Collectors.joining(", "))
            ));


            List<ProteinSequence> peptides = new ArrayList<>(alignmentsByRange.size());
            // TODO don't just use score.
            PeptideMatch match;
            for (
                    List<PeptideMatch> matches : alignmentsByRange.values())

            {
                match = matches.stream().max(Comparator.comparing(p -> p.alignment.getScore())).get();

                LOGGER.debug("With score {} %identity {} returning range {}\nS>{}\nQ>{}\nP>{}",
                        match.alignment.getScore(),
                        match.alignment.getPercentIdentity(),
                        getRangeString(match),
                        match.alignment.getGappedSubjectAlignment(),
                        match.alignment.getGappedQueryAlignment(),
                        match.peptide.getSequence());
                peptides.add(match.alignment.getGappedSubjectAlignment());
            }
            return peptides;
        } catch (IOException e) {
            throw new ServiceException(String.format("Problem finding peptide matches for sequence %s in database %s", protein, peptideDatabase), e);
        }
    }

    private static String getRangeString(PeptideMatch match) {
        Range queryRange = match.alignment.getQueryRange().asRange();
        Range subjectRange = match.alignment.getSubjectRange().asRange();
        return String.format("%s%s-%s%s",
                queryRange.getBegin() > 0 ? "<" : "", subjectRange.getBegin(Range.CoordinateSystem.RESIDUE_BASED),
                queryRange.getEnd() < match.peptide.getLength() - 1 ? ">" : "", subjectRange.getEnd(Range.CoordinateSystem.RESIDUE_BASED));
    }

    private double computeSimilarity(PeptideMatch match) {
        ProteinSequence querySequence = match.alignment.getGappedQueryAlignment();
        return   ((double)(match.alignment.getAlignmentLength() - querySequence.getNumberOfGaps() - match.alignment.getNumberOfMismatches()) / (double) match.peptide.getSequence().getLength());
    }

    private double computeCoverage(PeptideMatch match) {
        return Math.max(
                match.alignment.getQueryRange().getLength()/match.protein.getSequence().getLength(),
                match.alignment.getSubjectRange().getLength() / match.peptide.getSequence().getLength()
        );
    }


    private double scorePeptideByProfile(ProteinSequence peptide, PeptideProfile referenceProfile) {
        return scoreProfile(PeptideProfile.profileFromSequence(peptide), referenceProfile);
    }
    
    private double scoreProfile(PeptideProfile subjectProfile, PeptideProfile referenceProfile) {
        Set<String> peptideKeys = subjectProfile.keySet();
        Set<String> profileKeys = referenceProfile.keySet();
        long mismatches = 0;
        long matches = 0;
        long peptideCount;
        long profileCount;
        Set<String> commonKeys = Sets.intersection(peptideKeys, profileKeys);
        for (String key: commonKeys) {
            peptideCount = subjectProfile.get(key);
            profileCount = referenceProfile.get(key);
            matches += Math.min(peptideCount, profileCount);
            mismatches += Math.abs(peptideCount - profileCount);
        }

        //TODO account for gaps

        for (String key: Sets.symmetricDifference(peptideKeys, profileKeys)) {
            mismatches += subjectProfile.getOrDefault(key, referenceProfile.get(key));
        }

        double score = (referenceProfile.getPeptideLength() * (2 * matches)) / Math.abs(subjectProfile.getPeptideLength() - referenceProfile.getPeptideLength()) + ((2 * matches) + mismatches);
        return score;
    }

    Stream<PeptideMatch> getAlignments(ViralProtein protein, File peptideDatabase) throws IOException {

        ProteinFastaFileDataStore peptideDataStore = ProteinFastaFileDataStore.fromFile(peptideDatabase);

        return peptideDataStore.records()
                               .map(record -> PeptideMatch.of(record,
                                       protein,
                                       PairwiseAlignmentBuilder.createProtienAlignmentBuilder(
                                               record.getSequence(),
                                               protein.getSequence(),
                                               BlosumMatrices.blosum40())
                                                               .useLocalAlignment(true)
                                                               .gapPenalty(-16F,-8F)
                                                               .build())
                               );
    }


}
