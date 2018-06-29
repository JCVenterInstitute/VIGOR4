package org.jcvi.vigor.service;

import com.google.common.collect.Sets;
import com.google.common.graph.GraphBuilder;
import com.google.common.graph.MutableGraph;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.align.pairwise.PairwiseAlignmentBuilder;
import org.jcvi.jillion.align.pairwise.ProteinPairwiseSequenceAlignment;
import org.jcvi.jillion.core.DirectedRange;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStore;
;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.jcvi.vigor.component.MaturePeptideMatch;
import org.jcvi.vigor.component.PartialProteinSequence;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.SequenceUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Service
public class PeptideService implements PeptideMatchingService {

    /**
     * Maximum difference in amino acids length at beginning or end
     * for a peptide to be considered situated in the same spot on
     * the protein
     */
    private static final long PROXIMITY_MAX = 10L;
    private static final long MAX_GAP = 5L;
    private static Pattern productPattern = Pattern.compile("product\\s*=\\s*\"(?<product>[^\"]+)\"");
    private static Logger LOGGER = LogManager.getLogger(PeptideService.class);

    public static class PeptideMatch {
        public final ProteinFastaRecord peptide;
        public final PartialProteinSequence protein;
        public final ProteinPairwiseSequenceAlignment alignment;
        private Scores scores = Scores.of(-1,-1,-1);

        public PeptideMatch(ProteinFastaRecord peptide, PartialProteinSequence protein, ProteinPairwiseSequenceAlignment alignment) {
            this.peptide = peptide;
            this.protein = protein;
            this.alignment = alignment;
        }

        public static PeptideMatch of(ProteinFastaRecord peptide, PartialProteinSequence protein, ProteinPairwiseSequenceAlignment alignment) {
            return new PeptideMatch(peptide, protein, alignment);
        }

        public Scores getScores() {
            return scores;
        }

        public void setScores(Scores scores) {
            this.scores = scores;
        }

        public String toString() {
            Range r = alignment.getSubjectRange().getRange();
            return String.format("%s-%s %%i %.04f %%s %.04f %%c %.04f %s",
                    r.getBegin(), r.getEnd(),
                    scores.identity,
                    scores.similarity,
                    scores.identity,
                    peptide.getId(), peptide.getComment());
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
        // Return histogram of Kmers of length 1,2 and 3 to occurrences in the sequence.
        private static Map<String,Long> getKmerProfile(ProteinSequence sequence){
            String sequenceAsString = sequence.toString();
            int sequenceLength = sequenceAsString.length();

            // If the sequence is not evenly partitionable, we can ignore the extra at the end of the sequence,
            // because it would have already been counted by a smaller kmer.
            return Stream.of(1,2,3)
                  .flatMap(x->Stream.iterate(0, y -> y+x)
                                    .limit(sequenceLength/x )
                                    .map(i->sequenceAsString.substring(i,i+x)))
                  .collect(Collectors.groupingBy(String::valueOf, Collectors.counting()));
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


    static final Comparator<PeptideMatch> bySubjectRange = (a, b) -> {
        DirectedRange aDirectedRange = a.alignment.getSubjectRange();
        DirectedRange bDirectedRange = b.alignment.getSubjectRange();

        int result = aDirectedRange.getDirection().compareTo(bDirectedRange.getDirection());
        if (result == 0) {
            result = Range.Comparators.ARRIVAL.compare(aDirectedRange.getRange(),
                    bDirectedRange.getRange());
        }
        return result;
    };

    @Autowired
    public PeptideService() {

    }

    private static String formatMatchForLogging(PeptideMatch match) {
        Scores matchScores = match.getScores();
        ProteinPairwiseSequenceAlignment alignment = match.alignment;

        return String.format("%-20s     %04d    %04d-%04d     %04d-%04d   %.04f   %.04f   %.04f   %s",
                match.peptide.getId(),
                match.peptide.getLength(),
                alignment.getSubjectRange().getBegin(),
                alignment.getSubjectRange().getEnd(),
                alignment.getQueryRange().getBegin(),
                alignment.getQueryRange().getEnd(),
                matchScores.identity,
                matchScores.similarity,
                matchScores.coverage,
                match.peptide.getComment()
                );
    }

    private static String formatMatchForLogging(MaturePeptideMatch match) {
        return String.format("%-20s     %04d    %04d-%04d    %04d-%04d   %s",
                match.getReference().getProteinID(),
                match.getProtein().getLength(),
                match.getProteinRange().getBegin(),
                match.getProteinRange().getEnd(),
                match.getReferenceRange().getBegin(),
                match.getReferenceRange().getEnd(),
                match.getReference().getDefline()
        );

    }
    @Override
    public List<MaturePeptideMatch> findPeptides(PartialProteinSequence partialProtein, File peptideDatabase, Scores minscores) throws ServiceException {

        long lastAAIndex = partialProtein.getSequence().getLength() -1;
        if (partialProtein.getSequence().get(lastAAIndex) == AminoAcid.STOP) {
            partialProtein = PartialProteinSequence.of(partialProtein.getSequence().toBuilder().delete(Range.of(lastAAIndex)).build(),
                                                       partialProtein.isPartial3p(),
                                                       partialProtein.isPartial5p());
        }

        Predicate<PeptideMatch> filterByScore = match -> {

            LOGGER.debug(formatMatchForLogging(match));

            Scores matchScores = match.getScores();
            return matchScores.identity >= minscores.identity &&
                    matchScores.similarity >= minscores.similarity &&
                    matchScores.coverage >= minscores.coverage;
        };


        try (Stream<PeptideMatch> alignments = getAlignments(partialProtein, peptideDatabase).peek(m -> m.setScores(getMatchScores(m)))) {
            LOGGER.info(String.format("%-20s     %-4s    %-9s     %-9s   %-6s   %-6s   %-6s   %s",
                                       "id","len","sub","qry","%id","%sim","%cov","comment"));
            List<PeptideMatch> bestMatches = getBestMatches(alignments.sorted(bySubjectRange::compare)
                                                                      .filter(filterByScore));

            LOGGER.info( bestMatches.stream().map(m -> formatMatchForLogging(m)).collect(Collectors.joining("\n","Best matches:\n","\n")));
            List<MaturePeptideMatch> peptides = new ArrayList<>(bestMatches.size());
            MaturePeptideMatch prev = null;
            MaturePeptideMatch current = null;
            Range currentRange;
            Range previousRange = null;
            for (PeptideMatch currentMatch: bestMatches) {
                // handle first one with a previous range before the start
                if (previousRange == null) {
                    previousRange = Range.of(Range.CoordinateSystem.RESIDUE_BASED, 0);
                }

                currentRange = currentMatch.alignment.getSubjectRange().getRange();
                current = peptideFromMatch(currentMatch);

                if (currentRange.getBegin() - previousRange.getEnd() > MAX_GAP) {
                    LOGGER.debug("difference between range [{}:{}] and [{}:{}] > max gap {}. Setting fuzzy edges",
                            previousRange.getBegin(), previousRange.getEnd(),
                            currentRange.getBegin(), currentRange.getEnd(),
                            MAX_GAP);
                    if (prev != null) {
                        prev.setFuzzyEnd(true);
                    }
                    current.setFuzzyBegin(true);
                } else if (currentRange.getBegin() != previousRange.getEnd() + 1) {
                    LOGGER.debug("Range [{}:{}] doesn't align to following range [{}:{}]. Adjusting edges",
                            previousRange.getBegin(), previousRange.getEnd(),
                            currentRange.getBegin(), currentRange.getEnd());
                    // the first one we just set to the beginning
                    if (prev == null) {
                        current.setProteinRange(currentRange.toBuilder().setBegin(0).build());
                    } else {
                        adjustPeptideEdges(partialProtein, prev, current);
                        // TODO this removes zero length, but what about short ranges > 0?
                        if (prev.getProteinRange().getLength() == 0) {
                            peptides.remove(peptides.size() - 1);
                        }
                        if (current.getProteinRange().getLength() == 0) {
                            continue;
                        }
                    }
                }
                if (partialProtein.isPartial5p() && current.getProteinRange().getBegin() == 0) {
                    current.setFuzzyBegin(true);
                }
                if (partialProtein.isPartial3p() && current.getProteinRange().getEnd() == partialProtein.getSequence().getLength() -1) {
                    current.setFuzzyEnd(true);
                }
                peptides.add(current);
                prev = current;
                previousRange = current.getProteinRange();
            }
            if (! peptides.isEmpty()) {
                // if last peptide is close to end, adjust to end
                MaturePeptideMatch lastPeptide = peptides.get(peptides.size() -1);
                if (lastPeptide.getProteinRange().getEnd() + MAX_GAP >= partialProtein.getSequence().getLength()) {
                    lastPeptide.setProteinRange(lastPeptide.getProteinRange()
                                                           .toBuilder()
                                                           .setEnd(partialProtein.getSequence().getLength() -1)
                                                           .build());
                }
            }
            LOGGER.debug(peptides.stream()
            .map(m -> formatMatchForLogging(m))
            .collect(Collectors.joining("\n","After adjusting edges:\n","\n")));

            return peptides;
        } catch (IOException e) {
            throw new ServiceException(String.format("Problem finding peptide matches for sequence %s in database %s. got %s: %s",
                    partialProtein, peptideDatabase, e.getClass().getSimpleName(), e.getMessage()), e);
        }
    }

    private MutableGraph<PeptideMatch> matchesToGraph(List<PeptideMatch> matches) {

        MutableGraph<PeptideMatch> graph = GraphBuilder.directed().build();
        PeptideMatch current;
        PeptideMatch next;
        Range currentRange;
        Range nextRange;
        for (int i=0; i < matches.size(); i++) {
            current = matches.get(i);
            graph.addNode(current);
            currentRange = current.alignment.getSubjectRange().asRange();
            for (int j=i+1; j < matches.size(); j++) {
                next = matches.get(j);
                graph.addNode(next);
                nextRange = next.alignment.getSubjectRange().asRange();
                if ((double) nextRange.intersection(currentRange).getLength()/ (double)Math.min(nextRange.getLength(), currentRange.getLength()) <= .1d &&
                        Math.abs(nextRange.getBegin() - currentRange.getEnd()) < PROXIMITY_MAX) {
                    graph.putEdge(current, next);
                }
                // stop if the ranges are too far away from the current range.
                if (nextRange.getBegin() > currentRange.getEnd() + MAX_GAP) {
                    break;
                }
            }
        }
        return graph;
    }
    private List<PeptideMatch> getBestMatches(Stream<PeptideMatch> alignments) {
        List<PeptideMatch> matches = alignments.collect(Collectors.toList());
        MutableGraph<PeptideMatch> graph = matchesToGraph(matches);

        List<List<PeptideMatch>> paths = new ArrayList<>();
        List<PeptideMatch> starts = matches.stream().filter(m -> graph.inDegree(m) == 0).collect(Collectors.toList());
        for (PeptideMatch start: starts) {
            paths.addAll(findPaths(graph, start));
        }
        Function<Scores,Double> sumScores = (s) -> s.identity + s.similarity + s.coverage;
        Function<List<PeptideMatch>, Double> scorePath = a -> {
            //
            if (a.stream().map(m -> extractProduct(m.peptide.getComment())).distinct().count() != a.size()) {
                return 0d;
            }
            double scoreSum = a.stream()
             .collect(Collectors.summingDouble(m->sumScores.apply(m.getScores())));

            double coverageScore = 0;
            if (! a.isEmpty()) {
                coverageScore = (double) a.stream().collect(Collectors.summingLong(m -> m.alignment.getSubjectRange().getLength())) / (double) a.get(0).protein.getSequence().getLength();
            }
            return scoreSum +  coverageScore;
        };

        return paths.stream()
                    .max( Comparator.comparing( (a) -> scorePath.apply(a)))
                    .orElse(Collections.EMPTY_LIST);
    }

    private List<List<PeptideMatch>> findPaths(MutableGraph<PeptideMatch> graph, PeptideMatch node) {
        List<PeptideMatch> path = new ArrayList<>();
        path.add(node);
        List<PeptideMatch> tempPath;
        List<List<PeptideMatch>> paths = new ArrayList<>();
        for (PeptideMatch nextNode: graph.successors(node)) {
            for (List<PeptideMatch> nextPath: findPaths(graph, nextNode)) {
                tempPath = new ArrayList<>(path);
                tempPath.addAll(nextPath);
                paths.add(tempPath);
            }
        }
        if (paths.isEmpty()) {
            paths.add(path);
        }
        return paths;
    }

    // TODO return new MaturePeptideMatch objects rather than altering the existing ones.
    private void adjustPeptideEdges(PartialProteinSequence partialSubjectSequence, MaturePeptideMatch prev, MaturePeptideMatch current) {

        ProteinSequence subjectSequence = partialSubjectSequence.getSequence();
        Range previousRange = prev.getProteinRange();
        Range currentRange = current.getProteinRange();
        LOGGER.debug(() -> String.format("adjusting edges for\n[%s-%s] %s\n[%s-%s] %s",
                previousRange.getBegin(), previousRange.getEnd(),
                SequenceUtils.elipsedSequenceString(prev.getProtein().toBuilder().trim(prev.getProteinRange()).build(),30,30),
                currentRange.getBegin(),currentRange.getEnd(),
                SequenceUtils.elipsedSequenceString(current.getProtein().toBuilder().trim(current.getProteinRange()).build(), 30, 30)
                ));

        PeptideProfile previousReferenceProfile = PeptideProfile.profileFromSequence(prev.getReference().getSequence());
        PeptideProfile currentReferenceProfile = PeptideProfile.profileFromSequence(current.getReference().getSequence());

        long previousEnd = previousRange.getEnd();
        long prevReferenceLength = prev.getReference().getSequence().getLength();

        if (prev.getReferenceRange().getEnd() < prevReferenceLength -1) {
            previousEnd += ((prevReferenceLength - 1) - prev.getReferenceRange().getEnd());
        }

        long currentBegin = currentRange.getBegin();
        if (current.getReferenceRange().getBegin() > 0) {
            currentBegin -= current.getReferenceRange().getBegin();
        }
        long start, end;
        if (previousEnd >= currentBegin) {
            // Overlap
            start = currentBegin - 1;
            end = Math.min(currentRange.getEnd(), previousEnd);
        } else {
            // Gap
            start = previousEnd;
            end = currentBegin;
        }
        Range testPreviousRange;
        Range testCurrentRange;
        ProteinSequence testPrevious;
        ProteinSequence testCurrent;
        Range[] bestRange = {previousRange, currentRange};
        double bestScore = 0;
        double testScore = 0;
        double previousWeight = 1.0;
        double currentWeight = 1.0;
        for (; start < end; start++) {
            previousWeight = 1.0;
            currentWeight = 1.0;

            // profile and score new sequences
            testPreviousRange = previousRange.toBuilder().setEnd(start).build();
            testCurrentRange = currentRange.toBuilder().setBegin(start+1).build();
            // now we need the sequence for the new test ranges
            testPrevious = subjectSequence.toBuilder().trim(testPreviousRange).build();
            testCurrent = subjectSequence.toBuilder().trim(testCurrentRange).build();
            if (start > previousEnd) {
                previousWeight = .9;
            } else if (start == previousEnd) {
                previousWeight = 1.1;
            }
            if (start +1 < currentBegin) {
                currentWeight = .9;
            } else if (start + 1 == currentBegin) {
                currentWeight = 1.1;
            }
            testScore = scorePeptideByProfile(testPrevious, previousReferenceProfile) * previousWeight + scorePeptideByProfile(testCurrent, currentReferenceProfile) * currentWeight;
            LOGGER.trace("checking {}-{} and {}-{} got score {}",
                    previousRange.getBegin(),start,
                    start+1, currentRange.getEnd(),
                    testScore);
            if (testScore > bestScore) {
                bestRange[0] = testPreviousRange;
                bestRange[1] = testCurrentRange;
                bestScore = testScore;
            }
        }

        Range bestPreviousRange = bestRange[0];
        Range bestCurrentRange = bestRange[1];
        LOGGER.debug("best ranges after adjustment {}-{} and {}-{}",
                bestPreviousRange.getBegin(),bestPreviousRange.getEnd(),
                bestCurrentRange.getBegin(), bestCurrentRange.getEnd());

        // adjust previous and current
        prev.setProteinRange(bestPreviousRange);
        Range referenceRange = prev.getReferenceRange();
        prev.setReferenceRange(referenceRange.toBuilder().setEnd(referenceRange.getEnd() + (bestPreviousRange.getEnd() - previousEnd)).build());
        // TODO is this correct?
        prev.setFuzzyEnd(false);

        current.setProteinRange(bestCurrentRange);
        referenceRange = current.getReferenceRange();
        current.setReferenceRange(referenceRange.toBuilder().setBegin(referenceRange.getBegin() + (bestCurrentRange.getBegin() - currentBegin)).build());
        // TODO
        current.setFuzzyBegin(false);
    }

    private MaturePeptideMatch peptideFromMatch (PeptideMatch match) {

        ViralProtein referenceProtein = new ViralProtein();
        referenceProtein.setSequence(match.peptide.getSequence());
        referenceProtein.setDefline(String.join(" ", ">" + match.peptide.getId(), match.peptide.getComment()));
        referenceProtein.setProduct(extractProduct(referenceProtein.getDefline()));
        referenceProtein.setProteinID(match.peptide.getId());


        return MaturePeptideMatch.of(match.protein.getSequence(),
                                     referenceProtein,
                                     match.alignment.getSubjectRange().asRange(),
                                     match.alignment.getQueryRange().asRange(),
                                     false,
                                     false,
                                     match.scores.identity,
                                     match.scores.similarity,
                                     match.scores.coverage
                                     );
    }

    private static double computeSimilarity(double numberOfGaps, double peptideLength, double alignmentLength, double mismatchCount) {
           return   ((alignmentLength - numberOfGaps - mismatchCount) /  peptideLength);
    }

    private static long getComparisonPeptideLength(PartialProteinSequence protein, Range proteinRange, ProteinFastaRecord peptide, Range peptideRange) {
        long peptideLength = peptide.getSequence().getLength();
        if (protein.isPartial5p() || protein.isPartial3p()) {

            if (protein.isPartial5p() && peptideRange.getBegin() > proteinRange.getBegin()) {
                peptideLength -= Math.abs(proteinRange.getBegin() - peptideRange.getBegin());
            }
            if (protein.isPartial3p()) {
                long peptideUnmatchedAfter = ((peptide.getLength() - 1) - peptideRange.getEnd()) - ((protein.getSequence().getLength() - 1) - proteinRange.getEnd());
                if (peptideUnmatchedAfter > 0) {
                    peptideLength -= peptideUnmatchedAfter;
                }
            }
            if (peptideLength != peptide.getLength()) {
                LOGGER.trace("For peptide {} alignment length {}, using comparison length {} rather than {}", peptide.getId(), proteinRange.getLength(), peptideLength, peptide.getSequence().getLength());
            }
        }
        return peptideLength;
    }

    private static double computeCoverage(PartialProteinSequence protein, Range proteinRange, ProteinFastaRecord peptide, Range peptideRange) {
        long peptideLength = getComparisonPeptideLength(protein, proteinRange, peptide, peptideRange);
        return Math.max(
                (double) proteinRange.getLength() / (double) protein.getSequence().getLength(),
                (double) peptideRange.getLength() / (double) peptideLength);
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

        double score = (referenceProfile.getPeptideLength() * (2 * matches)) / ((Math.abs(subjectProfile.getPeptideLength() - referenceProfile.getPeptideLength())) + ((2 * matches) + mismatches));
        return score;
    }

    Stream<PeptideMatch> getAlignments(PartialProteinSequence protein, File peptideDatabase) throws IOException {

        LOGGER.info("finding alignments in {} for seq {}", peptideDatabase, SequenceUtils.elipsedSequenceString(protein.getSequence(), 40,20));

        ProteinFastaFileDataStore peptideDataStore = ProteinFastaFileDataStore.fromFile(peptideDatabase);
        // TODO configurable gap penalties and blosum matrix
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

    private static Scores getMatchScores(PeptideMatch match) {
        return Scores.of(match.alignment.getPercentIdentity(),
                         computeCoverage(match.protein,
                                         match.alignment.getSubjectRange().asRange(),
                                         match.peptide,
                                         match.alignment.getQueryRange().asRange()),
                         computeSimilarity(
                        match.alignment.getNumberOfGapOpenings(),
                        getComparisonPeptideLength(match.protein,
                                                   match.alignment.getSubjectRange().asRange(),
                                                   match.peptide,
                                                   match.alignment.getQueryRange().asRange()),
                        match.alignment.getAlignmentLength(),
                        SequenceUtils.computeMismatches(match.alignment.getGappedQueryAlignment(),
                                match.alignment.getGappedSubjectAlignment(),
                                BlosumMatrices.blosum40())
                )
        );

    }

    public static String extractProduct(String defline) {
        Matcher m = productPattern.matcher(defline);

        if (m.find()) {
            return m.group("product");
        }
        return "";
    }

}
