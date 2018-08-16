package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.utils.VigorUtils;

@Service
public class VirusGenomeService {


    /**
     * @param sequence
     * @return List<Range>: Only ranges having length greater or equal to
     * min_gap_length will be considered as a sequence gap
     */
    public static List<Range> findSequenceGapRanges ( Integer minGapLength, NucleotideSequence sequence ) {

        List<Range> rangesOfNs = sequence.getRangesOfNs();

        List<Range> filteredRangesOfNs = new ArrayList<Range>();
        Range previousRange = Range.of(0, 0);
        for (Range currentRange: rangesOfNs) {
            if (currentRange.getLength() >= minGapLength) {
                if (previousRange.getBegin() != 0 && previousRange.getEnd() != 0) {
                    if (previousRange.getEnd() <= currentRange.getBegin() + 6) {
                        previousRange = Range.of(previousRange.getBegin(), currentRange.getEnd());
                    }
                }
            }
            filteredRangesOfNs.add(previousRange);
        }
        return rangesOfNs;
        //return filteredRangesOfNs;
    }

    public static Map<Frame, List<Long>> findInternalStops ( NucleotideSequence NTSequence ) {

        Map<Frame, List<Long>> stops = IupacTranslationTables.STANDARD.findStops(NTSequence);
        stops = VigorFunctionalUtils.frameToSequenceFrame(stops);
        return stops;
    }

     public static VirusGenome fastaRecordToVirusGenome(NucleotideFastaRecord record, VigorConfiguration config) {
        VirusGenome virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(),
                                                  config.getOrDefault(ConfigurationParameters.CompleteGene, false),
                                                  config.getOrDefault(ConfigurationParameters.CircularGene, false));
        Integer min_gap_length = config.get(ConfigurationParameters.SequenceGapMinimumLength);
        virusGenome.setInternalStops(findInternalStops(virusGenome.getSequence()));
        virusGenome.setSequenceGaps(findSequenceGapRanges(min_gap_length,virusGenome.getSequence()));
        return virusGenome;
    }
}


