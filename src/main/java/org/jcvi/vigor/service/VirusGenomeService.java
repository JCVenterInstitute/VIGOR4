package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

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

@Service
public class VirusGenomeService {

    /**
     * @param sequence
     * @return List<Range>: Only ranges having length greater or equal to
     * min_gap_length will be considered as a sequence gap
     */
    public static List<Range> findSequenceGapRanges ( Integer minGapLength, NucleotideSequence sequence ) {

        List<Range> rangesOfNs = sequence.getRangesOfNs();
        rangesOfNs=rangesOfNs.stream().filter(x->x.getLength()>=minGapLength).collect(Collectors.toList());
        List<Range> filteredRangesOfNs = new ArrayList<Range>();
        if (!rangesOfNs.isEmpty()) {
            Range previousRange = Range.of(0,0);
            for (int i = 0; i < rangesOfNs.size(); i++) {
                if(i==0){
                    filteredRangesOfNs.add(rangesOfNs.get(0));
                    previousRange=rangesOfNs.get(0);
                }
                Range currentRange = rangesOfNs.get(i);
                    if (i!=0) {
                        if (Range.of(previousRange.getEnd(),currentRange.getBegin()).getLength()<=6) {
                             filteredRangesOfNs.set(i-1, Range.of(previousRange.getBegin(), currentRange.getEnd()));
                        }else{
                            filteredRangesOfNs.add(currentRange);
                        }
                    }
                    previousRange=filteredRangesOfNs.get(filteredRangesOfNs.size()-1);

            }
        }
        return filteredRangesOfNs;
    }

    public static Map<Frame, List<Long>> findInternalStops ( NucleotideSequence NTSequence ) {

        Map<Frame, List<Long>> stops = IupacTranslationTables.STANDARD.findStops(NTSequence);
        stops = VigorFunctionalUtils.frameToSequenceFrame(stops);
        return stops;
    }

    public static VirusGenome fastaRecordToVirusGenome( NucleotideFastaRecord record, VigorConfiguration config) {
        VirusGenome virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(),
                config.getOrDefault(ConfigurationParameters.CompleteGene, false),
                config.getOrDefault(ConfigurationParameters.CircularGene, false));
        Integer min_gap_length = config.get(ConfigurationParameters.SequenceGapMinimumLength);
        virusGenome.setInternalStops(findInternalStops(virusGenome.getSequence()));
        virusGenome.setSequenceGaps(findSequenceGapRanges(min_gap_length,virusGenome.getSequence()));
        return virusGenome;
    }

    public static VirusGenome reverseComplementVirusGenome(VirusGenome inputGenome,VigorConfiguration config){
        NucleotideSequence reverseCompGenome = inputGenome.getSequence().toBuilder().reverseComplement().build();
        VirusGenome virusGenome = new VirusGenome(reverseCompGenome, inputGenome.getDefline(), inputGenome.getId(),
                inputGenome.getIsComplete(), inputGenome.getIsCircular());
        Integer min_gap_length = config.get(ConfigurationParameters.SequenceGapMinimumLength);
        virusGenome.setInternalStops(findInternalStops(virusGenome.getSequence()));
        virusGenome.setSequenceGaps(findSequenceGapRanges(min_gap_length,virusGenome.getSequence()));
        return virusGenome;
    }

}


