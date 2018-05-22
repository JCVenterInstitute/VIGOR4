package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.align.pairwise.PairwiseAlignmentBuilder;
import org.jcvi.jillion.align.pairwise.ProteinPairwiseSequenceAlignment;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.springframework.stereotype.Service;

@Service
public class CheckCoverage implements EvaluateModel {

	@Override
	public Model evaluate(Model model,VigorForm form) {
		NucleotideSequence cds = determineCDS(model);
		List<Range> internalStops = getInternalStops(model);
		if(internalStops.size()>0){
		    model.setPseudogene(true);
        }
	    model = determineHomology(model,cds);
        model = addWeightageToScores(model,form);
		return model;
	}
	public Model addWeightageToScores(Model model,VigorForm form) {
        Map<String, Double> scores = model.getScores();
        VigorConfiguration vigorConfiguration = form.getConfiguration();
        //retrieve the parameters and multiply to the scores. define the parameters in the .ini files.
        return model;
    }
    public Model determineHomology(Model model, NucleotideSequence cds){
        long replacementOffset=0;
        if(model.getReplaceStopCodonRange()!=null){
            replacementOffset = model.getReplaceStopCodonRange().getBegin();
        }
        if(replacementOffset != 0){
            replacementOffset = getTranslatedProteinCooridnate(model.getExons(),replacementOffset,model.getInsertRNAEditingRange());
        }
        AminoAcid replacementAA = model.getAlignment().getViralProtein().getGeneAttributes().getStopTranslationException().getReplacementAA();
        ProteinSequence translatedSeq = IupacTranslationTables.STANDARD.translate(cds);
        ProteinSequenceBuilder proteinSeqBuilder = new ProteinSequenceBuilder(translatedSeq);
        if(replacementOffset !=0 && replacementAA !=null){
            proteinSeqBuilder.replace((int)replacementOffset,replacementAA);
        }
        ProteinSequence querySeq = proteinSeqBuilder.build();
        model.setTanslatedSeq(querySeq);
       	ProteinSequence subSeq = model.getAlignment().getViralProtein().getSequence();
    	AminoAcidSubstitutionMatrix blosom50 = BlosumMatrices.blosum50();
    	ProteinPairwiseSequenceAlignment actual = PairwiseAlignmentBuilder
				.createProtienAlignmentBuilder(querySeq,
						subSeq, blosom50).gapPenalty(-8, -8)
				.build();
   	    Map<String,Double> scores = new HashMap<String,Double>();
   	    if(model.getScores()!=null) {
            scores.putAll(model.getScores());
        }
        double percentIdentity = actual.getPercentIdentity()*100;
        int mismatches = actual.getNumberOfMismatches()+actual.getNumberOfGapOpenings();
        int matches = actual.getAlignmentLength()-mismatches;
        long maxSeqLength = Long.max(querySeq.getLength(),subSeq.getLength());
        double percentSimilarity = ((double)matches/maxSeqLength)*100;
        long coverage = Long.max(actual.getQueryRange().getLength(),actual.getSubjectRange().getLength());
        double percentCoverage = ((double)coverage/querySeq.getLength())*100;
        scores.put("%identity",percentIdentity);
        scores.put("%similarity",percentSimilarity);
        scores.put("%coverage",percentCoverage);
        double modelScore = percentIdentity+percentSimilarity+percentCoverage;
        scores.put("modelScore",modelScore);
        model.setScores(scores);
   	    return model;
    }

    public NucleotideSequence determineCDS(Model model){
        model.getExons().sort(Exon.Comparators.Ascending);
        List<Exon> exons = model.getExons();
        NucleotideSequence virusGenomeSeq = model.getAlignment().getVirusGenome().getSequence();
        boolean inserted=false;
        NucleotideSequenceBuilder NTSeqBuilder=new NucleotideSequenceBuilder("");
        NucleotideSequence NTSeq;
        for(Exon exon : exons){
            if(!inserted && model.getInsertRNAEditingRange()!=null && model.getInsertRNAEditingRange().getBegin()==(exon.getRange().getEnd()+1)){
                NTSeqBuilder.append(virusGenomeSeq.toBuilder(exon.getRange()));
                String insertionString = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing().getInsertionString();
                    NTSeqBuilder.append(insertionString);
                    inserted=true;
            }else {
                NTSeqBuilder.append(virusGenomeSeq.toBuilder(exon.getRange()));
            }
        }
        NTSeq = NTSeqBuilder.build();
        model.setCds(NTSeq);
        return NTSeq;
    }

	public long getTranslatedProteinCooridnate(List<Exon> exons,long NTOffset,Range insertionRange){
		long translatedProteinLength=0;
		long proteinCoordinate=0;
		Frame adjustedFrame=null;
		for(Exon exon:exons){
            long insertionLength=0;
            Range exonRange = exon.getRange();
            long length = exonRange.getLength();
		    if(insertionRange!=null ){
		        if(insertionRange.getBegin()==exon.getRange().getEnd()+1){
		            insertionLength=insertionRange.getLength();
                   long  tempLength = length-exon.getFrame().getFrame()+1+insertionLength;
                    long leftOvers = tempLength % 3;
                    if(leftOvers==1){
                        adjustedFrame=Frame.THREE;
                    }else if(leftOvers==2){
                        adjustedFrame=Frame.TWO;
                    }else if(leftOvers==0){
                        adjustedFrame=Frame.ONE;
                    }

                }
            }
		    if(exonRange.intersects(Range.of(NTOffset))){
		        Frame exonFrame = exon.getFrame();
		        if(adjustedFrame!=null){
		            exonFrame=adjustedFrame;
                }
		        long startCoordinate = exonRange.getBegin()+exonFrame.getFrame()-1;
		        long difference = NTOffset-startCoordinate;
		        if(difference<0){
		            proteinCoordinate=proteinCoordinate+translatedProteinLength;
		            break;
                }else {
		            float temp = ((float)difference)/3;
                    proteinCoordinate = (long)Math.ceil(temp);
                    proteinCoordinate = proteinCoordinate + translatedProteinLength;
                    break;
                }
		    }

		    length = length-exon.getFrame().getFrame()+1+insertionLength;
		    float temp = ((float)length)/3;
		    translatedProteinLength = (long)(translatedProteinLength+Math.ceil(temp));
		}

		return proteinCoordinate;
	}
	

	public List<Range> getInternalStops(Model model){
	    List<Range> internalStops = new ArrayList<Range>();
        NucleotideSequence virusGenomeSeq = model.getAlignment().getVirusGenome().getSequence();
        NucleotideSequenceBuilder NTSeqBuilder=new NucleotideSequenceBuilder("");
        for(Exon exon : model.getExons()){
            NTSeqBuilder.append(virusGenomeSeq.toBuilder(exon.getRange()));
        }
		NucleotideSequence cds = NTSeqBuilder.build();
        Map<Frame,List<Long>> stops = IupacTranslationTables.STANDARD.findStops(cds);
        Frame fFrame = model.getExons().get(0).getFrame();
		for(Map.Entry<Frame,List<Long>> pair :stops.entrySet()){
			if(pair.getKey().equals(fFrame)){
			List<Long> cdsStops = pair.getValue();
			for(Long stop : cdsStops){
				long NTStop = VigorFunctionalUtils.getNTRange(model.getExons(), stop);
				Range NTStopRange= Range.of(NTStop,NTStop+2);
				//if(model.getAlignment().getViralProtein().getProteinID().equals("399240871_NSP")){
				/*System.out.println(virusGenomeSeq.toBuilder().trim(NTStopRange).build());
				System.out.println(NTStopRange);
				System.out.println(model.getReplaceStopCodonRange());//}*/
				if(model.getReplaceStopCodonRange()!=null&&!NTStopRange.equals(model.getReplaceStopCodonRange()) && !Range.of(stop).equals(Range.of(cds.getLength()-3))){
				        internalStops.add(NTStopRange);
                }else if(!Range.of(stop).equals(Range.of(cds.getLength()-3))){
                    internalStops.add(NTStopRange);
                }
			}
		}
		}
		return internalStops;
	}
}
