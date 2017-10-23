package org.jcvi.vigor.service;

import org.jcvi.vigor.component.*;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.FormatVigorOutput;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.jcvi.vigor.utils.VigorUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by snettem on 5/16/2017.
 */
@Service
public class ModelGenerationService {

	private static boolean isDebug = false;

	@Autowired
	private GeneModelGenerationService geneModelGenerationService;
	private static final Logger LOGGER = LogManager.getLogger(ModelGenerationService.class);

	public void generateModels(List<Alignment> alignments, VigorForm form) {
		isDebug = form.isDebug();
		List<Model> candidateModels = determineCandidateModels(alignments, form);
		geneModelGenerationService.generateGeneModel(candidateModels, form);
	}

	/**
	 * 
	 * @param alignments
	 * @param form
	 * @return all the possible models of each alignment and after splitting
	 *         models at the sequence gaps
	 */
	public List<Model> determineCandidateModels(List<Alignment> alignments, VigorForm form) {
		//System.out.println("Number of alignments are " + alignments.size());
		List<Model> initialModels = new ArrayList<Model>();
		List<Model> candidateModels = new ArrayList<Model>();
		for (int i = 0; i < alignments.size(); i++) {
			initialModels.addAll(alignmentToModels(alignments.get(i), alignments.get(i).getAlignmentTool_name()));
		}
		if (isDebug) {
			System.out.println("************Initial Models*************");
			FormatVigorOutput.printModels2(initialModels);
		}
        
		List<Range> sequenceGaps=new ArrayList<Range>();
		// get sequence gaps
		if(initialModels.get(0).getAlignment().getVirusGenome().getSequenceGaps() !=null){
		sequenceGaps = initialModels.get(0).getAlignment().getVirusGenome().getSequenceGaps();
		}
		List<Range> validSequenceGaps = new ArrayList<Range>();
		String minGapLenString = "";
		
		if (sequenceGaps !=null && sequenceGaps.size() > 0) {
			minGapLenString = form.getVigorParametersList().get("min_gap_length");
			long minGapLength = 20;
			if (VigorUtils.is_Integer(minGapLenString)) {
				minGapLength = Long.parseLong(minGapLenString);
			}

			for (Range gapRange : sequenceGaps) {
				if (gapRange.getLength() >= minGapLength) {
					validSequenceGaps.add(gapRange);
				}
			}
		}

		// split models at sequence gaps
        System.out.println("Count of initial models"+ initialModels.size());
		for (Model model : initialModels) {
			Range diffRange=null;
			for(int i=0;i<model.getExons().size()-1;i++){
			 diffRange = Range.of(model.getExons().get(i).getRange().getEnd(), model.getExons().get(i+1).getRange().getBegin());
			}
			if(diffRange!=null && diffRange.getLength()>=20){
			
			List<Model> newModelsreturned = splitModelAtSequenceGaps(model,validSequenceGaps);
			//candidateModels.addAll(splitModelAtSequenceGaps(model, validSequenceGaps));
			if(newModelsreturned.size()>1){
			System.out.println("Before Splitting");
			model.getExons().stream().forEach(System.out::println);
			System.out.println("After splitting at the sequence gaps");
			newModelsreturned.stream().forEach(System.out::println);
			}
			candidateModels.addAll(newModelsreturned);
			}else
			{
				candidateModels.add(model);
			}
		}
       System.out.println("Count after splitting"+candidateModels.size());
		if (isDebug) {
			System.out.println("********After splitting models at the Genome sequence gaps**********");
			FormatVigorOutput.printModels2(candidateModels);
		}

		return candidateModels;
	}

	/**
	 * 
	 * @param alignment
	 * @param alignmentTool
	 * @return Models of each alignment.
	 */
	public List<Model> alignmentToModels(Alignment alignment, String alignmentTool) {
		Map<Direction, List<AlignmentFragment>> alignmentFragsGroupedList = alignment.getAlignmentFragments().stream()
				.collect(Collectors.groupingBy(w -> w.getDirection()));
		Set<Direction> keyset = alignmentFragsGroupedList.keySet();
		Iterator<Direction> iter = keyset.iterator();
		List<Model> models = new ArrayList<Model>();
		for (Direction direction : keyset) {
			List<List<AlignmentFragment>> ListOfCompatibleFragsList = generateCompatibleFragsChains(
					alignmentFragsGroupedList.get(iter.next()), alignmentTool);
			Iterator<List<AlignmentFragment>> iter1 = ListOfCompatibleFragsList.iterator();
			while (iter1.hasNext()) {
				List<AlignmentFragment> compatibleFragsList = (List<AlignmentFragment>) iter1.next();
				List<Exon> exons = determineVirusGenomeExons(compatibleFragsList);
				Model model = new Model();
				model.setExons(exons);
				model.setAlignment(alignment);
				model.setGeneSymbol(alignment.getViralProtein().getProteinID());
				model.setDirection(direction);
				model = generateScores(model, alignment);
				List<String> statusList = new ArrayList<String>();
				statusList.add("Initial Model");
				model.setStatus(statusList);
				models.add(model);
			}

		}

		return models;

	}

	public Model generateScores(Model model, Alignment alignment) {

		Map<String, Float> alignmentScores = alignment.getAlignmentScore();
		model.setScores(alignmentScores);

		return model;
	}

	/**
	 * 
	 * @param initModels
	 * @param form
	 * @param genome
	 * @return Models are split at sequence gaps and new list of models are
	 *         returned
	 */
	public List<Model> splitModelAtSequenceGaps(Model initModel, List<Range> validSequenceGaps) {

		List<Model> newModels = new ArrayList<Model>();
		Model model = new Model();
		model = Model.deepClone(initModel);
		if (validSequenceGaps.size() > 0) {
			Exon nextExon;
			Exon currentExon;
			Range diffRange;
			List<Exon> firstGroup = new ArrayList<Exon>();
			List<Exon> secondGroup = new ArrayList<Exon>();
			Model firstModel;
			Model secondModel=null;
				List<Exon> modelExons = new ArrayList<Exon>();
				modelExons = model.getExons();
				secondGroup.addAll(modelExons);
				boolean startExist = true;
				
				for (int j = 0; j < modelExons.size(); j++) {
										
					if (j != modelExons.size() - 1) {
						
						
						nextExon = modelExons.get(j + 1);
						currentExon = modelExons.get(j);
						diffRange = Range.of(currentExon.getRange().getEnd(), nextExon.getRange().getBegin());
						if (diffRange.getLength() >= 20) {
							firstGroup.add(modelExons.get(j));
							
							boolean temp = true;
							for (int k = 0; k < validSequenceGaps.size(); k++) {
								if (diffRange.intersects(validSequenceGaps.get(k)) && temp) {
									secondModel = new Model();
									secondModel = Model.deepClone(model);
									firstModel = new Model();
									firstModel = Model.deepClone(model);;
									if(currentExon.getRange().getEnd()<validSequenceGaps.get(k).getBegin()){
									currentExon.set_3p_adjusted(true);
									}
								    if(nextExon.getRange().getBegin()>validSequenceGaps.get(k).getEnd()){
								    	nextExon.set_5p_adjusted(true);
								    }
									List<Exon> tempFirst = new ArrayList<>();
									tempFirst.addAll(firstGroup);
									List<Exon> tempSecond = new ArrayList<>();
									tempSecond.addAll(secondGroup);
									firstModel.setExons(tempFirst);
									firstModel.setPartial3p(true);
									firstModel.getStatus().add("Model splitted at sequence gaps");
									if (!startExist) {
										firstModel.setPartial5p(true);
									}
									secondGroup.removeAll(tempFirst);
									tempSecond.removeAll(tempFirst);
									secondModel.setExons(tempSecond);
									secondModel.setPartial5p(true);
									secondModel.getStatus().add("Model splitted at sequence gaps");
									newModels.add(firstModel);
									startExist = false;
									firstGroup.clear();
									temp = false;

								}
							}
						}
					}
				}
				if (secondModel !=null && secondModel.getExons().size() > 0) {
					newModels.add(secondModel);
				}
				if(newModels.size()==0){
					newModels.add(initModel);
				}

				return newModels;
			} 
			
		 else {
			newModels.add(initModel);
			return newModels;
		}

	}

	/**
	 *
	 * @param alignmentFragments
	 *            : alignment fragments that are grouped based on direction is
	 *            input to the function
	 * @return List<List<AlignmentFragment>>: compatible(check for overlap) list
	 *         of alignment fragments and their permutations and combinations is
	 *         grouped
	 */

	public List<List<AlignmentFragment>> generateCompatibleFragsChains(List<AlignmentFragment> alignmentFragments,
			String alignmentTool) {
		List<List<AlignmentFragment>> ListOfCompatibleFragsList = new ArrayList<List<AlignmentFragment>>();
		int AAOverlapOffset = 10;
		int NTOverlapOffset = 30;
		boolean tempFlag = true;
		for (int i = 0; i < alignmentFragments.size(); i++) {
			if (alignmentFragments.get(i).isSubChain()) {
				tempFlag = generateSubChains(alignmentTool);
			}
			if (tempFlag) {
				long NTEnd = alignmentFragments.get(i).getNucleotideSeqRange().getEnd();
				long AAEnd = alignmentFragments.get(i).getProteinSeqRange().getEnd();
				List<AlignmentFragment> compatibleFragsList = new ArrayList<AlignmentFragment>();
				compatibleFragsList.add(AlignmentFragment.deepClone(alignmentFragments.get(i)));
				if (i != alignmentFragments.size() - 1) {
					for (int j = i + 1; j < alignmentFragments.size(); j++) {
						long nextNTStart = alignmentFragments.get(j).getNucleotideSeqRange().getBegin();
						long nextAAStart = alignmentFragments.get(j).getProteinSeqRange().getBegin();
						if (nextNTStart >= NTEnd - NTOverlapOffset && nextAAStart >= AAEnd - AAOverlapOffset) {
							alignmentFragments.get(j).setSubChain(true);
							compatibleFragsList.add(AlignmentFragment.deepClone(alignmentFragments.get(j)));
						}
					}
				}
				ListOfCompatibleFragsList.add(compatibleFragsList);
			
				if (compatibleFragsList.size() > 2) {
					int temp = 1;
					for (int k = 0; k < compatibleFragsList.size() - 2; k++) {
						List<AlignmentFragment> subChain = new ArrayList<AlignmentFragment>();
						for(AlignmentFragment alignFrag : compatibleFragsList){
							subChain.add(AlignmentFragment.deepClone(alignFrag));
						}
						for (int j = 1; j <= temp; j++) {
							subChain.remove(1);
						}
						ListOfCompatibleFragsList.add(subChain);
						temp++;
					}
				}
			
			}
		}
		return ListOfCompatibleFragsList;
	}

	/**
	 *
	 * @return returns true if subChains has to be generated. eg: for exonerate
	 *         this function returns false as subChains are not required to be
	 *         generated
	 */

	public boolean generateSubChains(String alignmentTool) {

		if (alignmentTool.equals("blast")) {
			return true;
		}
		return false;
	}

	/**
	 *
	 * @param compatibleFragsList
	 * @return List<Exon>: nucleotideSequence range of alignment fragment will
	 *         be set as range of exon. 5' and 3' edges of exon are not
	 *         modified.
	 *
	 */
	public List<Exon> determineVirusGenomeExons(List<AlignmentFragment> compatibleFragsList) {
		Iterator<AlignmentFragment> iter = compatibleFragsList.iterator();
		List<Exon> exons = new ArrayList<Exon>();
		while (iter.hasNext()) {
			Exon exon = new Exon();
			AlignmentFragment alignmentFragment = (AlignmentFragment) iter.next();
			exon.setRange(alignmentFragment.getNucleotideSeqRange());
			exon.setAlignmentFragment(alignmentFragment);
			exon.setFrame(alignmentFragment.getFrame());
			Frame sequenceFrame= VigorFunctionalUtils.getSequenceFrame(exon.getRange().getBegin()+exon.getFrame().getFrame()-1);
			exon.setSequenceFrame(sequenceFrame);
			exons.add(exon);
		}
		return exons;
	}

}
