package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.jillion.core.residue.nt.Triplet;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.jcvi.vigor.utils.VigorUtils;

public class DetermineStop implements EvaluateModel {

	@Override
	public List<Model> determine(Model model, VigorForm form) {
		List<Model> models = null;
		
		
		return models;
	}
	
	public void findInternalStops(Model model){
		List<Exon> exons = model.getExons();
		for(Exon exon: exons){
		   NucleotideSequence CDSSequence = new NucleotideSequenceBuilder(model.getAlignment().getVirusGenome().getSequence()).build();
		   
		}
	}
	
	/*public List<Model> findStop(Model model,
			String startCodonWindowParam) {
			List<Model> newModels = new ArrayList<Model>();
			long start;
			long end;
			Range startSearchRange;
			int windowSize = 50;
			boolean isSequenceMissing = false;
			if (startCodonWindowParam != null) {
				if (VigorUtils.is_Integer(startCodonWindowParam)) {
					windowSize = Integer.parseInt(startCodonWindowParam);
				}
			}
			Exon firstExon = model.getExons().get(0);
			Frame frame = VigorFunctionalUtils.getFrame(firstExon.getRange().getBegin());
			String spliceform = "";
			if (model.getAlignment().getViralProtein().getGeneAttributes()
					.getSplicing().getSpliceform() != null) {
				spliceform = model.getAlignment().getViralProtein()
						.getGeneAttributes().getSplicing().getSpliceform();
			}
			if (spliceform.equals("")) {
				long expectedStart = firstExon.getRange().getBegin()
						- ((firstExon.getAlignmentFragment().getProteinSeqRange()
								.getBegin()) * 3);
				start = expectedStart - windowSize;
				end = expectedStart + windowSize;
			} else {
				start = firstExon.getRange().getBegin() - windowSize;
				end = firstExon.getRange().getBegin() + windowSize;
			}
			if (start < 0) {
				isSequenceMissing = true;
				start = 0;
			}
			startSearchRange = Range.of(start, end);
			final long tempStart = start;
			NucleotideSequence NTSequence = model.getAlignment().getVirusGenome()
					.getSequence().toBuilder(startSearchRange).build();
			Map<Range, Float> rangeScoreMap = new HashMap<Range, Float>();
			for (Triplet triplet : startCodons) {
				Stream<Range> stream = NTSequence.findMatches(triplet.toString());
				List<Range> rangesInFrame = new ArrayList<Range>();
				List<Range> foundRanges = stream.map(
						range -> {
							range = Range.of(range.getBegin() + tempStart,
									range.getEnd() + tempStart);
							return range;
						}).collect(Collectors.toList());
				foundRanges.stream().forEach(x -> {
					if (VigorFunctionalUtils.getFrame(x.getBegin()).compareTo(frame) == 0) {
						rangesInFrame.add(x);
					}
				});
				for(Range range:rangesInFrame){
				rangeScoreMap.put(range,VigorFunctionalUtils.generateScore(firstExon.getRange().getBegin(), range.getBegin()));
				}
			}
			rangeScoreMap
					.entrySet()
					.stream()
					.sorted(Map.Entry.comparingByValue())
					.collect(
							Collectors.toMap(Map.Entry::getKey,
									Map.Entry::getValue, (e1, e2) -> e2,
									LinkedHashMap::new));
			if (!(rangeScoreMap.isEmpty())) {
				Set<Range> keys = rangeScoreMap.keySet();
				for (Range range : keys) {
					Model newModel = new Model();
					newModel = Model.deepClone(model);
					newModel.setStartCodon(range);
					if (model.getScores() != null) {
						newModel.getScores().put("startCodonScore",
								rangeScoreMap.get(range));
					} else {
						Map<String, Float> scores = new HashMap<String, Float>();
						scores.put("startCodonScore", rangeScoreMap.get(range));
						newModel.setScores(scores);
					}
					newModels.add(newModel);
				}
			}
			if (rangeScoreMap.isEmpty() && isSequenceMissing) {
				Model newModel = new Model();
				newModel = Model.deepClone(model);
				newModel.setPartial5p(true);
				newModels.add(newModel);

			} else if (rangeScoreMap.isEmpty()) {
				Model newModel = new Model();
				newModel = Model.deepClone(model);
				newModels.add(newModel);
			}

			return newModels;
		}
*/

}
