package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.Triplet;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.*;

@Service
public class DetermineStart implements EvaluateModel {

	private static final Logger LOGGER = LogManager
			.getLogger(DetermineStart.class);

	@Override
	public List<Model> determine(Model model, VigorForm form) {
		List<Model> models = null;
		try {
			List<Triplet> startCodons = LoadStartCodons(model.getAlignment()
					.getViralProtein().getGeneAttributes()
					.getStartTranslationException().getAlternateStartCodons(),
					form.getVigorParametersList());
			String startCodonWindowParam = form.getVigorParametersList().get(
					"start_codon_search_window");
			models = findStart(startCodons, model, startCodonWindowParam);

			for (Model mymodel : models) {
				System.out.println("Finding start for below model");
				mymodel.getExons().stream().forEach(System.out::println);
				System.out.println(mymodel.getGeneSymbol() + "    "
						+ mymodel.getStartCodon());
			}

		}

		catch (Exception e) {
			LOGGER.error(e.getMessage(), e);
		}

		return models;
	}

	/**
	 * 
	 * @param alternateStartCodons
	 *            defined in defline of reference protein
	 * @param vigorParameters
	 *            : Default start codons defined in vigor.ini file
	 * @return List of all the possible start codons
	 */
	public List<Triplet> LoadStartCodons(List<String> alternateStartCodons,
			Map<String, String> vigorParameters) {
		String startCodonsParam;
		List<String> startCodonStrings = new ArrayList<String>();
		if (vigorParameters.containsKey("StartCodons")) {
			startCodonsParam = vigorParameters.get("StartCodons");
			startCodonStrings = Arrays.asList(StringUtils.normalizeSpace(
					startCodonsParam).split(","));
		} else {

			startCodonStrings.add("ATG");
		}
		if (alternateStartCodons != null) {
			startCodonStrings.addAll(alternateStartCodons);
		}
		List<Triplet> startCodons = new ArrayList<Triplet>();
		for (String startCodonString : startCodonStrings) {
			if (startCodonString.length() == 3) {
				Triplet triplet = Triplet.create(startCodonString.charAt(0),
						startCodonString.charAt(1), startCodonString.charAt(2));
				startCodons.add(triplet);
			}

		}
		return startCodons;
	}

	public List<Model> findStart(List<Triplet> startCodons, Model model,
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
		Frame frame = VigorFunctionalUtils.getFrame(firstExon.getRange());
		String spliceform = "";
		if (model.getAlignment().getViralProtein().getGeneAttributes()
				.getSplicing().getSpliceform() != null) {
			spliceform = model.getAlignment().getViralProtein()
					.getGeneAttributes().getSplicing().getSpliceform();
		}
		if (spliceform.equals("")) {
			long expectedStart = firstExon.getRange().getBegin()
					- ((firstExon.getAlignmentFragment().getProteinSeqRange()
							.getBegin() - 1) * 3);
			start = expectedStart - windowSize;
			end = expectedStart - windowSize;
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
				if (VigorFunctionalUtils.getFrame(x).compareTo(frame) == 0) {
					rangesInFrame.add(x);
				}
			});
			rangeScoreMap.putAll(VigorFunctionalUtils.generateScore(firstExon
					.getRange().getBegin(), rangesInFrame));
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

}
