package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
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
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.springframework.stereotype.Service;

import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.VigorUtils;

@Service
public class DetermineStart implements EvaluateModel {

	private static final Logger LOGGER = LogManager.getLogger(DetermineStart.class);

	@Override
	public Model determine(Model model, VigorForm form) {
		try {
			List<String> startCodons = LoadStartCodons(model.getAlignment().getViralProtein().getGeneAttributes()
					.getStartTranslationException().getAlternateStartCodons(), form.getVigorParametersList());
			String startCodonWindowParam = form.getVigorParametersList().get("start_codon_search_window");
            
			findStart(startCodons, model, startCodonWindowParam);

		}
        catch(CloneNotSupportedException e){
        	LOGGER.error(e.getMessage(), e);
        }
		catch (Exception e) {
			LOGGER.error(e.getMessage(), e);
		}

		return model;
	}

	/**
	 * 
	 * @param alternateStartCodons
	 *            defined in defline of reference protein
	 * @param vigorParameters
	 *            : Default star codons defined in vigor.ini file
	 * @return List of all the possible start codons
	 */
	public List<String> LoadStartCodons(List<String> alternateStartCodons, Map<String, String> vigorParameters) {
		String startCodonsParam;
		List<String> startCodons = new ArrayList<String>();
		if (vigorParameters.containsKey("StartCodons")) {
			startCodonsParam = vigorParameters.get("StartCodons");
			startCodons = Arrays.asList(StringUtils.normalizeSpace(startCodonsParam).split(","));
		} else {
			startCodons.add("ATG");
		}
		startCodons.addAll(alternateStartCodons);
		return startCodons;
	}

	public List<Model> findStart(List<String> startCodons, Model model, String startCodonWindowParam) throws CloneNotSupportedException {

		List<Model> newModels = new ArrayList<Model>();
		int windowSize = 5;
		boolean missingSequence = false;
		if (startCodonWindowParam != null) {
			if (VigorUtils.is_Integer(startCodonWindowParam)) {
				windowSize = Integer.parseInt(startCodonWindowParam);
			}
		}
		Exon firstExon = model.getExons().get(0);
		
		long start = firstExon.getRange().getBegin() - windowSize;
		long centroid = firstExon.getRange().getBegin()-start;
		if (start < 0) {
			missingSequence = true;
			start = 0;
		}

		long end = firstExon.getRange().getBegin() + windowSize;
		Range startSearchRange = Range.of(start, end);
		NucleotideSequence NTSequence = model.getAlignment().getVirusGenome().getSequence().toBuilder(startSearchRange)
				.build();
		long difference = 0;
	
		Map<Range,Long> rangeScoreMap = new HashMap<Range,Long>(); 
		for (String codon : startCodons) {
			Stream<Range> stream = NTSequence.findMatches(codon);
			List<Range> foundRanges = stream.collect(Collectors.toList());
			for (Range range : foundRanges) {
						
				if (range.getBegin() - centroid < 0) {
					difference = centroid - range.getBegin();
					rangeScoreMap.put(range, difference);
				} else{
					difference = range.getBegin() - centroid;}
				  rangeScoreMap.put(range, difference);				
			}				
		}		
		rangeScoreMap.entrySet().stream().sorted(Map.Entry.comparingByValue()).collect(Collectors.toMap(Map.Entry::getKey,Map.Entry::getValue,(e1, e2) -> e2,LinkedHashMap::new));
		Set<Range> keys = rangeScoreMap.keySet();
		float score = 10;
		for(Range range : keys ){
			Model newModel = new Model();
			newModel = (Model)model.clone();
			newModel.setStartCodon(range);
			Map<String,Float> scores = new HashMap<String,Float>();
			scores.putAll(newModel.getScores());
			scores.put("startCodonOccurence", score);
			newModel.setScores(scores);
			newModels.add(newModel);
		}
			
		
		
		return newModels;
	}

}
