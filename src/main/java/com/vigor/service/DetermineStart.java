package com.vigor.service;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.springframework.stereotype.Service;

import com.vigor.component.Exon;
import com.vigor.component.Model;
import com.vigor.forms.VigorForm;
import com.vigor.utils.VigorUtils;

@Service
public class DetermineStart implements EvaluateModel {

	private static final Logger LOGGER = LogManager.getLogger(DetermineStart.class);

	@Override
	public Model determine(Model model, VigorForm form) {
		try{
		List<String> startCodons = LoadStartCodons(model.getAlignment().getViralProtein().getGeneAttributes()
				.getStartTranslationException().getAlternateStartCodons(), form.getVigorParametersList());
		}
		
		catch(Exception e){
			LOGGER.error(e.getMessage(),e);
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

	public List<Model> findStart(List<String> startCodons, Model model, Map<String,String> vigorParameters){
		List<Model> newModels = new ArrayList<Model>();
		int windowSize=5;
	    if(vigorParameters.containsKey("startCodonWindowSize")){
	      	if(VigorUtils.is_Integer(vigorParameters.get("startCodonWindowSize"))){
	    		windowSize = Integer.parseInt(vigorParameters.get("startCodonWindowSize"));
	    	}
	    }
		Exon firstExon = model.getExons().get(0);
		long start = firstExon.getRange().getBegin()-windowSize;
		if(start<0){
			start=0;
		}
		long end = firstExon.getRange().getBegin()+windowSize;
		Range range = Range.of(start,end);
		
		
		
		
		
		
		return newModels;
	}
	
		

}
