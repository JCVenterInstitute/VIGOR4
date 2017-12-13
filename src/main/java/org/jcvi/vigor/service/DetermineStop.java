package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.jcvi.vigor.utils.VigorUtils;
import org.springframework.stereotype.Service;

@Service
public class DetermineStop implements DetermineGeneFeatures {

	private static final Logger LOGGER = LogManager
			.getLogger(DetermineStop.class);
	long stopCodonWindow = 50;

	@Override
	public List<Model> determine(Model model, VigorForm form) {
		List<Model> models = null;
		String stopCodonWindowParam = form.getVigorParametersList().get("stop_codon_search_window");
		if (VigorUtils.is_Integer(stopCodonWindowParam)) {
			stopCodonWindow = Integer.parseInt(stopCodonWindowParam);
		}
		try{
		models = findStop(model);
		}
		catch(Exception e){
			LOGGER.error(e.getMessage(),e);
		}
		return models;
	}
	
	public List<Model> findStop(Model model)
		 throws CloneNotSupportedException {
		List<Model> newModels = new ArrayList<Model>();
		model.getExons().sort(Exon.Comparators.Ascending);
		List<Exon> exons = model.getExons();
		long start;
		long end;
		Range stopSearchRange;
		NucleotideSequence seq = model.getAlignment().getVirusGenome().getSequence();
		long proteinSeqLength = model.getAlignment().getViralProtein().getSequence().getLength();
		boolean isSequenceMissing = false;
		Exon lastExon = exons.get(exons.size()-1);
		Range lastExonAARange = lastExon.getAlignmentFragment().getProteinSeqRange();
		long expectedStart=0;
		if(lastExonAARange.getEnd()< proteinSeqLength-1){
			long difference = proteinSeqLength-1-lastExonAARange.getEnd();
			expectedStart = lastExon.getRange().getEnd()+difference*3;
		}else{
			expectedStart = lastExon.getRange().getEnd();
		}
		Frame lastExonFrame = VigorFunctionalUtils.getSequenceFrame(lastExon.getRange().getBegin()+(lastExon.getFrame().getFrame()-1));
			start = expectedStart - stopCodonWindow;
			end = expectedStart + stopCodonWindow;
		if (end > seq.getLength()-1 || start > seq.getLength()-1 ) {
			isSequenceMissing = true;
			end = seq.getLength()-1;
			start=end-stopCodonWindow;
		}	
		
		stopSearchRange = Range.of(start, end);	
		NucleotideSequence NTSequence = model.getAlignment().getVirusGenome()
				.getSequence().toBuilder(stopSearchRange).build();
		Map<Range, Double> rangeScoreMap = new HashMap<Range, Double>();
		Map<Frame,List<Long>> stops = IupacTranslationTables.STANDARD.findStops(NTSequence);
		for(Frame frame : stops.keySet()){
			List<Long> temp = new ArrayList<Long>();
			for(Long x : stops.get(frame)){
				x=x+start;
				temp.add(x);
			}
		stops.put(frame, temp);
		}
		stops = VigorFunctionalUtils.frameToSequenceFrame(stops);
		List<Long> stopsInFrame=null;
		for(Frame frame: stops.keySet()){
			if(frame.equals(lastExonFrame)){
				stopsInFrame = stops.get(frame);
			}
		}
		if(stopsInFrame!=null && stopsInFrame.size()>0){
		for(Long stop:stopsInFrame){
			rangeScoreMap.put(Range.of(stop,stop+2),VigorFunctionalUtils.generateScore(lastExon.getRange().getEnd(),stop));
		}
		rangeScoreMap
				.entrySet()
				.stream()
				.sorted(Map.Entry.comparingByValue())
				.collect(
						Collectors.toMap(Map.Entry::getKey,
								Map.Entry::getValue, (e1, e2) -> e2,
								LinkedHashMap::new));
		if (!(rangeScoreMap.isEmpty() && rangeScoreMap.size()>0)) {
			Set<Range> keys = rangeScoreMap.keySet();
			for (Range range : keys) {
				Model newModel = new Model();
				newModel = model.clone();
				List<Exon> newExons = newModel.getExons();
				Exon lExon = newExons.get(newExons.size()-1);
				lExon.setRange(Range.of(lExon.getRange().getBegin(),range.getBegin()-1));
				if (model.getScores() != null) {
					newModel.getScores().put("stopCodonScore",
							rangeScoreMap.get(range));
				} else {
					Map<String, Double> scores = new HashMap<String, Double>();
					scores.put("stopCodonScore", rangeScoreMap.get(range));
					newModel.setScores(scores);
				}
				
				newModels.add(newModel);
			}
		}
		}
		if (rangeScoreMap.isEmpty() && isSequenceMissing) {
			Model newModel = new Model();
			newModel = model.clone();
			newModel.setPartial3p(true);
			newModels.add(newModel);
			System.out.println("Sequence is missin. No stop found. Partial gene "+newModel.getAlignment().getViralProtein().getProteinID());

		} else if (rangeScoreMap.isEmpty()) {
			Model newModel = new Model();
			newModel = model.clone();
			newModel.setPseudogene(true);
			newModels.add(newModel);
			System.out.println("Pseudogene. No stop found. "+newModel.getAlignment().getViralProtein().getProteinID());
		}

		return newModels;
		}


}
