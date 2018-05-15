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
import org.jcvi.vigor.service.exception.ServiceException;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.*;

@Service
public class DetermineStart implements DetermineGeneFeatures {

	private static final Logger LOGGER = LogManager
			.getLogger(DetermineStart.class);

	@Override
	public List<Model> determine(Model model, VigorForm form) throws ServiceException {
		List<Triplet> startCodons = loadStartCodons(model.getAlignment()
														 .getViralProtein().getGeneAttributes()
														 .getStartTranslationException().getAlternateStartCodons(),
				form.getConfiguration());
		String startCodonWindowParam = form.getConfiguration().get(
				ConfigurationParameters.StartCodonSearchWindow);
		try {
			List<Model> models = findStart(startCodons, model, startCodonWindowParam);
			/*LOGGER.debug("Models after determining start: {}", () ->
				models.stream().map(String::valueOf).collect(Collectors.joining("\n")));*/
			return models;
		} catch (CloneNotSupportedException e) {
			LOGGER.error("for model {} problem finding start using codons {} and search window {}",
					model,
					startCodons.stream().map(String::valueOf).collect(Collectors.joining(",")),
					startCodonWindowParam);
			throw new ServiceException(e);
		}

	}

	/**
	 * 
	 * @param alternateStartCodons
	 *            defined in defline of reference protein
	 * @param vigorParameters
	 *            : Default start codons defined in vigor.ini file
	 * @return List of all the possible start codons
	 */
	public List<Triplet> loadStartCodons(List<String> alternateStartCodons,
			VigorConfiguration vigorParameters) {
		String startCodonsParam;
		List<String> startCodonStrings = new ArrayList<String>();
		if (vigorParameters.containsKey(ConfigurationParameters.StartCodons)) {
			startCodonsParam = vigorParameters.get(ConfigurationParameters.StartCodons);
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
		String startCodonWindowParam) throws CloneNotSupportedException {
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
		Frame firstExonFrame = VigorFunctionalUtils.getSequenceFrame(firstExon.getRange().getBegin());
		long expectedStart = firstExon.getRange().getBegin()
					- ((firstExon.getAlignmentFragment().getProteinSeqRange()
							.getBegin()) * 3);
		expectedStart = expectedStart >= 0 ? expectedStart : 0;

		start = expectedStart - windowSize;
		if (start < 0) {
			isSequenceMissing = true;
			start = 0;
		}
		end = expectedStart+windowSize;
		if(end>firstExon.getRange().getEnd()){
		    end = firstExon.getRange().getEnd();
        }
		startSearchRange = Range.of(start, end);
		//find any internal stops
        Map<Frame,List<Long>> internalStops = model.getAlignment().getVirusGenome().getInternalStops();
        List<Long> stopsInFrame=new ArrayList<Long>();

        if(internalStops!=null) {
            for (Frame frame : internalStops.keySet()) {
                if (frame.equals(firstExonFrame)) {
                    List<Long> tempList = internalStops.get(firstExonFrame);
                    for(Long stop : tempList){
                        Range tempRange = Range.of(stop);
                        if(tempRange.isSubRangeOf(startSearchRange)){
                            stopsInFrame.add(stop);
                        }
                    }
                }
            }
        }


		final long tempStart = start;
		NucleotideSequence NTSequence = model.getAlignment().getVirusGenome()
				.getSequence().toBuilder(startSearchRange).build();
		Map<Range, Double> rangeScoreMap = new HashMap<Range, Double>();
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
				if (VigorFunctionalUtils.getSequenceFrame(x.getBegin()).compareTo(firstExonFrame) == 0) {
					rangesInFrame.add(x);
				}
			});
			for(Range range:rangesInFrame){
			    boolean isValid = true;
			    for(Long stop: stopsInFrame){
			        if(range.endsBefore(Range.of(stop,stop+2))){
			            isValid=false;
                    }
                }
                if(isValid) {
                    rangeScoreMap.put(range, VigorFunctionalUtils.generateScore(firstExon.getRange().getBegin(), range.getBegin()));
                }
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

				//Model newModel = new Model();
				Model newModel = model.clone();
				Exon fExon = newModel.getExons().get(0);
				fExon.setRange(Range.of(range.getBegin(),fExon.getRange().getEnd()));
				newModel.getScores().put("startCodonScore",
                           rangeScoreMap.get(range));
				   newModels.add(newModel);
			}
		}
		//set 5' partial and extend start of the first exon to the beginning of the sequence
		if (rangeScoreMap.isEmpty() && isSequenceMissing) {
			//Model newModel = new Model();
			Model newModel = model.clone();
			newModel.setPartial5p(true);
			Exon fExon = newModel.getExons().get(0);
			Range fExonRange = fExon.getRange();
			long bases = fExonRange.getBegin()-0;
			int frameshift = (int)bases % 3;
			if(frameshift>0) fExon.setFrame(fExon.getFrame().shift(frameshift));
			fExon.setRange(Range.of(0,fExonRange.getEnd()));
			newModels.add(newModel);
			//System.out.println("Sequence is missin. No Start found. Partial gene "+newModel.getAlignment().getViralProtein().getProteinID());

		} else if (rangeScoreMap.isEmpty()) {
			//Model newModel = new Model();
			Model newModel = model.clone();
			newModel.setPseudogene(true);
			newModels.add(newModel);	
			//System.out.println("Pseudogene. No Start found. "+newModel.getAlignment().getViralProtein().getProteinID());
		}

		return newModels;
	}

 

	
	
}



