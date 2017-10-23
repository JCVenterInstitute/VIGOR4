package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.FormatVigorOutput;
import org.jcvi.vigor.utils.VigorFunctionalUtils;

@Service
public class GeneModelGenerationService {
	@Autowired
	private DetermineMissingExons determineMissingExons;
	@Autowired
	private DetermineStart determineStart;
	boolean isDebug = true;
	List<Model> partialGeneModels = new ArrayList<Model>();
	List<Model> pseudoGenes = new ArrayList<Model>();

	public void generateGeneModel(List<Model> models, VigorForm form) {
		
		evaluateGeneFeatures(models,form);
		isDebug = form.isDebug();
		/*Determine missing exons*/
		
	}
	
	public List<Model> evaluateGeneFeatures(List<Model> models,VigorForm form){
			
		List<Model> processedModels = new ArrayList<Model>();
		List<Model> modelsWithMissingExonsDetermined = new ArrayList<Model>();
	    List<Model> modelsAfterDeterminingStart = new ArrayList<Model>();
		
	    /* Determine Missing Exons */
		models.stream().forEach(model -> { 
			modelsWithMissingExonsDetermined.addAll(determineMissingExons.determine(model, form));
		});			
		if(isDebug)FormatVigorOutput.printModels(modelsWithMissingExonsDetermined,"After determining missing exons");
				
		/* Determine Start */
		modelsWithMissingExonsDetermined.stream().forEach(x -> { 
			if(!x.isPartial5p()){
			List<Model> outputModels = determineStart.determine(x, form);
			outputModels.stream().forEach(y->{
				if((y.getStartCodon()==null || y.getStartCodon().isEmpty())&& !(y.isPartial5p())) {
					pseudoGenes.add(y);
				}else if (y.getStartCodon().isEmpty() && y.isPartial5p()){
					partialGeneModels.add(y);
				}else{
					modelsAfterDeterminingStart.add(y);
				}
				
			});
			}});
		
		/*AdjustExonBoundaries*/
		
		
	
		
		return processedModels;
	}
	
	public NucleotideSequence determineCDS(Model model){
		
		model.getExons().sort(Exon.Comparators.Ascending);
		List<Exon> exons = model.getExons();
		Exon replaceExon = exons.stream().filter(thisExon-> thisExon.getReplacementString()!="").findAny().map(thisExon -> thisExon).orElse(null);
		NucleotideSequence virusGenomeSeq = new NucleotideSequenceBuilder(model.getAlignment().getVirusGenome().getSequence()).replace(replaceExon.getRange(), replaceExon.getReplacementString()).build();
		NucleotideSequenceBuilder CDSBuilder = new NucleotideSequenceBuilder("");
		for(Exon exon : exons){
			if(exon.getInsertionString()!=""){
			CDSBuilder.append(exon.getInsertionString());
			}else if(exon.getReplacementString().equals("") || exon.getReplacementString()==null) {
			CDSBuilder.append(virusGenomeSeq.toBuilder(exon.getRange()));
			}
			
		}
		NucleotideSequence cds = CDSBuilder.build();
		return cds;		
				
	}
	
	public List<Range> getInternalStops(Model model, NucleotideSequence cds){
		//dont have to translate, just find stops by calling a function
		List<Range> internalStops = new ArrayList<Range>();
		Map<Frame,List<Long>> stops = IupacTranslationTables.STANDARD.findStops(cds);
		for(Map.Entry<Frame,List<Long>> pair :stops.entrySet()){
			if(pair.getKey().equals(Frame.ONE)){
			List<Long> cdsStops = (List<Long>)pair.getValue();
			for(Long stop : cdsStops){
				Range NTStopRange = VigorFunctionalUtils.getNTRange(model.getExons(), Range.of(stop));
				internalStops.add(NTStopRange);			
			}
		}
		}
		return internalStops;
	}
	
	
	
	
	}

