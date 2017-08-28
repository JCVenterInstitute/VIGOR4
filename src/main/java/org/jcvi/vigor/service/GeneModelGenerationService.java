package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.List;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.FormatVigorOutput;

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
	
	}

