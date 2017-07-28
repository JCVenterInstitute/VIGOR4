package com.vigor.service;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import com.vigor.component.Model;
import com.vigor.forms.VigorForm;
import com.vigor.utils.FormatVigorOutput;

@Service
public class GeneModelGenerationService {

	@Autowired
	private DetermineMissingExonsService determineMissingExons;
	boolean isDebug = false;
	

	public void generateGeneModel(List<Model> models, VigorForm form) {

		List<Model> partialGeneModels = new ArrayList<Model>();
		isDebug = form.isDebug();
		
		/*Determine missing exons*/
		List<Model> modelsAfterMissingExonsDetermined = new ArrayList<Model>();
		for(Model model : models){
			Model outputModel = determineMissingExons.determine(model, form);
			if(outputModel.getStatus().contains("partial genes")){
				partialGeneModels.add(outputModel);
			}
			else{
				modelsAfterMissingExonsDetermined.add(outputModel);
			}
		}
		if (isDebug) {
			System.out.println("**********After determining the missing exons************");
			FormatVigorOutput.printModels2(modelsAfterMissingExonsDetermined);
		}
		if (isDebug) {
			System.out.println("**********Partial Gene Models************");
			FormatVigorOutput.printModels2(partialGeneModels);
		}

		
		/*Determine start*/
		
		
		
		
		

	}

}
