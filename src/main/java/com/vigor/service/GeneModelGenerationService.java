package com.vigor.service;

import java.util.List;

import org.springframework.beans.factory.annotation.Autowired;

import com.vigor.component.Model;
import com.vigor.forms.VigorForm;

public class GeneModelGenerationService {
	
	@Autowired 
	private DetermineMissingExonsService determineMissingExons;
	

	public void generateGeneModel(List<Model> models,VigorForm form){
						
		for(Model model : models){
	    determineMissingExons.determine(model);
		
		
		}
	}
	
	
	
}
