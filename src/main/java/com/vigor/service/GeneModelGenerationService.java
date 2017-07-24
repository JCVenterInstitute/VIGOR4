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

	public void generateGeneModel(List<Model> models, VigorForm form) {

		// Determine missing exons
	    List<Model> modelsAfterMissingExonsDetermined = new ArrayList<Model>();
		modelsAfterMissingExonsDetermined = models.stream().map(model -> determineMissingExons.determine(model,form)).collect(Collectors.toList());

		// Determine start
		
		
		
		
		
		
		

	}

}
