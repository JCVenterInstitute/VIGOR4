package com.vigor.service;


import java.io.IOException;
import java.util.ArrayList;

import java.util.List;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import com.vigor.component.Alignment;

import com.vigor.component.Model;



/**
 * Created by snettem on 5/19/2017.
 */


public class ModelGenerationServiceTest {
	
	private List<Alignment> alignments;

	@Before
	public void getAlignments() throws IOException{
		AlignmentGenerationServiceTest alignmentGenerationServiceTest = new AlignmentGenerationServiceTest();
		alignmentGenerationServiceTest.generateAlignmentTest();
		alignments = alignmentGenerationServiceTest.alignments;
	}
  
 
	
	
	@Test
	public void alignmentToModelsTest() throws IOException{
		
		
		List<Model> models = new ArrayList<Model>();
		ModelGenerationService modelGenerationService = new ModelGenerationService();
		for(int i=0;i<alignments.size();i++){
			models.addAll(modelGenerationService.alignmentToModels(alignments.get(i),"Exonerate"));
			
	    }
		Assert.assertEquals(7, models.size());
				  	
	}
	
	

	
	
	
	
	
    }