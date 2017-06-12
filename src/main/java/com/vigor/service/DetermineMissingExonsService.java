package com.vigor.service;

import java.util.ArrayList;
import java.util.List;

import org.jcvi.jillion.core.Range;

import com.vigor.component.Exon;
import com.vigor.component.Model;
import com.vigor.component.ViralProtein;


public class DetermineMissingExonsService implements EvaluateModel{

	@Override
	public Model determine(Model model) {
		findMissingExons(model);
		return model;
	}
	public List<Range> findMissingExons(Model model){
		List<Range> missingExonRanges = new ArrayList<Range>();
		ViralProtein viralProtein = model.getAlignment().getViralProtein();
		List<Exon> referenceExons = viralProtein.getGeneStructure().getExons();
		List<Exon> modelExons = model.getExons();
		for(int i=0;i<referenceExons.size();i++){
			Exon refExon = referenceExons.get(i);
			boolean found=false;
			for(int j=0;j<modelExons.size();j++){
				Exon modelExon = modelExons.get(j);
				if(refExon.getRange().isSubRangeOf(modelExon.getRange())){
					modelExons.get(j).setRange(refExon.getRange());
					found=true;
				}
				else if(refExon.getRange().intersects(modelExon.getRange())){
					found=true;
				}
			}
			if(!found){
				if(i>0 && i<referenceExons.size()-1){
				Range range = Range.of(referenceExons.get(i+1).getRange().getBegin(),referenceExons.get(i-1).getRange().getEnd());
				missingExonRanges.add(range);
				}
				
			}
		}
		
		return missingExonRanges;
	}
	
	
	
}
