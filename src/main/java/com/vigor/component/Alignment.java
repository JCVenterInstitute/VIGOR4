package com.vigor.component;

import java.util.List;
import java.util.Map;

import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class Alignment {

	private String alignmentTool_name;
	private List<AlignmentFragment> alignmentFragments;
	private Map<String,Float> alignmentScore;
	private VirusGenome virusGenome;
	private ViralProtein viralProtein;
	private AlignmentEvidence alignmentEvidence;



}
