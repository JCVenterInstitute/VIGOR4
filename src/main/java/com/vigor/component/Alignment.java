package com.vigor.component;

import java.util.ArrayList;
import java.util.Map;

import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class Alignment {

	private String alignmentTool_name;
	private ArrayList<AlignmentFragment> alignmentFragments;
	private Map<String,Float> alignemntScore;
	private VirusGenome virusGenome;
	private ViralProtein viralProtein;
	private AlignmentEvidence alignmentEvidence;


}
