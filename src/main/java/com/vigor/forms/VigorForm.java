package com.vigor.forms;

import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import com.vigor.component.AlignmentEvidence;
import com.vigor.component.VirusGenome;

import lombok.Data;

@Data
public class VigorForm {

	private List<VirusGenome> inputGenomeSequencesList;
	private AlignmentEvidence alignmentEvidence;
    private HashMap<String,String> vigorParametersList;
	

}
