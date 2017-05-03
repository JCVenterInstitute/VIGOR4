package com.vigor.forms;

import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import com.vigor.component.AlignmentEvidence;
import com.vigor.component.VirusGenome;

import lombok.Data;
import org.springframework.context.annotation.Scope;

@Data
public class VigorForm {

	private AlignmentEvidence alignmentEvidence;
    private HashMap<String,String> vigorParametersList;
	private VirusGenome virusGenome;

}
