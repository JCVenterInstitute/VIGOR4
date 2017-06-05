package com.vigor.forms;

import com.vigor.component.AlignmentEvidence;
import com.vigor.component.Model;
import lombok.Data;

import java.util.List;
import java.util.Map;

@Data
public class VigorForm {

	private Map<String,String> vigorParametersList;
	private AlignmentEvidence alignmentEvidence;
	private String alignmentTool;
	private List<Model> models;
	private boolean debug=true;


}
