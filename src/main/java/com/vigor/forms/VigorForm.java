package com.vigor.forms;

import com.vigor.component.AlignmentEvidence;
import com.vigor.utils.LoadDefaultParameters;
import com.vigor.utils.VigorUtils;

import lombok.Data;
import java.util.Map;

@Data
public class VigorForm {

	private Map<String,String> vigorParametersList;
	private AlignmentEvidence alignmentEvidence;
	private String alignmentTool;
	private boolean debug=true;

public VigorForm()
{
	vigorParametersList = LoadDefaultParameters
		.loadVigorParameters(VigorUtils.getVigorParametersPath());
}
	
}
