package org.jcvi.vigor.forms;


import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.utils.LoadDefaultParameters;
import org.jcvi.vigor.utils.VigorUtils;
import org.jcvi.vigor.component.AlignmentEvidence;
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
