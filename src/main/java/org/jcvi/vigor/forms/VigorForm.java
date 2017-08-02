package org.jcvi.vigor.forms;

<<<<<<< HEAD:src/main/java/com/vigor/forms/VigorForm.java
import com.vigor.component.AlignmentEvidence;
import com.vigor.utils.LoadDefaultParameters;
import com.vigor.utils.VigorUtils;

=======
import org.jcvi.vigor.component.AlignmentEvidence;
>>>>>>> 4aafe0f0632a193b0d31d36b3f8efd0b9439888a:src/main/java/org/jcvi/vigor/forms/VigorForm.java
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
