package org.jcvi.vigor.forms;


import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.AlignmentTool;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.utils.LoadDefaultParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorUtils;
import lombok.Data;


@Data
public class VigorForm {

	private VigorConfiguration configuration;
	private AlignmentEvidence alignmentEvidence;
	private AlignmentTool alignmentTool;
	private boolean debug;
	private String AlignmentOutputTempFile;
	private String tempDirectoryPath;

	public VigorForm() throws VigorException
	{
		configuration = LoadDefaultParameters
				.loadVigorConfiguration("defaults",Thread.currentThread().getContextClassLoader().getResource(VigorUtils.getDefaultConfigurationPath()));
	}

	public VigorForm(VigorConfiguration configuration) {
		this.configuration = new VigorConfiguration(configuration);
	}
	
}
