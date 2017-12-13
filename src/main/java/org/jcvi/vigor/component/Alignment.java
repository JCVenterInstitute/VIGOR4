package org.jcvi.vigor.component;
import java.io.Serializable;
import java.util.List;
import java.util.Map;

import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
@SuppressWarnings("serial")
public class Alignment implements Serializable {

	private transient String alignmentTool_name;
	private List<AlignmentFragment> alignmentFragments;
	private Map<String,Double> alignmentScore;
	private transient VirusGenome virusGenome;
	private transient ViralProtein viralProtein;
	private transient AlignmentEvidence alignmentEvidence;

}
