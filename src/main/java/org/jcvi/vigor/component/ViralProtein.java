package org.jcvi.vigor.component;

import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import java.util.List;
import org.jcvi.jillion.core.Range;
import lombok.Data;

@Component
@Scope("prototype")
@Data
public class ViralProtein {
	
	/**
	 * 
	 */
	private ProteinSequence sequence;
    private GeneAttributes geneAttributes;
    private String proteinID;
    private String geneSymbol;
    private String defline;
    private String product;


    
}
