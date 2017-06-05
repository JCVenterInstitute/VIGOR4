package org.jcvi.vigor.component;

import java.util.List;
import java.util.Map;


import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class Model implements Cloneable{
	
	//These exons are the Hsps converted to exons
	private List<Exon> exons;
	private Alignment alignment;
	private Map<String,Float> scores;
	private String geneSymbol;
	private List<String> status;
	private Direction direction;
    private boolean partial5p=false;
    private boolean partial3p=false;
    private Range startCodon;
    
   public Object clone() throws CloneNotSupportedException{
	   return super.clone();
   }


}
