package org.jcvi.vigor.component;


import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class Model implements Cloneable{
	
	private List<Exon> exons;
	private Alignment alignment;
	private Map<String,Float> scores;
	private String geneSymbol;
	private List<String> status;
	private Direction direction;
    private boolean partial5p=false;
    private boolean partial3p=false;
    private Range startCodon;
    private boolean isPseudogene=false;
    
   public Model clone() throws CloneNotSupportedException {
	   Model model = (Model) super.clone();
	   model.exons = exons.stream().map(x->x.clone()).collect(Collectors.toList());
	
		
	/*
	   
	   
	   try {
	     ByteArrayOutputStream baos = new ByteArrayOutputStream();
	     ObjectOutputStream oos = new ObjectOutputStream(baos);
	     oos.writeObject(model);
	     ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
	     ObjectInputStream ois = new ObjectInputStream(bais);
	     return (Model)(ois.readObject());
	   }
	   catch (Exception e) {
	     e.printStackTrace();
	     return null;
	   }
	 }*/
return model;

}
}
