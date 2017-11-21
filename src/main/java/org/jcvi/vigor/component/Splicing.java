package org.jcvi.vigor.component;

import lombok.Data;

import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import java.util.List;



/**
 * Created by snettem on 5/8/2017.
 */
@Component
@Scope("prototype")
@Data
public class Splicing {

    private boolean isSpliced=false;
    private List<SpliceSite> nonCanonical_spliceSites; 
    private String spliceform;

    public class SpliceSite{
    	public String donor;
    	public String acceptor;
    	public SpliceSite(String donor,String acceptor){
    		this.donor=donor;
    		this.acceptor=acceptor;
    	}
    	public SpliceSite(){
    		
    	}
    	    	
    	
    }
  

}

