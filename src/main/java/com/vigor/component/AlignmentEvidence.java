package com.vigor.component;

import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class AlignmentEvidence {
	
	private String reference_db;
	private String matpep_db;
	
	

}
