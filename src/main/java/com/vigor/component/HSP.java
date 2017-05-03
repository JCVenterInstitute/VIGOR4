package com.vigor.component;

import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class HSP {

	private double score;
	private String evalue;
	private int strand;
	private float positve;
	private float identity;
	
	public HSP(double score,String evalue, int strand,float positive,float identity)
	{
		this.score=score;
		this.evalue=evalue;
		this.strand=strand;
		this.positve=positive;
		this.identity=identity;
	}
	public HSP()
	{

	}
	
	
}
