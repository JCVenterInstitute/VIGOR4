package com.vigor.component;

import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Data
public class Exon {

	private long begin;
	private long end;
	private int strand;
	private long size;
	
	
	
}
