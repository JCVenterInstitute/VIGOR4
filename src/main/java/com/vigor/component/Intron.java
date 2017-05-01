package com.vigor.component;

import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Data
public class Intron {

	private long begin;
	private long end;
	private long size;
	
}
