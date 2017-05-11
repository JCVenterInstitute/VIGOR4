package com.vigor.component;

import org.jcvi.jillion.core.Range;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class Intron {

	private Range range;

	
}
