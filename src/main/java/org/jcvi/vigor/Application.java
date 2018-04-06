package org.jcvi.vigor;

import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.context.annotation.Configuration;

@Configuration
@ComponentScan("org.jcvi.vigor")
public class Application {

	public static void main(String... args)  {
    	try (AnnotationConfigApplicationContext ctx = new AnnotationConfigApplicationContext(Application.class); ) {
			ctx.getBean(Vigor.class).run(args);
		}
	}
}
