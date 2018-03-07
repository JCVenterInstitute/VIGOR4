package org.jcvi.vigor;

import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.ConfigurationSource;
import org.apache.logging.log4j.core.config.Configurator;
import org.jcvi.jillion.core.datastore.DataStoreException;
import org.jcvi.vigor.service.VigorInputValidationService;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.context.annotation.Configuration;




@Configuration
@ComponentScan("org.jcvi.vigor")
public class AppConfig {
	private static final Logger LOGGER = LogManager.getLogger(AppConfig.class);

	public static void main(String... args) throws IOException, DataStoreException {
		LOGGER.debug("******** VIGOR4 APPLICATION STARTED ************");
		try (AnnotationConfigApplicationContext ctx = new AnnotationConfigApplicationContext(AppConfig.class); ) {
			ctx.getBean(VigorInputValidationService.class).processInput(args);
		}
		LOGGER.debug("******** VIGOR4 APPLICATION STOPPED ************");
	}
	
	static {
		try {
			ConfigurationSource source = new ConfigurationSource(ClassLoader.getSystemResourceAsStream("log4j2.xml"));
			Configurator.initialize(null, source);
		} catch (FileNotFoundException e) {
			LOGGER.error(e.getMessage(),e);
		} catch (IOException e) {
			LOGGER.error(e.getMessage(),e);
		}
	}

}
