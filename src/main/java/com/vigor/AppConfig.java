package com.vigor;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.ConfigurationSource;
import org.apache.logging.log4j.core.config.Configurator;
import org.jcvi.jillion.core.datastore.DataStoreException;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.core.io.IOUtil;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.jillion.core.util.iter.StreamingIterator;
import org.jcvi.jillion.fasta.nt.*;
import org.springframework.context.ApplicationContext;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;
import org.springframework.context.annotation.ComponentScan;
import org.springframework.context.annotation.Configuration;
import com.vigor.service.VigorInputValidationService;



@Configuration
@ComponentScan("com.vigor")
public class AppConfig {
	private static final Logger LOGGER = LogManager.getLogger(AppConfig.class);

	public static void main(String... args) throws IOException, DataStoreException {
		LOGGER.debug("******** VIGOR4 APPLICATION STARTED ************");
		ApplicationContext ctx = new AnnotationConfigApplicationContext(AppConfig.class);
		ctx.getBean(VigorInputValidationService.class).processInput(args);
		((AnnotationConfigApplicationContext)ctx).close();
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
