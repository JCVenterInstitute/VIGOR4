package org.jcvi.vigor.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Stream;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.jcvi.vigor.component.VirusGenome;


public class GenerateExonerateOutput {
	
		
	public static String queryExonerate(VirusGenome virusGenome, String referenceDB,
			String workspace,String proteinID,String exoneratePath) {

        File file = new File(workspace + File.separator + "sequence_temp.fasta");
		Path path = Paths.get(file.getAbsolutePath());
		List<String> sequence = java.util.Arrays.asList(virusGenome
				.getSequence().toString().split("(?<=\\G.{70})"));
		try (BufferedWriter writer = Files.newBufferedWriter(path)) {
			writer.write(">" + virusGenome.getId() + " "
					+ virusGenome.getDefline());
			sequence.stream().forEach(line -> {
				try {
					writer.write("\n");
					writer.write(line);
				} catch (IOException e) {
					System.out.println(e.getMessage());
				}
			});
			writer.close();

		} catch (IOException e) {
			System.out.println(e.getMessage());
		}
		String fileName = "";
		fileName = virusGenome.getId().replaceAll("\\|", "") + ".txt";
		File dbFile = new File(referenceDB);
		String refDBFolder = dbFile.getName().replaceAll("_db", "");
		String dbPath=referenceDB;
		try {

		  	if(proteinID!=null){
		  	  File dbFileTemp = new File(workspace+File.separator+"db_temp.fasta");
		  	  ProteinFastaDataStore dataStore = new ProteinFastaFileDataStoreBuilder(
						dbFile).hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
								.build();
				Stream<ProteinFastaRecord> records = dataStore.records();
				ProteinFastaRecord record = records.filter(r->r.getId().equals(proteinID)).findFirst().get();
				Path dbpath = Paths.get(dbFileTemp.getAbsolutePath());
				List<String> Psequence = java.util.Arrays.asList(record
						.getSequence().toString().split("(?<=\\G.{70})"));
				try (BufferedWriter writer = Files.newBufferedWriter(dbpath)) {
					writer.write(">" + record.getId() + " "
							+ record.getComment());
					Psequence.stream().forEach(line -> {
						try {
							writer.write("\n");
							writer.write(line);
						} catch (IOException e) {
							System.out.println(e.getMessage());
						}
					});
					writer.close();

				} catch (IOException e) {
					System.out.println(e.getMessage());
				}
				dbPath = dbFileTemp.getAbsolutePath();
			}
			String command2 = workspace + File.separator + refDBFolder;
      
			Process p1 = new ProcessBuilder("mkdir", command2).start();
			p1.waitFor();
			p1.destroy();
			Process p2 = new ProcessBuilder(exoneratePath, "--model",
					"protein2genome", "-q", dbPath, "-t",
					file.getAbsolutePath(), "--showcigar", "true")
					.redirectOutput(
							new File(workspace + File.separator + refDBFolder
									+ File.separator + fileName)).start();
			p2.waitFor();
			p2.destroy();

		} catch (Exception e) {
			e.printStackTrace();
		}

	 return workspace +File.separator + refDBFolder + File.separator + fileName;
	}

}
