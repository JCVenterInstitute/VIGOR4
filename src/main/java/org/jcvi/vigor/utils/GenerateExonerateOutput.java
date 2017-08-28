package org.jcvi.vigor.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import org.jcvi.vigor.component.VirusGenome;

public class GenerateExonerateOutput {
	
		
	public static String queryExonerate(VirusGenome virusGenome, String referenceDB,
			String workspace) {
  
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

		} catch (IOException e) {
			System.out.println(e.getMessage());
		}

		String fileName = "";
		fileName = virusGenome.getId().replaceAll("\\|", "") + ".txt";
		String refDBFolder = referenceDB.replaceAll("_db", "");
		try {
			String dbPath = VigorUtils.getVirusDatabasePath() + File.separator
					+ referenceDB;
			String command2 = workspace + File.separator + refDBFolder;

			Process p1 = new ProcessBuilder("mkdir", command2).start();
			p1.waitFor();
			p1.destroy();
			Process p2 = new ProcessBuilder("exonerate", "--model",
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
