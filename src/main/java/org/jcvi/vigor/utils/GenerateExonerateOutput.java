package org.jcvi.vigor.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jcvi.jillion.core.Sequence;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.exception.VigorException;


public class GenerateExonerateOutput {
	
		
	public static String queryExonerate(VirusGenome virusGenome, String referenceDB,
			String workspace,String proteinID,String exoneratePath) throws VigorException{

        File file = new File(workspace + File.separator + "sequence_temp.fasta");
		Path path = Paths.get(file.getAbsolutePath());
		Iterator<String> sequence = SequenceUtils.steamOf(virusGenome.getSequence(), 70).iterator();
		try (BufferedWriter writer = Files.newBufferedWriter(path)) {
			writer.write(">" + virusGenome.getId() + " "
					+ virusGenome.getDefline());
			writer.newLine();
			while(sequence.hasNext()) {
				writer.write(sequence.next());
				writer.newLine();
			};

		} catch (IOException e) {
			throw new VigorException(String.format("problem creating input file %s", file), e);
		}

		String fileName = virusGenome.getId().replaceAll("\\|", "") + ".txt";
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
				Iterator<String> pSequence = SequenceUtils.steamOf( record.getSequence(), 70).iterator();
				try (BufferedWriter writer = Files.newBufferedWriter(dbpath)) {
					writer.write(">" + record.getId() + " "
							+ record.getComment());
					writer.newLine();
					while (pSequence.hasNext()) {
						writer.write(pSequence.next());
						writer.newLine();
					}
				} catch (IOException e) {
					throw new VigorException(String.format("problem writing reference file %s for protein %s", dbFileTemp, proteinID), e);
				}
				dbPath = dbFileTemp.getAbsolutePath();
			}

			File refDBFolderFile = Paths.get(workspace, refDBFolder).toFile();
		  	if (! (refDBFolderFile.exists() || refDBFolderFile.mkdirs() )) {
		  		throw new VigorException(String.format("unable to create folder %s", refDBFolderFile.getAbsolutePath()));
			}

			List<String> exonerateCommand = Arrays.asList(exoneratePath, "--model",
					"protein2genome", "-q", dbPath, "-t",
					file.getAbsolutePath(), "--showcigar", "true");

			Process p2 = new ProcessBuilder(exonerateCommand)
					.redirectOutput(Paths.get(workspace, refDBFolder, fileName).toFile())
					.start();
			int result = p2.waitFor();
			p2.destroy();
			if (result != 0) {
				throw new VigorException(String.format("exonerate process %s returned with non-zero exit code %s",
						String.join(" ", exonerateCommand), result));
			}

		} catch (Exception e) {
			throw new VigorException(String.format("Exception running exonerate. got %s: %s", e.getClass().getSimpleName(), e.getMessage()), e);
		}

	 return Paths.get(workspace , refDBFolder , fileName).toString();
	}

}
