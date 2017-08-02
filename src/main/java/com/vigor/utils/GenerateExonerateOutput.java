package com.vigor.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;

import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;

import com.jcraft.jsch.Channel;
import com.jcraft.jsch.ChannelExec;
import com.jcraft.jsch.JSch;
import com.jcraft.jsch.Session;
import com.vigor.component.VirusGenome;

public class GenerateExonerateOutput {
	
		
	public File queryExonerate(VirusGenome virusGenome, String referenceDB) {

		File file = new File(VigorUtils.getVigorWorkSpace() + "/sequence_temp.fasta");
		Path path = Paths.get(file.getAbsolutePath());
		System.out.println("Length is" + virusGenome.getSequence().getLength());
		List<String> sequence = java.util.Arrays.asList(virusGenome.getSequence().toString().split("(?<=\\G.{70})"));
		try (BufferedWriter writer = Files.newBufferedWriter(path)) {
			writer.write(">" + virusGenome.getId() + " " + virusGenome.getDefline());
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

		// connect to lserver1 and run exonerate
		String fileName = virusGenome.getId().replaceAll("\\|", "");

		try {
			Runtime.getRuntime().exec(
					"pscp -pw \"Ilovemyself@3\" -scp C:/git/VIGOR4/VigorWorkSpace/sequence_temp.fasta lserver1:/home/snettem/vigor4_workspace/");

			String dbPath = "/home/snettem/vigor4_workspace/data3/" + referenceDB;
			String host = "lserver1";
			String user = "snettem";
			String password = "Ilovemyself@3";
			String command1 = "chmod 777 /home/snettem/vigor4_workspace/sequence_temp.fasta && exonerate --model protein2genome -q "
					+ dbPath
					+ " -t /home/snettem/vigor4_workspace/sequence_temp.fasta --showcigar true > /home/snettem/vigor4_workspace/"+fileName+".txt && chmod 777 /home/snettem/vigor4_workspace/exonerate.txt";
           
			java.util.Properties config = new java.util.Properties();
			config.put("StrictHostKeyChecking", "no");
			JSch jsch = new JSch();
			Session session = jsch.getSession(user, host, 22);
			session.setPassword(password);
			session.setConfig(config);
			session.connect();
			System.out.println("Connected to unix server");

			Channel channel = session.openChannel("exec");
			((ChannelExec) channel).setCommand(command1);
			channel.setInputStream(null);
			((ChannelExec) channel).setErrStream(System.err);

			InputStream in = channel.getInputStream();
			channel.connect();
			byte[] tmp = new byte[1024];
			while (true) {
				while (in.available() > 0) {
					int i = in.read(tmp, 0, 1024);
					if (i < 0)
						break;
					System.out.print(new String(tmp, 0, i));
				}
				if (channel.isClosed()) {
					System.out.println("exit-status: " + channel.getExitStatus());
					break;
				}
				try {
					Thread.sleep(1000);
				} catch (Exception ee) {
				}
			}
			channel.disconnect();
			session.disconnect();
			
					Runtime.getRuntime().exec(
					"pscp -pw \"Ilovemyself@3\" -scp lserver1:/home/snettem/vigor4_workspace/"+fileName+".txt C:/git/VIGOR4/src/test/resources/VigorRegressionTestOutput/flu/"+fileName+".txt");
			
			System.out.println("DONE");
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}

		// *********************************
			
		File outputFile = new File(VigorUtils.getVigorWorkSpace() + "/exonerate.txt");
		return outputFile;
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
    GenerateExonerateOutput generateExonerateOutput = new GenerateExonerateOutput();
		NucleotideFastaDataStore dataStore;
		try {
			dataStore = new NucleotideFastaFileDataStoreBuilder(
					new File("VigorWorkSpace/sequence.fasta")).hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
							.build();
		 
		Stream<NucleotideFastaRecord> records = dataStore.records();
		VirusGenome virusGenome;
		Iterator<NucleotideFastaRecord> i = records.iterator();
		while (i.hasNext()) {
			NucleotideFastaRecord record = i.next();
			 virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(),
					false, false);
			 File file = generateExonerateOutput.queryExonerate(virusGenome,"flua_db");
			 			          
		}
		}
		
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}

}
