package com.vigor.service;

import com.vigor.forms.VigorForm;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.stereotype.Service;


@Service
public class ReferenceDBGenerationService {

	private static final Logger LOGGER = LogManager.getLogger(ReferenceDBGenerationService.class);

/*	public void blastVirusGenome(VirusGenome virusGenome, VigorForm form) {
		String vigorspace = VigorUtils.getVigorWorkSpace();
		File file = new File(vigorspace + "/chooseref.fasta");
		File file1 = new File(vigorspace + "/chooseref.xml");
		File file2 = new File(vigorspace + "/chooseref.log");
		try {
			file.createNewFile();
			file1.createNewFile();
			file2.createNewFile();
		} catch (IOException e) {
			VigorException.printExceptionMessage(e.getMessage());
		}
		String inputFilePath = file.getAbsolutePath();
		String outputXMLPath = file1.getAbsolutePath();
		String logFilePath = file2.getAbsolutePath();
		BufferedWriter out = null;
		try {
			FileWriter fstream = new FileWriter(new File(inputFilePath));
			int seqLength = virusGenome.getSeqLength();
			int noOfLines = (int) Math.ceil(seqLength / 60);
			System.out.println("Number of lines" + noOfLines);
			out = new BufferedWriter(fstream);
			out.write(virusGenome.getSeqHeader().toString());
			out.newLine();
			int offset = 0;
			System.out.println(virusGenome.getNTSequence().toString().length());
			System.out.println(virusGenome.getNTSequence().toString());
			for (int i = 0; i <= noOfLines; i++) {
				if (offset + 60 < seqLength) {
					out.write(virusGenome.getNTSequence().toString(), offset, 60);
					offset = offset + 60;
				} else {
					out.write(virusGenome.getNTSequence().toString(), offset, seqLength - offset);
				}

				out.newLine();
			}
			out.flush();

			String blastFilePath = resourceLoader.getResource("classpath:" + VigorUtils.getBlastFilePath()).getFile()
					.getAbsolutePath();
			String reference_dbPath = resourceLoader.getResource("classpath:" + VigorUtils.getVirusDatabasePath()
					+ File.separator + form.getAlignmentEvidence().getReference_db()).getFile().getAbsolutePath();
			System.out.println(VigorUtils.getBlastCommand(blastFilePath, inputFilePath, reference_dbPath, outputXMLPath,
					logFilePath));
			String shellCommand = VigorUtils.getBlastCommand(blastFilePath, inputFilePath, reference_dbPath,
					outputXMLPath, logFilePath);
			runblastCommand(shellCommand, outputXMLPath);
			System.out.println("Finally");

		} catch (IOException | IllegalArgumentException e) {
			VigorException.printExceptionMessage(e.getMessage());
			LOGGER.error(e.getMessage(), e);
		} catch (Exception e) {
			VigorException.printExceptionMessage(e.getMessage());
			LOGGER.error(e.getMessage(), e);
		} finally {
			if (out != null) {
				try {
					out.close();
				} catch (IOException e) {
					VigorException.printExceptionMessage(e.getMessage());
				}
			}
		}

	}

	public void runblastCommand(String shellCommand, String outputXMLPath) throws IOException, InterruptedException {
		Runtime runtime = Runtime.getRuntime();
		Process process = runtime.exec(shellCommand);
		BlastParser parser = XmlFileBlastParser.create(process.getInputStream());
		System.out.println(parser.canParse());
		BlastVistorImpl blastImpl = new BlastVistorImpl();
		parser.parse(blastImpl);
		System.out.println(blastImpl.getBlastHitList());

	}*/

    /**
     * Note: Now we continue developing code assuming that we have reference database
     * @param form : Universal data holder
     */
	public void chooseAlignmentEvidence(VigorForm form) {




    }


}
