package org.jcvi.vigor.utils.TBLParser;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.util.iter.StreamingIterator;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.jcvi.vigor.component.Exon;

public class TBLFileParser {

	public List<TBLModel> getModels(String TBLFilePath, String PEPFilePath) {
		List<TBLModel> TBLModels = parseFile(TBLFilePath);
		TBLModels = setReferenceViralProteinID(TBLModels, PEPFilePath);
		return TBLModels;
	}

	public List<TBLModel> parseFile(String TBLFilePath) {
		List<TBLModel> models = new ArrayList<TBLModel>();
		try {
			Stream<String> tblFile = Files.lines(Paths.get(TBLFilePath));
			Pattern pattern;
			Matcher matcher;
			TBLModel model = null;
			List<Exon> exons = null;
			Exon exon = null;
			String virusGenomeID = "";
			Pattern startPattern = Pattern
					.compile("(\\s*)?([\\d*]+)([\\s*]+)(>)?([\\d*]+)(\\s*)?(CDS)");
			Pattern pattern10 = Pattern.compile("(\\s*)?([\\d*]+)([\\s*]+)([\\d*]+)(\\s*)?(misc_feature)");
			Pattern pattern2 = Pattern
					.compile("(<)(\\s*)?([\\d*]+)([\\s*]+)(>)([\\d*]+)(\\s*)?(CDS)");
			Pattern pattern3 = Pattern
					.compile("(<)(\\s*)?([\\d*]+)([\\s*]+)([\\d*]+)(\\s*)?(CDS)");
			Pattern pattern4 = Pattern.compile("(protein_id)([\\s*]+)(.*)");
			Pattern pattern5 = Pattern.compile("(product)([\\s*]+)(.*)");
			Pattern pattern6 = Pattern.compile("(note)([\\s*]+)(.*)");
			Pattern pattern7 = Pattern.compile("^(\\s*)?(gene)(\\s*)(.*)$");
			Pattern pattern9 = Pattern
					.compile("(\\s*)?([\\d*]+)([\\s*]+)([\\d*]+)(\\s*)?$");
			Pattern pattern8 = Pattern
					.compile("(\\s*)?([\\d*]+)([\\s*]+)(>)?([\\d*]+)(\\s*)?(gene)");
			// TBLFragment fragment = new TBLFragment();
			boolean isPseudoGene=false;
			for (String s : (Iterable<String>) tblFile::iterator) {
				if (s.startsWith(">")) {
					if (model != null && exons != null && exons.size() > 0) {
						model.setExons(exons);
						model.setPseudoGene(isPseudoGene);
						models.add(model);
						isPseudoGene=false;
					}
					model = null;
					pattern = Pattern.compile("Features[\\s](\\S*)");
					matcher = pattern.matcher(s);
					if (matcher.find()) {
						virusGenomeID = matcher.group(1).toString();
					}

				} else {

					matcher = startPattern.matcher(s);
					Matcher matcher2 = pattern2.matcher(s);
					Matcher matcher3 = pattern3.matcher(s);

					Matcher matcher4 = pattern4.matcher(s);

					Matcher matcher5 = pattern5.matcher(s);

					Matcher matcher6 = pattern6.matcher(s);

					Matcher matcher7 = pattern7.matcher(s);

					Matcher matcher8 = pattern8.matcher(s);
					Matcher matcher10 = pattern10.matcher(s);

					Matcher matcher9 = pattern9.matcher(s);

					if (matcher4.find() && model != null) {
						model.setViralProteinID((matcher4.group(3)));
					} else if (matcher5.find() && model != null) {
						model.setProduct(matcher5.group(3));
					} else if (matcher6.find() && model != null) {
						model.setNote(matcher6.group(3));
					} else if (matcher7.find() && model != null) {

						model.setGene(matcher7.group(4));

					} else if (matcher8.find()) {
						if (model != null && exons != null && exons.size() > 0) {
							model.setExons(exons);
							model.setPseudoGene(isPseudoGene);
							models.add(model);
							isPseudoGene=false;
						}
						model = new TBLModel();
						model.setVirusGenomeID(virusGenomeID);
						exons = new ArrayList<Exon>();
					} else if (matcher.matches()) {
						exon = new Exon();
						Range range = null;
						range = Range.of(Long.parseLong(matcher.group(2)),
								Long.parseLong(matcher.group(5)));
						exon.setRange(range);
						exons.add(exon);
					} else if (matcher10.matches()){
						exon = new Exon();
						Range range = null;
						range = Range.of(Long.parseLong(matcher10.group(2)),
								Long.parseLong(matcher10.group(4)));
						exon.setRange(range);
						exons.add(exon);	
						isPseudoGene = true;
					}
					  else if (matcher2.matches()) {
						exon = new Exon();
						Range range = null;
						range = Range.of(Long.parseLong(matcher2.group(3)),
								Long.parseLong(matcher2.group(6)));
						exon.setRange(range);
						exons.add(exon);
					} else if (matcher3.matches()) {
						exon = new Exon();
						Range range = null;
						range = Range.of(Long.parseLong(matcher3.group(3)),
								Long.parseLong(matcher3.group(5)));
						exon.setRange(range);
						exons.add(exon);
					} else if (matcher9.matches()) {
						exon = new Exon();
						Range range = null;
						range = Range.of(Long.parseLong(matcher9.group(2)),
								Long.parseLong(matcher9.group(4)));
						exon.setRange(range);
						exons.add(exon);
					}
				}
			}
			if (model != null && exons != null && exons.size() > 0) {
				model.setExons(exons);
				model.setPseudoGene(isPseudoGene);
				models.add(model);
				isPseudoGene = false;
			}
			tblFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		// map protein id's from pep file
		// models.stream().forEach(System.out::println);
		return models;
	}

	public static List<TBLModel> setReferenceViralProteinID(
			List<TBLModel> models, String pepFilePath) {
		try {
			File fastaFile = new File(pepFilePath);
			ProteinFastaDataStore dataStore = new ProteinFastaFileDataStoreBuilder(
					fastaFile).build();
			StreamingIterator<String> ids = dataStore.idIterator();
			HashMap<String, String> map = new HashMap<String, String>();
			while (ids.hasNext()) {
				ProteinFastaRecord record = dataStore.get(ids.next());
				String defline = record.getComment();
				Pattern pattern = Pattern.compile("ref_id=\"(\\S*)\"");
				Matcher matcher = pattern.matcher(defline);
				if (matcher.find()) {
					map.put(record.getId(), matcher.group(1));
				}
			}

			for (int i = 0; i < models.size(); i++) {

				String refProteinID = models.get(i).getViralProteinID();
				if (map.containsKey(refProteinID)) {
					models.get(i).setViralProteinID(map.get(refProteinID));
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		return models;
	}
}
