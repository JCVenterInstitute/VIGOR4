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

public class TBLFileParser {

	public List<TBLModel> vigor3Models() {
		String fileName = "C:/Users/snettem/Desktop/TEST.tbl";
		List<TBLModel> models = new ArrayList<TBLModel>();
		try {
			Stream<String> tblFile = Files.lines(Paths.get(fileName));
			Pattern pattern;
			Matcher matcher;
			TBLModel model = null;
			List<TBLFragment> listOfFragments = null;
			TBLFragment fragment = null;
			String virusGenomeID = "";
			boolean tempFlag = false;
			for (String s : (Iterable<String>) tblFile::iterator) {
				if (s.startsWith(">")) {
					model = new TBLModel();

					pattern = Pattern.compile("Features[\\s](\\S*)");
					matcher = pattern.matcher(s);
					if (matcher.find()) {
						model.setVirusGenomeID(matcher.group(1).toString());
						virusGenomeID = matcher.group(1).toString();
					}

				} else {                
					pattern = Pattern.compile("(\\s*)?(\\d*)(\\s*)(>)?(\\d*)(\\s*)?(CDS)");
					Pattern pattern2 = Pattern.compile("(<)(\\s*)?(\\d*)(\\s*)(>)(\\d*)(\\s*)?(CDS)");
					Pattern pattern3 = Pattern.compile("(<)(\\s*)?(\\d*)(\\s*)?(\\d*)(\\s*)?(CDS)");
					
					matcher = pattern.matcher(s);
					Matcher matcher2 = pattern2.matcher(s);
					Matcher matcher3 = pattern3.matcher(s);

					Pattern pattern4 = Pattern.compile("(protein_id)(\\s{1})(.*)");
					Matcher matcher4 = pattern4.matcher(s);
					
					Pattern pattern5 = Pattern.compile("(product)(\\s{1})(.*)");
					Matcher matcher5 = pattern5.matcher(s);

					Pattern pattern6 = Pattern.compile("(note)(\\s{1})(.*)");
					Matcher matcher6 = pattern6.matcher(s);
					
					Pattern pattern7 = Pattern.compile("^(\\s*?)(gene)(\\s{1})(.*)$");
					Matcher matcher7 = pattern7.matcher(s);
					
					Pattern pattern8 = Pattern.compile("(\\s*)?(\\d*)(\\s*)(\\d*)(\\s*)?(gene)");
					Matcher matcher8 = pattern8.matcher(s);
					
					Pattern pattern9 = Pattern.compile("(\\s*)?(\\d*)(\\s*)(\\d*)(\\s*)?$");
					Matcher matcher9 = pattern9.matcher(s);

					if (matcher4.find() && model != null) {
						fragment.setProteinID(matcher4.group(3));
					}
					else if (matcher5.find() && model != null) {
						fragment.setProduct(matcher5.group(3));
					}
					else if (matcher6.find() && model != null) {
						fragment.setNote(matcher6.group(3));
					}
					else if (matcher7.find() && model != null && fragment != null) {
						fragment.setGene(matcher7.group(4));
					}
					else if (matcher8.matches()) {
						if (tempFlag) {
							model.setTBLFragments(listOfFragments);
							models.add(model);
							model = new TBLModel();
							model.setVirusGenomeID(virusGenomeID);
							fragment = new TBLFragment();
						} else {
							tempFlag = true;
						}
					}
					else if (matcher.matches()) {
						listOfFragments = new ArrayList<TBLFragment>();
						fragment = new TBLFragment();
						Range range = null;
						range = Range.of(Long.parseLong(matcher.group(2)), Long.parseLong(matcher.group(5)));
						fragment.setRange(range);
						listOfFragments.add(fragment);
					} else if (matcher2.matches()) {
						fragment = new TBLFragment();
						Range range = null;
						range = Range.of(Long.parseLong(matcher2.group(3)), Long.parseLong(matcher2.group(6)));
						fragment.setRange(range);
						listOfFragments.add(fragment);
					} else if (matcher3.matches()) {
						fragment = new TBLFragment();
						Range range = null;
						range = Range.of(Long.parseLong(matcher3.group(3)), Long.parseLong(matcher3.group(5)));
						fragment.setRange(range);
						listOfFragments.add(fragment);
					}
					else if (matcher9.matches()){
						fragment = new TBLFragment();
						Range range = null;
						range = Range.of(Long.parseLong(matcher3.group(3)), Long.parseLong(matcher3.group(5)));
						fragment.setRange(range);
						listOfFragments.add(fragment);
					}
				}
			}
			model.setTBLFragments(listOfFragments);
			models.add(model);
			model = new TBLModel();
			model.setVirusGenomeID(virusGenomeID);
			tblFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		// map protein id's from pep file
		//models.stream().forEach(System.out::println);
		return models;
	}
	public static List<TBLModel> setReferenceViralProteinID(List<TBLModel> models, String pepFilePath) {
		try {
			File fastaFile = new File(pepFilePath);
			ProteinFastaDataStore dataStore = new ProteinFastaFileDataStoreBuilder(fastaFile).build();
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

				for (int j = 0; j < models.get(i).getTBLFragments().size(); j++) {

					String refProteinID = models.get(i).getTBLFragments().get(j).getProteinID();
					if (map.containsKey(refProteinID)) {
						models.get(i).getTBLFragments().get(j).setProteinID(map.get(refProteinID));
					}
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		return models;
	}
}
