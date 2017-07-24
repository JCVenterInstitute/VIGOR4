package com.vigor.utils.TBLParser;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.jcvi.jillion.core.util.iter.StreamingIterator;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;

public class TBLVigorParser {
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		TBLFileParser fileparser = new TBLFileParser();
		List<TBLModel> models = setReferenceViralProteinID(fileparser.vigor3Models());
		models.stream().forEach(System.out::println);

	}
	
	
	

	public static List<TBLModel> setReferenceViralProteinID(List<TBLModel> models) {
		try {
			File fastaFile = new File("C:/Users/snettem/Desktop/TEST.pep");
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
