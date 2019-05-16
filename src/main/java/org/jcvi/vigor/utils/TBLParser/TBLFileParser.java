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
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.util.iter.StreamingIterator;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.jcvi.vigor.component.Exon;

public class TBLFileParser {

    public List<TBLModel> getModels ( String TBLFilePath ) {

        List<TBLModel> TBLModels = parseFile(TBLFilePath);
        //TBLModels = setReferenceViralProteinID(TBLModels, PEPFilePath);
        return TBLModels;
    }

    public List<TBLModel> parseFile ( String TBLFilePath ) {

        List<TBLModel> models = new ArrayList<>();
        try {
            Stream<String> tblFile = Files.lines(Paths.get(TBLFilePath));
            Pattern pattern;
            Matcher matcher;
            TBLModel model = null;
            List<Exon> exons = null;
            Exon exon = null;
            String virusGenomeID = "";
            Pattern startPattern = Pattern.compile("(<)?(\\s*)?([\\d*]+)([\\s*]+)(>)?([\\d*]+)(\\s*)?(CDS)(\\s*)?$");
            Pattern proteinIdPattern = Pattern.compile("(\\s*)?(protein_id)([\\s*]+)(.*)(\\s*)?$");
            Pattern productPattern = Pattern.compile("(\\s*)?(product)([\\s*]+)(.*)(\\s*)?$");
            Pattern notePattern = Pattern.compile("(\\s*)?(note)([\\s*]+)(.*)(\\s*)?$");
            Pattern geneNamePattern = Pattern.compile("^(\\s*)?(gene)(\\s*)(.*)$");
            Pattern codonStartPattern = Pattern.compile("^(\\s*)?(codon_start)(\\s*)?([0-9])(\\s*)?$");
            Pattern nextFragment = Pattern
                    .compile("(\\s*)?([\\d*]+)([\\s*]+)(>)?([\\d*]+)(\\s*)?$");
            Pattern pseudogenePattern = Pattern.compile("(\\s*)?(pseudogene)(.*)(\\s*)?$");
            Pattern riboSlippagePattern = Pattern.compile("(\\s*)?(ribosomal(_?\\s?)slippage)(.*)(\\s*)?$");
            Pattern stopReadThroughPattern = Pattern.compile("(\\s*)?(transl_except)(\\s*)(\\((pos:)(\\d*)(..)(\\d*)(,.*)\\))(\\s*)?$");
            Pattern geneLinePattern = Pattern
                    .compile("(<)?(\\s*)?([\\d*]+)([\\s*]+)(>)?([\\d*]+)(\\s*)?(gene)(\\s*)?$");
            Pattern miscFeaturePattern = Pattern.compile("(<)?(\\s*)?([\\d*]+)([\\s*]+)(>)?([\\d*]+)(\\s*)?(misc_feature)(\\s*)?$");
            boolean isPseudoGene = false;
            boolean is5Partial = false;
            boolean is3Partial = false;
            boolean isRiboSlippage = false;
            Range stopCodonReadThrough = null;
            for (String s : (Iterable<String>) tblFile:: iterator) {
                Matcher geneLineMatcher = geneLinePattern.matcher(s);
                if (s.startsWith(">") || geneLineMatcher.matches()) {
                    if (model != null && exons != null && exons.size() > 0) {
                        model.setExons(exons);
                        model.setPseudoGene(isPseudoGene);
                        model.set5Partial(is5Partial);
                        model.set3Partial(is3Partial);
                        model.setRiboSlippage(isRiboSlippage);
                        model.setStopCodonReadThrough(stopCodonReadThrough);
                        models.add(model);
                        isPseudoGene = false;
                        is5Partial = false;
                        is3Partial = false;
                        isRiboSlippage = false;
                        stopCodonReadThrough = null;
                    }
                    model = null;
                    pattern = Pattern.compile("Feature(s?)[\\s](\\S*)");
                    matcher = pattern.matcher(s);
                    if (matcher.find()) {
                        virusGenomeID = matcher.group(2);
                    }
                    if (geneLineMatcher.matches()) {
                        model = new TBLModel();
                        model.setVirusGenomeID(virusGenomeID);
                        exons = new ArrayList<Exon>();
                    }
                } else {
                    matcher = startPattern.matcher(s);
                    Matcher miscMatcher = miscFeaturePattern.matcher(s);
                    Matcher proteinIDMatcher = proteinIdPattern.matcher(s);
                    Matcher productMatcher = productPattern.matcher(s);
                    Matcher noteMatcher = notePattern.matcher(s);
                    Matcher geneNameMatcher = geneNamePattern.matcher(s);
                    Matcher riboSlippageMatcher = riboSlippagePattern.matcher(s);
                    Matcher nextFragMatcher = nextFragment.matcher(s);
                    Matcher pseudogeneMatcher = pseudogenePattern.matcher(s);
                    Matcher stopReadThroughMatcher = stopReadThroughPattern.matcher(s);
                    Matcher codonStartMatcher = codonStartPattern.matcher(s);
                    if (proteinIDMatcher.find() && model != null) {
                        model.setViralProteinID(( proteinIDMatcher.group(4) ));
                        model.setGeneID(proteinIDMatcher.group(4));
                    } else if (productMatcher.find() && model != null) {
                        model.setProduct(productMatcher.group(4));
                    } else if (noteMatcher.find() && model != null) {
                        model.setNote(noteMatcher.group(4));
                    } else if (geneNameMatcher.find() && model != null) {
                        model.setGene(geneNameMatcher.group(4));
                    } else if (matcher.matches()) {
                        exon = new Exon();
                        Range range;
                        exon.setFrame(Frame.ONE);
                        range = Range.of(Long.parseLong(matcher.group(3)),
                                Long.parseLong(matcher.group(6)));
                        if (matcher.group(1) != null && matcher.group(1).equals("<")) {
                            is5Partial = true;
                        }
                        if (matcher.group(5) != null && matcher.group(5).equals(">")) {
                            is3Partial = true;
                        }
                        exon.setRange(range);
                        exons.add(exon);
                    } else if (codonStartMatcher.matches()) {
                        Frame frame = Frame.parseFrame(Integer.parseInt(codonStartMatcher.group(4)));
                        exons.get(0).setFrame(frame);
                    } else if (miscMatcher.matches()) {
                        isPseudoGene = true;
                        if (isPseudoGene) {
                            Range range = Range.of(Long.parseLong(miscMatcher.group(3)),
                                    Long.parseLong(miscMatcher.group(6)));
                            exon = new Exon();
                            exon.setRange(range);
                            exons.add(exon);
                            if (miscMatcher.group(1) != null && miscMatcher.group(1).equals("<")) {
                                is5Partial = true;
                            }
                            if (miscMatcher.group(5) != null && miscMatcher.group(5).equals(">")) {
                                is3Partial = true;
                            }
                        }
                    } else if (nextFragMatcher.matches()) {
                        exon = new Exon();
                        Range range;
                        range = Range.of(Long.parseLong(nextFragMatcher.group(2)),
                                Long.parseLong(nextFragMatcher.group(5)));
                        if (nextFragMatcher.group(4) != null && nextFragMatcher.group(4).equals(">")) {
                            is3Partial = true;
                        }
                        exon.setRange(range);
                        exons.add(exon);
                    } else if (pseudogeneMatcher.matches()) {
                        isPseudoGene = true;
                    } else if (riboSlippageMatcher.find()) {
                        isRiboSlippage = true;
                    } else if (stopReadThroughMatcher.matches()) {
                        stopCodonReadThrough = Range.of(Long.parseLong(stopReadThroughMatcher.group(6)),
                                Long.parseLong(stopReadThroughMatcher.group(8)));
                    }
                }
            }
            if (model != null && exons != null && exons.size() > 0) {
                model.setExons(exons);
                model.setPseudoGene(isPseudoGene);
                model.set5Partial(is5Partial);
                model.set3Partial(is3Partial);
                model.setRiboSlippage(isRiboSlippage);
                model.setStopCodonReadThrough(stopCodonReadThrough);
                models.add(model);
            }
            tblFile.close();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return models;
    }

    public static List<TBLModel> setReferenceViralProteinID (
            List<TBLModel> models, String pepFilePath ) {

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
