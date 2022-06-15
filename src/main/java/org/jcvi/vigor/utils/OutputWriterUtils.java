package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.Sequence;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.MaturePeptideMatch;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.ViralProtein;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class OutputWriterUtils {

    private static Logger LOGGER = LogManager.getLogger(OutputWriterUtils.class);

    // Just for holding static functions, not to instantiate
    private OutputWriterUtils() {
    }

    public static String getGeneCoordinatesString (Model model, List<Model> geneModels ) {

        long seqLength = model.getAlignment().getVirusGenome().getSequence().getLength();
        long start = VigorFunctionalUtils.getDirectionBasedCoordinate(model.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED),
                                                                      seqLength,
                                                                      model.getDirection());
        /*Model endGeneModel = geneModels.stream()
                                       .filter(m -> model.getProteinID().equals(m.getProteinID()))
                                       .findFirst().orElse(model);*/
        long end = VigorFunctionalUtils.getDirectionBasedCoordinate(model.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED),
                                                                    seqLength,
                                                                    model.getDirection());
        return String.join("\t",
                           ( model.isPartial5p() ? "<" : "" ) + start,
                           ( model.isPartial3p() ? ">" : "" ) + end);
    }

    public static String getDefline (Model model ) {

        long seqLength = model.getAlignment().getVirusGenome().getSequence().getLength();
        ViralProtein refProtein = model.getAlignment().getViralProtein();
        StringBuilder defline = new StringBuilder();
        defline.append(">" + model.getGeneID());
        if (model.isPseudogene()) {
            defline.append(" pseudogene");
        }
        int codon_start = model.getExons().get(0).getFrame().getFrame();
        defline.append(" location=");
        for (int i = 0; i < model.getExons().size(); i++) {
            Exon exon = model.getExons().get(i);
            String start = Long.toString(VigorFunctionalUtils.getDirectionBasedCoordinate(
                    exon.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED),
                    seqLength,
                    model.getDirection()));
            String end = Long.toString(VigorFunctionalUtils.getDirectionBasedCoordinate(
                    exon.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED),
                    seqLength,
                    model.getDirection()));
            if (model.isPartial5p() && i == 0) {
                start = "<" + start;
            }
            if (model.isPartial3p() && i == model.getExons().size() - 1) {
                end = ">" + end;
            }
            if (i != 0) defline.append(",");
            defline.append(String.format("%s..%s", start, end));
        }
        defline.append(" codon_start=" + codon_start);
        defline.append(String.format(" gene=\"%s\"", refProtein.getGeneSymbol()));
        String product = refProtein.getProduct();
        if (! NullUtil.isNullOrEmpty(product)) {
            defline.append(String.format(" product=\"%s\"", VigorUtils.putativeName(product, model.isPartial3p(), model.isPartial5p())));
        } else {
            LOGGER.warn("Missing product for protein {}", refProtein.getProteinID());
        }
        String reference_db = model.getAlignment().getAlignmentEvidence().getReference_db();
        if (! NullUtil.isNullOrEmpty(reference_db) ) {
            defline.append(String.format(" ref_db=\"%s\"", Paths.get(reference_db).getFileName().toString()));
        }
        defline.append(String.format(" ref_id=\"%s\"", refProtein.getProteinID()));
        return defline.toString();
    }

    public static String formatMaturePeptideRange ( Model model,
                                                    MaturePeptideMatch match,
                                                    List<Range> ranges,
                                                    Range.CoordinateSystem coordinateSystem,
                                                    String rangeDelimiter,
                                                    long startCoordinate,
                                                    long endCoordinate, boolean CDSRanges) {
        long seqLength = CDSRanges ? model.getAlignment().getVirusGenome().getSequence().getLength()
                : model.getAlignment().getViralProtein().getSequence().getLength();
        List<String> rangeStrings = new ArrayList<>(ranges.size());
        Exon initialExon = model.getExons().get(0);
        Exon lastExon = model.getExons().get(model.getExons().size() - 1);
        for (Range range : ranges) {
            LOGGER.trace("range {}-{} start {} sframe {} end {} eframe {} 5p {} 3p {}",
                         range.getBegin(coordinateSystem),
                         range.getEnd(coordinateSystem),
                         startCoordinate,
                         initialExon.getFrame().getFrame(),
                         endCoordinate,
                         lastExon.getFrame().getFrame(),
                         model.isPartial5p(),
                         model.isPartial3p());
            long start = VigorFunctionalUtils.getDirectionBasedCoordinate(range.getBegin(coordinateSystem),seqLength,model.getDirection());
            String startStr = String.valueOf(start);
            if (match.isFuzzyBegin() || model.isPartial5p() && start == startCoordinate) {
                startStr = "<" + startStr;
            }
            long end = VigorFunctionalUtils.getDirectionBasedCoordinate(range.getEnd(coordinateSystem),seqLength,model.getDirection());
            String endStr = String.valueOf(end);
            if (match.isFuzzyEnd() || ( model.isPartial3p() && end == endCoordinate )) {
                endStr = ">" + endStr;
            }
            rangeStrings.add(String.format("%s%s%s", startStr, rangeDelimiter, endStr));
        }
        return String.join(",", rangeStrings);
    }

        public static void writeSequence (WriterBundle bw, Sequence seq ) throws IOException {

        Iterator<String> sequenceLineIter = SequenceUtils.steamOf(seq, 60).iterator();
        while (sequenceLineIter.hasNext()) {
            bw.write(sequenceLineIter.next());
            bw.newLine();
        }
    }


    public static String getSequenceFilePath(String sequenceID) {
        String fileName = sequenceID.replaceAll("\\.", "-");
        fileName = fileName.split("-", 2)[0];
        fileName = fileName.replaceAll("[:/]", "_");
        LOGGER.trace("filename for sequence ID {} is {}", sequenceID, fileName);
        return fileName;
    }

    public static String getGeneFilePath(String geneID) {
        String fileName = geneID.replaceAll("[.:/]", "_");
        LOGGER.trace("filename for geneID {} is {}", geneID, fileName);
        return fileName;
    }


}
