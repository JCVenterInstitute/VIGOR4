package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.Sequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.service.Scores;
import org.springframework.stereotype.Service;

import java.io.*;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

@Service
public class GenerateVigorOutput {

    public enum Outfile {
        TBL("tbl"),
        CDS("cds"),
        GFF3("gff3"),
        PEP("pep"),
        ALN("aln"),
        SUM("sum");

        final public String extension;

        Outfile ( String extension ) {

            this.extension = extension;
        }
    }

    public static class Outfiles extends EnumMap<Outfile, BufferedWriter> implements AutoCloseable {

        public Outfiles () {

            super(Outfile.class);
        }

        public void close () throws IOException {

            List<IOException> exceptions = new ArrayList<>();
            for (BufferedWriter writer : values()) {
                try {
                    writer.close();
                } catch (IOException e) {
                    exceptions.add(e);
                }
            }
            if (!exceptions.isEmpty()) {
                // TODO join all exceptions somehow meaningfully
                throw exceptions.get(0);
            }
        }

        public void flush () throws IOException {

            for (BufferedWriter writer : values()) {
                writer.flush();
            }
        }
    }

    public static class WriterBundle {
        private final BufferedWriter[] writers;

        public WriterBundle(BufferedWriter... writers) {
            this.writers = writers;
        }

        public void write(String value) throws IOException {
            for (BufferedWriter w : writers) {
                w.write(value);
            }
        }

        public void flush() throws IOException {
            for (BufferedWriter w : writers) {
                w.flush();
            }
        }

        public void write(char c) throws IOException {
            for (BufferedWriter w: writers) {
                w.write(c);
            }
        }

        public void newLine() throws IOException {
            for (BufferedWriter w: writers) {
                w.newLine();
            }
        }
    }

    private static final Logger LOGGER = LogManager.getLogger(GenerateVigorOutput.class);
    private static Range.CoordinateSystem oneBased = Range.CoordinateSystem.RESIDUE_BASED;

    public void generateOutputFiles ( VigorConfiguration config, Outfiles outfiles, List<Model> geneModels ) throws IOException {

        generateTBLReport(config, outfiles.get(Outfile.TBL), geneModels);
        generateCDSReport(config, outfiles.get(Outfile.CDS), geneModels);
        generatePEPReport(config, outfiles.get(Outfile.PEP), geneModels);
        generateSUMReport(config, outfiles.get(Outfile.SUM), geneModels);
    }

    public void generateSUMReport(VigorConfiguration config, BufferedWriter bw, List<Model> geneModels) throws IOException {

        if (geneModels.isEmpty()) {
            LOGGER.warn("no gene models to write to file");
            return;
        }
        String baseSequenceFileName = getSequenceFilePath(geneModels.get(0).getAlignment().getVirusGenome().getId());
        String outputDir = config.get(ConfigurationParameters.OutputDirectory);

        String sequenceSummaryReport = Paths.get(outputDir, baseSequenceFileName + ".sum").toString();

        try (FileWriter fw = new FileWriter(sequenceSummaryReport);
             BufferedWriter bwSequence = new BufferedWriter(fw)) {

            double identityAvg = 0;
            double similarityAvg = 0;
            double coverageAvg = 0;
            long totalCDSBases = 0;
            long totalPepBases = 0;
            VirusGenome virusGenome = geneModels.get(0).getAlignment().getVirusGenome();
            long seqLength = virusGenome.getSequence().getLength();
            StringBuffer content = new StringBuffer("");
            content.append("gene_id\t%identity\t%similarity\t%coverage\tstart..stop\tpep_size\tref_size\tref_id\tgene\tgene_product\n");
            for (Model model : geneModels) {
                ViralProtein viralProtein = model.getAlignment().getViralProtein();
                Map<String, Double> scores = model.getScores();
                long cdsBases = 0;
                content.append(model.getGeneID());
                content.append("\t" + String.format("%.02f", scores.get(Scores.IDENTITY_SCORE)));
                content.append("\t" + String.format("%.02f", scores.get(Scores.SIMILARITY_SCORE)));
                content.append("\t" + String.format("%.02f", scores.get(Scores.COVERAGE_SCORE)) + "\t");
                for (int i = 0; i < model.getExons().size(); i++) {
                    Exon exon = model.getExons().get(i);
                    String start = Long.toString(exon.getRange().getBegin(oneBased));
                    String end = Long.toString(exon.getRange().getEnd(oneBased));
                    if (i == 0 && model.isPartial5p()) {
                        start = "<" + start;
                    }
                    if (i == model.getExons().size() - 1 && model.isPartial3p()) {
                        end = ">" + end;
                    }
                    if (model.getExons().size() > 1 && i != 0) {
                        content.append(", ");
                    }
                    content.append(start + ".." + end);
                    cdsBases = cdsBases + exon.getRange().getLength();
                }
                if (!model.isPartial3p()) {
                    cdsBases = cdsBases - 3;
                }
                content.append("\t" + (cdsBases) / 3);
                content.append("\t" + viralProtein.getSequence().getLength());
                content.append("\t" + viralProtein.getProteinID());
                content.append("\t" + model.getGeneSymbol());
                content.append("\t" + viralProtein.getProduct());

                content.append(System.lineSeparator());
                totalCDSBases = totalCDSBases + cdsBases;
                totalPepBases = cdsBases + totalPepBases;
                identityAvg = identityAvg + scores.get(Scores.IDENTITY_SCORE);
                similarityAvg = similarityAvg + scores.get(Scores.SIMILARITY_SCORE);
                coverageAvg = coverageAvg + scores.get(Scores.COVERAGE_SCORE);

                IDGenerator idGenerator = IDGenerator.of(model.getGeneID());
                for (MaturePeptideMatch match : model.getMaturePeptides()) {
                    content.append(idGenerator.next());
                    content.append("\t" + String.format("%.02f", match.getIdentity() * 100));
                    content.append("\t" + String.format("%.02f", match.getSimilarity() * 100));
                    content.append("\t" + String.format("%.02f", match.getCoverage() * 100) + "\t");
                    List<Range> cdRanges = VigorFunctionalUtils.proteinRangeToCDSRanges(model, match.getProteinRange());
                    long start = VigorFunctionalUtils.getDirectionBasedCoordinate(
                            model.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqLength, model.getDirection());
                    long end = VigorFunctionalUtils.getDirectionBasedCoordinate(
                            model.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqLength, model.getDirection());

                    content.append(String.format(
                            GenerateVigorOutput.formatMaturePeptideRange(model, match, cdRanges, Range.CoordinateSystem.RESIDUE_BASED,
                                                                         "..", start + model.getExons().get(0).getFrame().getFrame() - 1, end, true)));
                    content.append("\t" + match.getProteinRange().getLength());
                    content.append("\t" + match.getReference().getSequence().getLength());
                    content.append("\t" + match.getReference().getProteinID());
                    content.append("\t" + model.getGeneSymbol());
                    content.append("\t"
                                           + VigorUtils.putativeName(match.getReference().getProduct(), match.isFuzzyEnd(), match.isFuzzyBegin()));
                    content.append(System.lineSeparator());
                }
            }
            bw.write(content.toString());
            bwSequence.write(content.toString());
        }

    }

    public void generateTBLReport ( VigorConfiguration config, BufferedWriter tblWriter, List<Model> geneModels ) throws IOException {

        if (geneModels.isEmpty()) {
            LOGGER.warn("no gene models to write to file");
            return;
        }

        String outputDir = config.get(ConfigurationParameters.OutputDirectory);

        VirusGenome virusGenome = geneModels.get(0).getAlignment().getVirusGenome();
        String genomeID = virusGenome.getId();

        String locusPrefix = config.get(ConfigurationParameters.Locustag);
        boolean writeLocus = ! NullUtil.isNullOrEmpty(locusPrefix);
        long seqlength = virusGenome.getSequence().getLength();

        String geneFileBase = getSequenceFilePath(genomeID);
        String geneFilePath = Paths.get(outputDir, geneFileBase + ".tbl").toString();
        try (FileWriter geneFileWriter = new FileWriter(geneFilePath);
             BufferedWriter geneFileBWWriter = new BufferedWriter(geneFileWriter)) {

            WriterBundle bw = new WriterBundle(tblWriter, geneFileBWWriter);

            bw.write(">Features " + genomeID + "\n");
            String proteinID = "";
            for (int i = 0; i < geneModels.size(); i++) {
                Model model = geneModels.get(i);
                List<String> modelNotes = model.getNotes();
                Ribosomal_Slippage riboSlippage = model.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage();
                RNA_Editing rna_editing = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing();
                List<SpliceSite> spliceSites = model.getAlignment().getViralProtein().getGeneAttributes().getSpliceSites();
                StringBuilder notes = new StringBuilder("");
                List<Exon> exons = model.getExons();
                Collections.sort(exons, Comparator.comparing(e -> VigorFunctionalUtils.getDirectionBasedRange(e.getRange(), seqlength, model.getDirection()), Range.Comparators.ARRIVAL));
                Exon firstExon = exons.get(0);
                int codon_start = firstExon.getFrame().getFrame();
                if (!model.getAlignment().getViralProtein().getProteinID().equals(proteinID)) {
                    bw.write(getGeneCoordinatesString(model, geneModels));
                    bw.write("\tgene\n");
                    if (writeLocus) {
                        bw.write("\t\t\tlocus_tag\t" + VigorUtils.nameToLocus(model.getGeneSymbol(), locusPrefix, model.isPseudogene()) + "\n");
                    }
                    bw.write("\t\t\tgene\t" + model.getGeneSymbol() + "\n");
                }
                proteinID = model.getAlignment().getViralProtein().getProteinID();
                String geneSynonym = model.getAlignment().getViralProtein().getGeneSynonym();
                if (geneSynonym != null && geneSynonym != "") {
                    bw.write("\t\t\tgene_syn\t" + geneSynonym + "\n");
                }
                for (int j = 0; j < exons.size(); j++) {
                    Exon exon = exons.get(j);
                    String Cstart = Long.toString(VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getBegin(oneBased), seqlength, model.getDirection()));
                    String Cend = Long.toString(VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getEnd(oneBased), seqlength, model.getDirection()));
                    if (j == 0 && model.isPartial5p()) {
                        Cstart = "<" + Cstart;
                    }
                    if (j == exons.size() - 1 && model.isPartial3p()) {
                        Cend = ">" + Cend;
                    }
                    if (j == 0) {
                        bw.write(String.join("\t",
                                             Cstart,
                                             Cend,
                                             model.isPseudogene() ? "misc_feature" : "CDS"));
                        bw.newLine();
                    } else {
                        bw.write(Cstart + "\t" + Cend + "\n");
                    }
                }
                bw.write("\t\t\tcodon_start\t" + codon_start + "\n");
                if (model.getReplaceStopCodonRange() != null) {
                    long replaceStopBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(model.getReplaceStopCodonRange().getBegin(oneBased), seqlength, model.getDirection());
                    long replaceStopEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(model.getReplaceStopCodonRange().getEnd(oneBased), seqlength, model.getDirection());
                    bw.write("\t\t\ttransl_except\t" + String.format("(pos:%s..%s,aa:R)", replaceStopBegin, replaceStopEnd) + "\n");
                }
                bw.write("\t\t\tprotein_id\t" + model.getGeneID() + "\n");
                if (writeLocus) {
                    bw.write("\t\t\tlocus_tag\t" + VigorUtils.nameToLocus(model.getGeneSymbol(), locusPrefix, model.isPseudogene()) + "\n");
                }
                bw.write("\t\t\tgene\t" + model.getGeneSymbol() + "\n");
                String product = model.getAlignment().getViralProtein().getProduct();
                if (!NullUtil.isNullOrEmpty(product)) {
                    bw.write("\t\t\tproduct\t" + VigorUtils.putativeName(product, model.isPartial3p(), model.isPartial5p()) + "\n");
                } else {
                    LOGGER.warn("Missing product for {}", genomeID);
                }
                if (riboSlippage.isHas_ribosomal_slippage()) {
                    bw.write("\t\t\tribosomal_slippage\n");
                }
                if (rna_editing.isHas_RNA_editing()) {
                    bw.write("\t\t\texception\tRNA editing\n");
                    notes = notes.append(rna_editing.getNote() + ";");
                }
                if (spliceSites != SpliceSite.DEFAULT_SPLICE_SITES) {
                    notes.append("non-canonical splicing");
                }
                if (model.getInsertRNAEditingRange() != null) {
                    bw.write("\t\t\tnote\t" + notes + "\n");
                    long insertBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(model.getInsertRNAEditingRange().getBegin(oneBased), seqlength, model.getDirection());
                    long insertEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(model.getInsertRNAEditingRange().getEnd(oneBased), seqlength, model.getDirection());
                    // TODO coordinate system?
                    bw.write(insertBegin + "\t" + insertEnd + "\t" + "misc_feature\n");
                    NucleotideSequence subSeq = model.getAlignment().getVirusGenome().getSequence().toBuilder(model.getInsertRNAEditingRange()).build();
                    bw.write("\t\t\tnote\tlocation of RNA editing (" + subSeq + "," + rna_editing.getInsertionString() + ") in " + model.getAlignment().getViralProtein().getProduct() + "\n");
                }
                if (modelNotes.size() > 0) {
                    String notesText = String.join(",", modelNotes);
                    bw.write("\t\t\tnote\t" + notesText + "\n");
                }
            }
            for (Model model : geneModels) {
                if (model.getMaturePeptides() != null && !model.getMaturePeptides().isEmpty()) {
                    bw.write(">Features " + model.getGeneID());
                    long proteinLength = model.getAlignment().getViralProtein().getSequence().getLength();
                    bw.newLine();
                    String product;
                    for (MaturePeptideMatch match : model.getMaturePeptides()) {
                        long start = VigorFunctionalUtils.getDirectionBasedCoordinate(1, proteinLength, model.getDirection());
                        long end = VigorFunctionalUtils.getDirectionBasedCoordinate(proteinLength, proteinLength, model.getDirection());
                        bw.write(formatMaturePeptideRange(model,
                                                          match,
                                                          Arrays.asList(match.getProteinRange()),
                                                          Range.CoordinateSystem.RESIDUE_BASED,
                                                          "\t",
                                                          start,
                                                          end, false));
                        bw.write("\t");
                        product = match.getReference().getProduct();
                        if (!NullUtil.isNullOrEmpty(product)) {
                            if (product.contains("signal")) {
                                // TODO pre-classify type
                                bw.write("sig_peptide");
                            } else {
                                bw.write("mat_peptide");
                                bw.newLine();
                                bw.write("\t\t\tproduct\t");
                                // TODO check that there aren't other factors here.
                                bw.write(VigorUtils.putativeName(product, match.isFuzzyEnd(), match.isFuzzyBegin()));
                            }
                        } else {
                            LOGGER.warn("Missing product for mature peptide {}", match.getReference().getProteinID());
                        }
                        bw.newLine();
                        String geneSymbol = model.getGeneSymbol();
                        if (!NullUtil.isNullOrEmpty(geneSymbol)) {
                            if (writeLocus) {
                                bw.write("\t\t\tlocus_tag\t");
                                bw.write(VigorUtils.nameToLocus(geneSymbol, locusPrefix, model.isPseudogene()));
                                bw.newLine();
                            }
                            bw.write("\t\t\tgene\t");
                            bw.write(geneSymbol);
                            bw.newLine();
                        }
                    }
                }
            }
        }
    }

    private String getGeneCoordinatesString ( Model model, List<Model> geneModels ) {

        long seqLength = model.getAlignment().getVirusGenome().getSequence().getLength();
        long start = VigorFunctionalUtils.getDirectionBasedCoordinate(model.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED),
                                                                      seqLength,
                                                                      model.getDirection());
        Model endGeneModel = geneModels.stream()
                .filter(m -> model.getProteinID().equals(m.getProteinID()))
                .findFirst().orElse(model);
        long end = VigorFunctionalUtils.getDirectionBasedCoordinate(endGeneModel.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED),
                                                                    seqLength,
                                                                    model.getDirection());
        return String.join("\t",
                ( model.isPartial5p() ? "<" : "" ) + start,
                ( endGeneModel.isPartial3p() ? ">" : "" ) + end);
    }

    private String getDefline (Model model ) {

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
            String start = Long.toString(VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getBegin(oneBased),seqLength,model.getDirection()));
            String end = Long.toString(VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getEnd(oneBased),seqLength,model.getDirection()));
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

    private void writeSequence ( WriterBundle bw, Sequence seq ) throws IOException {

        Iterator<String> sequenceLineIter = SequenceUtils.steamOf(seq, 60).iterator();
        while (sequenceLineIter.hasNext()) {
            bw.write(sequenceLineIter.next());
            bw.newLine();
        }
    }

    public void generateCDSReport ( VigorConfiguration config, BufferedWriter cdsWriter, List<Model> geneModels ) throws IOException {

        if (geneModels.isEmpty()) {
            LOGGER.warn("No gene models to write CDS reports");
            return;
        }

        String outputDir = config.get(ConfigurationParameters.OutputDirectory);
        String sequenceFileBase = getSequenceFilePath(geneModels.get(0).getAlignment().getVirusGenome().getId());
        String sequenceFilePath = Paths.get(outputDir, sequenceFileBase + ".cds").toString();
        try (FileWriter sequenceFileWriter = new FileWriter(sequenceFilePath);
             BufferedWriter sequenceWriter = new BufferedWriter(sequenceFileWriter)) {

                for (Model model : geneModels) {
                    String geneFilePath = Paths.get(outputDir, getGeneFilePath(model.getGeneID()) + ".cds").toString();
                    try (FileWriter modelFileWriter = new FileWriter(geneFilePath);
                    BufferedWriter modelWriter = new BufferedWriter(modelFileWriter)) {
                    WriterBundle bw = new WriterBundle(cdsWriter, sequenceWriter, modelWriter);

                    bw.write(getDefline(model));
                    bw.newLine();
                    NucleotideSequenceBuilder builder = new NucleotideSequenceBuilder();
                    NucleotideSequence virusGenome = model.getAlignment().getVirusGenome().getSequence();
                    long seqLength = virusGenome.getLength();
                    List<Range> translatedRanges = model.getExons()
                                                        .stream()
                                                        .map(e -> VigorFunctionalUtils.getDirectionBasedRange(e.getRange(), seqLength, model.getDirection()))
                                                        .collect(Collectors.toList());
                    Collections.sort(translatedRanges, Range.Comparators.ARRIVAL);
                    for (Range exonRange : translatedRanges) {
                        builder.append(virusGenome.toBuilder(exonRange).build());
                    }

                    writeSequence(bw, builder.build());
                }
            }
            }
        }

    public void generatePEPReport ( VigorConfiguration config, BufferedWriter pepWriter, List<Model> geneModels ) throws IOException {

        StringBuilder defline;

        if (geneModels.isEmpty()) {
            LOGGER.warn("no geneModels to write to PEP report");
            return;
        }
        VirusGenome genome = geneModels.get(0).getAlignment().getVirusGenome();
        String genomeID = genome.getId();
        long seqLength = genome.getSequence().getLength();
        String outputDir = config.get(ConfigurationParameters.OutputDirectory);

        String sequenceFilePath = Paths.get(outputDir, getSequenceFilePath(genomeID) + ".pep").toString();
        try (FileWriter fw = new FileWriter(sequenceFilePath);
             BufferedWriter sequenceWriter = new BufferedWriter(fw)) {

            for (Model model : geneModels) {

                String geneFilePath = Paths.get(outputDir, getGeneFilePath(model.getGeneID()) + ".pep").toString();
                try (FileWriter modelFileWriter = new FileWriter(geneFilePath);
                     BufferedWriter modelWriter = new BufferedWriter(modelFileWriter)) {
                    WriterBundle fbw = new WriterBundle(pepWriter, sequenceWriter, modelWriter);

                    fbw.write(getDefline(model));
                    fbw.newLine();
                    writeSequence(fbw, model.getTranslatedSeq());
                    IDGenerator idGenerator = IDGenerator.of(model.getGeneID());
                    for (MaturePeptideMatch match : model.getMaturePeptides()) {
                        String pepID = idGenerator.next();
                        String pepFilePath = Paths.get(outputDir, getGeneFilePath(pepID) + ".pep").toString();
                        try (FileWriter pepFileWriter = new FileWriter(pepFilePath);
                             BufferedWriter pWriter = new BufferedWriter(pepFileWriter)) {
                            WriterBundle bw = new WriterBundle(pepWriter, sequenceWriter, pWriter);
                            defline = new StringBuilder();
                            defline.append(">" + pepID);
                            if (model.isPseudogene()) {
                                defline.append(" pseudogene");
                            }
                            defline.append(" mat_peptide");
                            List<Range> cdsRanges = VigorFunctionalUtils.proteinRangeToCDSRanges(model, match.getProteinRange());
                            // TODO handle truncation etc
                            Exon initialExon = model.getExons().get(0);
                            long startCoordinate = VigorFunctionalUtils.getDirectionBasedCoordinate(initialExon.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqLength, model.getDirection());
                            long endCoordinate = VigorFunctionalUtils.getDirectionBasedCoordinate(model.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqLength, model.getDirection());
                            defline.append(String.format(" location=%s", formatMaturePeptideRange(model,
                                                                                                  match,
                                                                                                  cdsRanges,
                                                                                                  Range.CoordinateSystem.RESIDUE_BASED,
                                                                                                  "..",
                                                                                                  // start_codon adjustment
                                                                                                  startCoordinate + initialExon.getFrame().getFrame() - 1,
                                                                                                  endCoordinate, true)));
                            defline.append(String.format(" gene=\"%s\"", model.getGeneSymbol()));
                            String product = match.getReference().getProduct();
                            if (!NullUtil.isNullOrEmpty(product)) {
                                defline.append(String.format(" product=\"%s\"", VigorUtils.putativeName(product, match.isFuzzyEnd(), match.isFuzzyBegin())));
                            } else {
                                LOGGER.warn("Missing product for {}", idGenerator.getCurrent());
                            }
                            String refDB = model.getAlignment().getAlignmentEvidence().getMatpep_db();
                            if (!NullUtil.isNullOrEmpty(refDB)) {
                                defline.append(String.format(" ref_db=\"%s\"", Paths.get(refDB).getFileName().toString()));
                            }
                            defline.append(String.format(" ref_id=\"%s\"", match.getReference().getProteinID()));
                            bw.write(defline.toString());
                            bw.newLine();
                            writeSequence(bw, match.getProtein().toBuilder().trim(match.getProteinRange()).build());
                        }
                    }
                }
            }
        }
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

    public static String getSequenceFilePath(String sequenceID) {
        String fileName = sequenceID.replaceAll("\\.", "-");
        fileName = fileName.split("-", 2)[0];
        return fileName.replaceAll("[:/]", "_");
    }

    public static String getGeneFilePath(String geneID) {
        return geneID.replaceAll("[.:/]", "_");
    }


}
