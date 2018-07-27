package org.jcvi.vigor.utils;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Stream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.AlignmentToolFactory;
import org.jcvi.vigor.service.ExonerateService;
import org.jcvi.vigor.service.ViralProteinService;
import org.jcvi.vigor.service.VirusGenomeService;
import org.jcvi.vigor.service.exception.ServiceException;

import static org.junit.Assume.assumeTrue;

public class VigorTestUtils {

    private static final Logger LOGGER = LogManager.getLogger(VigorTestUtils.class);

    public static List<Alignment> getAlignments ( File inputSeqFile, String refDB, File alignmentOutput, VigorConfiguration config ) throws VigorException {

        ExonerateService exonerateService = new ExonerateService();
        ViralProteinService viralProteinService = new ViralProteinService();
        try (NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(inputSeqFile).hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED).build();) {
            Stream<NucleotideFastaRecord> records = dataStore.records();
            Iterator<NucleotideFastaRecord> iter = records.iterator();
            NucleotideFastaRecord record = iter.next();
            VirusGenome virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(), false,
                    false);
            List<Range> sequenceGaps = VirusGenomeService.findSequenceGapRanges(20, virusGenome.getSequence());
            virusGenome.setSequenceGaps(sequenceGaps);
            // create alignment evidence
            AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
            VigorForm form = new VigorForm(config);
            alignmentEvidence.setReference_db(refDB);
            form.setAlignmentEvidence(alignmentEvidence);
            List<Alignment> alignments = exonerateService.parseExonerateOutput(alignmentOutput,
                    form, virusGenome, refDB);
            for (int i = 0; i < alignments.size(); i++) {
                alignments.set(i, viralProteinService
                        .setViralProteinAttributes(alignments.get(i), new VigorForm(config)));
            }
            return alignments;
        } catch (IOException e) {
            throw new ServiceException(String.format("Problem reading fasta file %s", inputSeqFile), e);
        }
    }

    public static void assumeReferenceDB ( String referenceDBPath ) {

        String skipReferenceDBTests = System.getProperty("vigor.skip_reference_db_tests");
        LOGGER.info("skip reference db tests '{}'", skipReferenceDBTests);
        LOGGER.info("reference database path is '{}'", referenceDBPath);
        if (!( skipReferenceDBTests == null || skipReferenceDBTests.isEmpty() )) {
            LOGGER.info("skipping test if reference database pathis not set");
            assumeTrue("Skipping test requiring reference database path", !( referenceDBPath == null || referenceDBPath.isEmpty() ));
        }
    }
}
