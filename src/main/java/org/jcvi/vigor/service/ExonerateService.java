package org.jcvi.vigor.service;

import org.jcvi.jillion.core.Direction;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.align.exonerate.Exonerate2;
import org.jcvi.jillion.align.exonerate.vulgar.VulgarProtein2Genome2;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by snettem on 5/17/2017.
 */
@Service
public class ExonerateService implements AlignmentService {

    private static final Logger LOGGER = LogManager.getLogger(ExonerateService.class);
    private final AlignmentTool alignmentTool = new Exonerate("exonerate");

     /**
     * @param config
     * @param virusGenome
     * @param referenceDB
     * @param workspace
     * @return
     * @throws ServiceException
     */
    @Override
    public List<Alignment> getAlignment (VigorConfiguration config, VirusGenome virusGenome, String referenceDB, String workspace ) throws ServiceException {

        try {
            String exoneratePathString = config.get(ConfigurationParameters.ExoneratePath);

            LOGGER.debug("Using exonerate path {}", exoneratePathString);

            VigorUtils.checkFilePath("exonerate path via config value " + ConfigurationParameters.ExoneratePath.configKey,
                                     exoneratePathString, VigorUtils.FileCheck.EXISTS, VigorUtils.FileCheck.EXECUTE);
            Path exoneratePath = Paths.get(exoneratePathString);
            String outputFilePath = GenerateExonerateOutput.queryExonerate(virusGenome, referenceDB, workspace, null, exoneratePath.toString());
            File outputFile = new File(outputFilePath);
            return parseExonerateOutput(outputFile, virusGenome, referenceDB,config);
        } catch (VigorException e) {
            throw new ServiceException(String.format("error getting alignment got %s: %s", e.getClass().getSimpleName(), e.getMessage()), e);
        }
    }

    @Override
    public AlignmentTool getAlignmentTool() {
        return alignmentTool;
    }

    /**
     * @param exonerateOutput
     * @param virusGenome
     * @param referenceDB
     * @return
     * @throws ServiceException
     */
    public List<Alignment> parseExonerateOutput ( File exonerateOutput, VirusGenome virusGenome, String referenceDB, VigorConfiguration config ) throws ServiceException {

        List<Alignment> alignments = new ArrayList<Alignment>();
        List<VulgarProtein2Genome2> Jalignments;
        long seqLength = virusGenome.getSequence().getLength();
        VirusGenome reverseCompGenome = VirusGenomeService.reverseComplementVirusGenome(virusGenome,config);
        AlignmentEvidence evidence = new AlignmentEvidence();
        evidence.setReference_db(referenceDB);
        // TODO results directory
        evidence.setRaw_alignment(exonerateOutput);
        evidence.setResults_directory(exonerateOutput.getParentFile());
        AlignmentTool alignmentTool = getAlignmentTool();
        try {
            Jalignments = Exonerate2.parseVulgarOutput(exonerateOutput);
        } catch (IOException e) {
            throw new ServiceException(String.format("Error parsing exonerate output %s", exonerateOutput.getName()));
        }
        try (ProteinFastaDataStore datastore = new ProteinFastaFileDataStoreBuilder(new File(referenceDB))
                .hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED).build();
        ) {
            for (VulgarProtein2Genome2 Jalignment : Jalignments) {
                Alignment alignment = new Alignment();
                Map<String, Double> alignmentScores = new HashMap<String, Double>();
                alignmentScores.put("alignmentScore", (double) Jalignment.getScore());
                alignment.setAlignmentScore(alignmentScores);
                alignment.setAlignmentTool(alignmentTool);
                List<AlignmentFragment> alignmentFragments = new ArrayList<>();
                for (VulgarProtein2Genome2.AlignmentFragment fragment : Jalignment.getAlignmentFragments()) {
                    if(fragment.getDirection().equals(Direction.FORWARD)) {
                        alignmentFragments.add(new AlignmentFragment(fragment.getProteinSeqRange(),
                                fragment.getNucleotideSeqRange(),
                                fragment.getDirection(),
                                fragment.getFrame()));
                    }else {

                            alignmentFragments.add(new AlignmentFragment(fragment.getProteinSeqRange(),
                                    VigorFunctionalUtils.convertToOppositeStrandRange(fragment.getNucleotideSeqRange(),seqLength),
                                    fragment.getDirection(),
                                    fragment.getFrame()));

                    }
                }
                ProteinFastaRecord fasta = datastore.get(Jalignment.getQueryId());
                ViralProtein viralProtein = new ViralProtein();
                viralProtein.setProteinID(fasta.getId());
                viralProtein.setDefline(fasta.getComment());
                viralProtein.setSequence(fasta.getSequence());
                alignment.setAlignmentFragments(alignmentFragments);
                alignment.setViralProtein(viralProtein);
                if(Jalignment.getTargetStrand().get().equals(Direction.REVERSE)){
                    alignment.setVirusGenome(reverseCompGenome);
                }else alignment.setVirusGenome(virusGenome);
                alignment.setAlignmentEvidence(evidence.copy());
                alignments.add(alignment);
            }
        } catch (IOException e) {
            LOGGER.error(String.format("Problem reading virus database file %s", referenceDB), e);
            throw new ServiceException(e);
        }
        return alignments;
    }
}
