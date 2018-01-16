package org.jcvi.vigor.RegressionTest;

import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.jcvi.vigor.AppConfig;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.*;
import org.jcvi.vigor.utils.GenerateExonerateOutput;
import java.util.Comparator;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class GenerateVigor4GeneModels {

    public Map<String, List<Model>> generateModels(String workspace,
                                                   String inputFilePath, String refDB)
    {
        Map<String, List<Model>> vigor4Models = new HashMap<String, List<Model>>();
        try{
            File file = new File(workspace);
            if (!file.isDirectory())
                file = file.getParentFile();
            if (file.exists()) {
                File inputFile = new File(inputFilePath);
                if (inputFile.exists()) {
                    AlignmentGenerationService alignmentGenerationService = new AlignmentGenerationService();
                    ExonerateService exonerateService = new ExonerateService();
                    ViralProteinService viralProteinService = new ViralProteinService();
                    ModelGenerationService modelGenerationService = new ModelGenerationService();
                    GeneModelGenerationService geneModelGenerationService = new GeneModelGenerationService();
                    CheckCoverage checkCoverage = new CheckCoverage();
                    EvaluateScores evaluateScores = new EvaluateScores();
                    NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(
                            inputFile).hint(
                            DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
                            .build();
                    Stream<NucleotideFastaRecord> records = dataStore.records();
                    Iterator<NucleotideFastaRecord> i = records.iterator();
                    while (i.hasNext()) {
                        NucleotideFastaRecord record = i.next();
                        VigorForm form = new VigorForm();
                        VirusGenome virusGenome = new VirusGenome(
                                record.getSequence(), record.getComment(),
                                record.getId(), false, false);
                        alignmentGenerationService.GenerateAlignment(virusGenome,new VigorForm());
                        String fileName = GenerateExonerateOutput.queryExonerate(
                                virusGenome, refDB, file.getAbsolutePath(),null);
                        File outputFile = new File(fileName);
                        AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
                        alignmentEvidence.setReference_db(refDB);
                        List<Alignment> alignments = exonerateService
                                .parseExonerateOutput(outputFile,
                                        alignmentEvidence, virusGenome);
                        alignments = alignments
                                .stream()
                                .map(alignment -> viralProteinService
                                        .setViralProteinAttributes(alignment,form))
                                .collect(Collectors.toList());
                        alignments =  modelGenerationService.mergeIdenticalProteinAlignments(alignments);
                        List<Model> candidateModels = modelGenerationService.determineCandidateModels(alignments,form);
                        List<Model> processedModels = geneModelGenerationService.determineGeneFeatures(candidateModels,form);
                        processedModels.stream().forEach(model-> checkCoverage.evaluate(model,form));
                        processedModels.stream().forEach(model-> evaluateScores.evaluate(model, form));

                        processedModels.sort(new Comparator<Model>(){
                            @Override
                            public int compare(Model m1, Model m2){
                                return Double.compare(m1.getScores().get("totalScore"), m2.getScores().get("totalScore"));
                            }
                        });
                        processedModels.forEach(System.out::println);
                        vigor4Models.put(virusGenome.getId(), candidateModels);
                    }

                } else {
                    System.out.println("InputFile does not exist");
                }
            }

            else {
                System.out.println("Workspace folder does not exit");
            }}
        catch(IOException e){
            System.out.println(e.getMessage());
        }
        catch(Exception e){
            System.out.println(e.getMessage());
        }

        return vigor4Models;

    }
}
