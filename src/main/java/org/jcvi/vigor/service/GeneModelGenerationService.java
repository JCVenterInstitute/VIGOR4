package org.jcvi.vigor.service;

import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Pseudogene;
import org.jcvi.vigor.component.StructuralSpecifications;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.*;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;

import static java.util.stream.Collectors.groupingBy;

@Service
public class GeneModelGenerationService {

    @Autowired
    private DetermineMissingExons determineMissingExons;
    @Autowired
    private DetermineStart determineStart;
    @Autowired
    private AdjustViralTricks adjustViralTricks;
    @Autowired
    private AdjustUneditedExonBoundaries adjustUneditedExonBounds;
    @Autowired
    private DetermineStop determineStop;
    @Autowired
    private CheckCoverage checkCoverage;
    @Autowired
    private EvaluateScores evaluateScores;

    private static final Logger LOGGER = LogManager.getLogger(GeneModelGenerationService.class);

    public List<Model> generateGeneModel ( List<Model> models, VigorForm form ) throws ServiceException {

        List<Model> pseudoGenes = new ArrayList<>();
        boolean isDebug = form.getConfiguration().getOrDefault(ConfigurationParameters.Verbose, false);
        int max_gene_overlap = form.getConfiguration().getOrDefault(ConfigurationParameters.MaxGeneOverlap, 0);

        List<Model> processedModels = determineGeneFeatures(models, form, isDebug);
        // TODO process pseudogenes/partial genes, Not included in initial release
        processedModels.stream().forEach(model -> checkCoverage.evaluate(model, form));
        processedModels.stream().forEach(model -> evaluateScores.evaluate(model, form));
        processedModels.stream().forEach(model -> {
            if (model.isPseudogene()) {
                pseudoGenes.add(model);
            }
        });
        processedModels.removeAll(pseudoGenes);
        List<Model> processedPseudoGenes = processPseudogenes(pseudoGenes);
        if (processedPseudoGenes.size() > 0 && isDebug) {
            FormatVigorOutput.printAllGeneModelsWithScores(processedPseudoGenes, "Pseudogenes");
        }
        processedModels = filterModelsWithStructuralSpecifications(processedModels, form.getConfiguration());
        processedPseudoGenes = filterModelsWithStructuralSpecifications(processedPseudoGenes, form.getConfiguration());
        if (processedModels.size() <= 0) {
            LOGGER.error("No gene models found. Currently Vigor4 does not support annotating Pseudogenes ");
            return Collections.EMPTY_LIST;
        }
        List<Model> geneModels = filterGeneModels(processedModels, processedPseudoGenes, max_gene_overlap, isDebug);
        return geneModels;
    }

    /**
     * @param pseudogenes
     * @return
     */
    private List<Model> processPseudogenes ( List<Model> pseudogenes ) {

        pseudogenes = filterModelsOfaGene(pseudogenes);
        pseudogenes = filterUnOverlappedCandidateModels(pseudogenes);
        return pseudogenes;
    }

    /**
     * @param models
     * @param scoreType
     * @return
     */
    private List<Model> sortModels ( List<Model> models, String scoreType ) {

        if (models.size() > 1) {
            Collections.sort(models, ( m1, m2 ) -> {
                if (Double.compare(m2.getScores().get(scoreType), m1.getScores().get(scoreType)) == 0) {
                    return Double.compare(m1.getExons().size(), m2.getExons().size());
                } else
                    return Double.compare(m2.getScores().get(scoreType), m1.getScores().get(scoreType));
            });
        }
        return models;
    }

    /**
     * @param m1
     * @param m2
     * @return
     */
    private boolean checkExonOverlap ( Model m1, Model m2 ) {

        for (Exon M1exon : m1.getExons()) {
            for (Exon M2exon : m2.getExons()) {
                Range intersection = M1exon.getRange().intersection(M2exon.getRange());
                if (intersection.getLength() != 0) return true;
            }
        }
        return false;
    }

    /**
     * @param models
     * @param model
     * @return
     */
    private boolean isUnoverlappedCandidateModel ( List<Model> models, Model model ) {

        boolean overlap = false;
        for (Model model1 : models) {
            if (checkExonOverlap(model1, model)) {
                overlap = true;
            }
        }
        if (overlap) {
            return false;
        } else return true;
    }

    /**
     * @param models
     * @return
     */
    private double calculateTotalModelsScore ( List<Model> models ) {

        double totalScore = 0;
        for (Model model : models) {
            double score = model.getScores().get("modelScore");
            totalScore = totalScore + score;
        }
        return totalScore;
    }

    /**
     * @param models
     * @return
     */
    private List<Model> filterModelsOfaGene ( List<Model> models ) {

        Map<String, List<Model>> groupedModels = models.stream().collect(groupingBy(model -> model.getGeneSymbol()));
        List<Model> filteredModels = new ArrayList<Model>();
        Set<String> geneIDs = groupedModels.keySet();
        for (String geneID : geneIDs) {
            List<Model> similarModels = groupedModels.get(geneID);
            similarModels = sortModels(similarModels, "totalScore");
            List<Model> unOverlappedModels = new ArrayList<>();
            Model highScoringModel = similarModels.get(0);
            unOverlappedModels.add(highScoringModel);
            similarModels.remove(highScoringModel);
            for (Model checkModel : similarModels) {
                if (isUnoverlappedCandidateModel(unOverlappedModels, checkModel)) unOverlappedModels.add(checkModel);
            }
            filteredModels.addAll(unOverlappedModels);
        }
        return filteredModels;
    }

    /**
     * @param models
     * @return
     */
    private List<Model> filterUnOverlappedCandidateModels ( List<Model> models ) {

        List<Model> candidateGenes = new ArrayList<>();
        Map<String, List<Model>> genewiseModels = models.stream().collect(Collectors.groupingBy(x -> x.getGeneSymbol()));
        for (Map.Entry<String, List<Model>> entry : genewiseModels.entrySet()) {
            List<Model> aGeneModels = entry.getValue();
            Map<String, List<Model>> aGeneModelsProteinWise = aGeneModels.stream().collect(Collectors.groupingBy(x -> x.getAlignment().getViralProtein().getProteinID()));
            double highScore = 0;
            List<Model> highScoredModels = new ArrayList<>();
            for (Map.Entry<String, List<Model>> aProteinModelsOfGene : aGeneModelsProteinWise.entrySet()) {
                List<Model> aProteinModels = aProteinModelsOfGene.getValue();
                double totalModelsScore = calculateTotalModelsScore(aProteinModels);
                if (totalModelsScore > highScore) {
                    highScoredModels = new ArrayList<>();
                    highScoredModels.addAll(aProteinModels);
                    highScore = totalModelsScore;
                }
            }
            if (highScoredModels.size() > 1) {
                highScoredModels = highScoredModels.stream()
                        .sorted(Comparator.comparing(g -> g.getRange(), Range.Comparators.ARRIVAL))
                        .collect(Collectors.toList());
                char alphabet = 'a';
                for (Model fragmentModel : highScoredModels) {
                    fragmentModel.setGeneID(Character.toString(alphabet));
                    alphabet++;
                }
            }
            candidateGenes.addAll(highScoredModels);
            aGeneModels.removeAll(highScoredModels);
            for (Model tempModel : aGeneModels) {
                if (isUnoverlappedCandidateModel(aGeneModels, tempModel)) {
                    candidateGenes.add(tempModel);
                }
            }
        }
        return candidateGenes;
    }

    /**
     * @param models
     * @return
     */
    private List<Model> filterModelsWithStructuralSpecifications ( List<Model> models, VigorConfiguration config ) {

        List<Model> filteredModels = new ArrayList<>();
        int min_gene_size = config.getOrDefault(ConfigurationParameters.GeneMinimumSize,0);
        int min_coverage = config.getOrDefault(ConfigurationParameters.GeneMinimumCoverage, 0);
        for (Model model : models) {
            StructuralSpecifications specs = model.getAlignment().getViralProtein().getGeneAttributes().getStructuralSpecifications();
            int minFunLength = specs.getMinFunctionalLength();
            if (( model.getTranslatedSeq().getLength() >= minFunLength ) && ( model.getScores().get("%coverage") >= min_coverage ) &&
                    ( model.getTranslatedSeq().getLength() >= min_gene_size )) {
                filteredModels.add(model);
            }
        }
        return filteredModels;
    }

    /**
     * @param models
     * @param pseudogenes
     * @return
     */
    private List<Model> filterGeneModels ( List<Model> models, List<Model> pseudogenes, int max_gene_overlap, boolean isDebug ) {

        List<Model> geneModels = new ArrayList<>();
        List<Model> filteredModels = filterModelsOfaGene(models);
        List<Model> candidateGenes = filterUnOverlappedCandidateModels(filteredModels);
        if (isDebug) {
            FormatVigorOutput.printAllGeneModelsWithScores(candidateGenes, "All Gene Models");
        }
        candidateGenes = sortModels(candidateGenes, "modelScore");
        Model geneModel = candidateGenes.get(0);
        candidateGenes.remove(geneModel);
        String genomeID = geneModel.getAlignment().getVirusGenome().getId();
        String[] genomeIDParts = genomeID.split(Pattern.quote("|"));
        String proteinIDOfGenome;
        if (genomeIDParts.length >= 2) {
            proteinIDOfGenome = genomeIDParts[ 0 ] + "_" + genomeIDParts[ 1 ];
        } else {
            proteinIDOfGenome = genomeIDParts[ 0 ];
        }
        IDGenerator idGenerator = IDGenerator.of(proteinIDOfGenome);
        geneModels.add(geneModel);
        List<String> sharedCDSList = geneModel.getAlignment().getViralProtein().getGeneAttributes().getStructuralSpecifications().getShared_cds();
        for (int j = candidateGenes.size() - 1; j >= 0; j--) {
            Model candidateGene = candidateGenes.get(j);
            List<Exon> candidateExons = candidateGene.getExons();
            boolean overlap = false;
            for (Model model : geneModels) {
                List<Exon> exons = model.getExons();
                for (Exon candidateExon : candidateExons) {
                    for (Exon exon : exons) {
                        Range intersection = candidateExon.getRange().intersection(exon.getRange());
                        if (intersection.getLength() > max_gene_overlap) {
                            overlap = true;
                        }
                    }
                }
            }
            if (!overlap) geneModels.add(candidateGene);
        }
        // add shared_cds model to gene models
        candidateGenes.removeAll(geneModels);
        if (sharedCDSList != null) {
            for (String sharedCDS : sharedCDSList) {
                boolean found = false;
                for (Model model : candidateGenes) {
                    if (model.getGeneSymbol().equals(sharedCDS)) {
                        geneModels.add(model);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    for (Model pseudoModel : pseudogenes) {
                        if (pseudoModel.getGeneSymbol().equals(sharedCDS)) {
                            geneModels.add(pseudoModel);
                            break;
                        }
                    }
                }
            }
        }
        geneModels = geneModels.stream()
                .sorted(Comparator.comparing(g -> g.getRange(), Range.Comparators.ARRIVAL))
                .collect(Collectors.toList());
        String proteinID = "";
        String id = "";
        for (Model tempGeneModel : geneModels) {
            if (proteinID.equals(tempGeneModel.getAlignment().getViralProtein().getProteinID())) {
                tempGeneModel.setGeneID(id + tempGeneModel.getGeneID());
            } else {
                id = idGenerator.next();
                if (tempGeneModel.getGeneID() != null) {
                    tempGeneModel.setGeneID(id + tempGeneModel.getGeneID());
                } else tempGeneModel.setGeneID(id);
                proteinID = tempGeneModel.getAlignment().getViralProtein().getProteinID();
            }
        }
        return geneModels;
    }

    /**
     * @param models
     * @param form
     * @return
     * @throws ServiceException
     */
    private List<Model> determineGeneFeatures ( List<Model> models, VigorForm form, boolean isDebug ) throws ServiceException {

        List<Model> modelsWithMissingExonsDetermined = new ArrayList<Model>();
        List<Model> modelsAfterDeterminingStart = new ArrayList<Model>();
        List<Model> modelsAfterDeterminingViralTricks = new ArrayList<Model>();
        List<Model> modelsAfterAdjustingBounds = new ArrayList<Model>();
        List<Model> modelsAfterDeterminingStop = new ArrayList<Model>();
        if (isDebug) {
            FormatVigorOutput.printModels(models, "Models after processing fragments");
        }

        /* Determine Start */
        for (Model model : models) {
            if (!model.isPartial5p()) {
                List<Model> outputModels = determineStart.determine(model, form);
                outputModels.stream().forEach(model1 -> {
                    modelsAfterDeterminingStart.add(model1);
                });
            } else modelsAfterDeterminingStart.add(model);
        }
        ;
        if (isDebug) {
            FormatVigorOutput.printModels(modelsAfterDeterminingStart, "After Determining Start");
        }

        /*Adjust RNAEditing, Ribosomal Slippage and find StopCodonReadThrough*/
        for (Model model : modelsAfterDeterminingStart) {
            List<Model> outputModels = adjustViralTricks.determine(model, form);
            modelsAfterDeterminingViralTricks.addAll(outputModels);
        }
        if (isDebug) {
            FormatVigorOutput.printModels(modelsAfterDeterminingViralTricks, "After determining viral tricks");
        }

        /*Adjust unedited Exon boundaries*/
        for (Model model : modelsAfterDeterminingViralTricks) {
            List<Model> outputModels = adjustUneditedExonBounds.determine(model, form);
            for (Model adjustedModel : outputModels) {
                int exonsCount = adjustedModel.getExons().size();
                Model missingExonsDeterminedModel = determineMissingExons.determine(adjustedModel, form).get(0);
                int afterExonsCount = missingExonsDeterminedModel.getExons().size();
                if (afterExonsCount > exonsCount) {
                    List<Model> revisitedModels = adjustUneditedExonBounds.determine(missingExonsDeterminedModel, form);
                    modelsWithMissingExonsDetermined.addAll(revisitedModels);
                } else {
                    modelsWithMissingExonsDetermined.add(adjustedModel);
                }
            }
            modelsAfterAdjustingBounds.addAll(outputModels);
        }
        if (isDebug) {
            FormatVigorOutput.printModels(modelsWithMissingExonsDetermined, "After determining missing exons");
        }

        /* Determine Stop */
        for (Model model : modelsWithMissingExonsDetermined) {
            if (!model.isPartial3p()) {
                List<Model> outputModels = determineStop.determine(model, form);
                outputModels.stream().forEach(m -> {
                    modelsAfterDeterminingStop.add(m);
                });
            } else modelsAfterDeterminingStop.add(model);
        }
        if (isDebug) {
            FormatVigorOutput.printModels(modelsAfterDeterminingStop, "Models after determining stop");
        }
        return modelsAfterDeterminingStop;
    }
}

