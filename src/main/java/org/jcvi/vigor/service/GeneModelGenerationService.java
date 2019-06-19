package org.jcvi.vigor.service;

import java.util.*;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.*;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.component.Model;

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

    private static Function<Model,Integer> getMinFunctionalLength = (model) -> model.getAlignment()
                                                                                    .getViralProtein()
                                                                                    .getGeneAttributes()
                                                                                    .getStructuralSpecifications()
                                                                                    .getMinFunctionalLength();

    public List<Model> generateGeneModel ( List<Model> models, VigorConfiguration configuration ) throws ServiceException {

        List<Model> pseudoGenes = new ArrayList<>();
        boolean isDebug = configuration.getOrDefault(ConfigurationParameters.Verbose, false);
        int max_gene_overlap = configuration.getOrDefault(ConfigurationParameters.MaxGeneOverlap, 0);

        List<Model> processedModels = determineGeneFeatures(models, configuration, isDebug);
        // TODO process pseudogenes/partial genes, Not included in initial release
        processedModels.stream().forEach(model -> checkCoverage.evaluate(model, configuration));
        processedModels.stream().forEach(model -> evaluateScores.evaluate(model, configuration));
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
        processedModels = filterModelsWithStructuralSpecifications(processedModels, configuration);
        processedPseudoGenes = filterModelsWithStructuralSpecifications(processedPseudoGenes, configuration);
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
                break;
            }
        }
        return ! overlap;
    }

    /**
     * @param models
     * @return
     */
    private double calculateTotalModelsScore ( List<Model> models ) {

        double totalScore = 0;
        for (Model model : models) {
            double score = model.getScores().get(Scores.MODEL_SCORE);
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
        for (String geneID : groupedModels.keySet()) {
            List<Model> similarModels = groupedModels.get(geneID);
            similarModels = sortModels(similarModels, Scores.TOTAL_SCORE);
            List<Model> unOverlappedModels = new ArrayList<>();
            unOverlappedModels.add(similarModels.remove(0));
            for (Model checkModel : similarModels) {
                if (isUnoverlappedCandidateModel(unOverlappedModels, checkModel)) {
                    unOverlappedModels.add(checkModel);
                }
            }
            filteredModels.addAll(unOverlappedModels);
        }
        return filteredModels;
    }

    /**
     * @param models
     * @return models of different genes where models belonging to a gene should not overlap(choose the best one if there is overlap)
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
     * @return models satisfying min gene size, coverage and functional length
     */
    private List<Model> filterModelsWithStructuralSpecifications(List<Model> models, VigorConfiguration config) {

        int min_coverage = config.getOrDefault(ConfigurationParameters.GeneMinimumCoverage, 0);

        return models.stream()
                     .filter(
                             m -> m.getTranslatedSeq().getLength() >= getMinFunctionalLength.apply(m) &&
                                     m.getScores().get(Scores.COVERAGE_SCORE) >= min_coverage
                     )
                     .collect(Collectors.toList());
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
        candidateGenes = sortModels(candidateGenes, Scores.MODEL_SCORE);
        Model geneModel = candidateGenes.remove(0);
        geneModels.add(geneModel);
        //generate geneID
        String genomeID = geneModel.getAlignment().getVirusGenome().getId();
        String[] genomeIDParts = genomeID.split(Pattern.quote("|"));
        String proteinIDOfGenome;
        if (genomeIDParts.length >= 2) {
            proteinIDOfGenome = genomeIDParts[ 0 ] + "_" + genomeIDParts[ 1 ];
        } else {
            proteinIDOfGenome = genomeIDParts[ 0 ];
        }
        IDGenerator idGenerator = IDGenerator.of(proteinIDOfGenome);
        ViralProtein protein = geneModel.getAlignment().getViralProtein();
        List<String> sharedCDSList = NullUtil.nullOrElse(protein.getGeneAttributes().getStructuralSpecifications().getShared_cds(), Collections.EMPTY_LIST);
        LOGGER.debug("for proteinID {}, gene {} checking for candidate models with shared cds: {}",
                     protein.getProteinID(),
                     protein.getGeneSymbol(),
                     String.join(",",sharedCDSList));

        LOGGER.debug("checking {} models for non overlapping models", candidateGenes.size());

        for (int j = candidateGenes.size() - 1; j >= 0; j--) {
            Model candidateGene = candidateGenes.get(j);
            List<Exon> candidateExons = candidateGene.getExons();
            boolean overlap = false;
            // Best model is picked as a gene model. Now compare other candidate models with genemodel. Check for overlap and do not add shared_CDS models at this step.
            CHECKOVERLAP:
            for (Model model : geneModels) {
                boolean isSharedCDS = false;
                List<String> tempSharedCDS = NullUtil.nullOrElse(model.getAlignment().getViralProtein().getGeneAttributes().getStructuralSpecifications().getShared_cds(), Collections.EMPTY_LIST);
               // below step is to not to add shared_CDS genes at this stage
                for(String sharedCDS:tempSharedCDS){
                    if(sharedCDS.equals(candidateGene.getGeneSymbol())) {
                        isSharedCDS = true;
                        overlap = true;
                    }
                }
                List<Exon> exons = model.getExons();
                if (!isSharedCDS && candidateGene.getDirection().equals(model.getDirection())){
                    for (Exon candidateExon : candidateExons) {
                        for (Exon exon : exons) {
                            Range intersection = candidateExon.getRange().intersection(exon.getRange());
                            if (intersection.getLength() > max_gene_overlap) {
                                overlap = true;
                                break CHECKOVERLAP;
                            }
                        }
                    }
            }
            }
            if (!overlap) {
                LOGGER.debug("for protein {}, gene {} adding non-overlapping candidate {}",
                             protein.getProteinID(),
                             protein.getGeneSymbol(),
                             candidateGene);
                geneModels.add(candidateGene);
                candidateGenes.remove(j);
            }
        }

        // add shared_cds model to gene models
        for (String sharedCDS : sharedCDSList) {
            Model sharedCDSModel = null;
            Optional<Model> model = candidateGenes.stream()
                    .filter(m -> m.getGeneSymbol().equals(sharedCDS))
                    .findAny();
            LOGGER.debug("for proteinID {} gene {} checking {} candidate models for shared CDs model {}: found {}",
                    protein.getProteinID(),
                    protein.getGeneSymbol(),
                    candidateGenes.size(),
                    sharedCDS,
                    model.isPresent());
            if (!geneModels.stream().filter(m -> m.getGeneSymbol().equals(sharedCDS)).findAny().isPresent()){
                if (model.isPresent()) {
                    sharedCDSModel = model.get();
                    geneModels.add(sharedCDSModel);
                } else {
                    model = pseudogenes.stream()
                            .filter(m -> m.getGeneSymbol().equals(( sharedCDS )))
                            .findAny();
                    LOGGER.debug("for proteinID {} gene {} checking {} pseudogenes for shared CDs model {}: found {}",
                            protein.getProteinID(),
                            protein.getGeneSymbol(),
                            pseudogenes.size(),
                            sharedCDS,
                            model.isPresent());
                    if (model.isPresent()) {
                        sharedCDSModel = model.get();
                        geneModels.add(sharedCDSModel);
                    }
                }
            }

        }
        geneModels.sort(Comparator.comparing(g -> g.getRange(), Range.Comparators.ARRIVAL));

        String proteinID = "";
        String id = "";
        for (Model tempGeneModel : geneModels) {
            if (! proteinID.equals(tempGeneModel.getAlignment().getViralProtein().getProteinID())) {
                id = idGenerator.next();
                proteinID = tempGeneModel.getAlignment().getViralProtein().getProteinID();
            }
            tempGeneModel.setGeneID(id + NullUtil.nullOrElse(tempGeneModel.getGeneID(),""));
        }
        return geneModels;
    }

    /**
     * @param models
     * @param configuration
     * @return
     * @throws ServiceException
     */
    private List<Model> determineGeneFeatures ( List<Model> models, VigorConfiguration configuration, boolean isDebug ) throws ServiceException {

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
                List<Model> outputModels = determineStart.determine(model);
                modelsAfterDeterminingStart.addAll(outputModels);
            } else {
                modelsAfterDeterminingStart.add(model);
            }
        }
        ;
        if (isDebug) {
            FormatVigorOutput.printModels(modelsAfterDeterminingStart, "After Determining Start");
        }

        /*Adjust RNAEditing, Ribosomal Slippage and find StopCodonReadThrough*/
        for (Model model : modelsAfterDeterminingStart) {
            List<Model> outputModels = adjustViralTricks.determine(model);
            modelsAfterDeterminingViralTricks.addAll(outputModels);
        }
        if (isDebug) {
            FormatVigorOutput.printModels(modelsAfterDeterminingViralTricks, "After determining viral tricks");
        }

        /*Adjust unedited Exon boundaries*/
        for (Model model : modelsAfterDeterminingViralTricks) {
            List<Model> outputModels = adjustUneditedExonBounds.determine(model);
            for (Model adjustedModel : outputModels) {
                int exonsCount = adjustedModel.getExons().size();
                Model missingExonsDeterminedModel = determineMissingExons.determine(adjustedModel).get(0);
                int afterExonsCount = missingExonsDeterminedModel.getExons().size();
                if (afterExonsCount > exonsCount) {
                    List<Model> revisitedModels = adjustUneditedExonBounds.determine(missingExonsDeterminedModel);
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
                List<Model> outputModels = determineStop.determine(model);
                modelsAfterDeterminingStop.addAll(outputModels);
            } else {
                modelsAfterDeterminingStop.add(model);
            }
        }
        if (isDebug) {
            FormatVigorOutput.printModels(modelsAfterDeterminingStop, "Models after determining stop");
        }
        return modelsAfterDeterminingStop;
    }
}

