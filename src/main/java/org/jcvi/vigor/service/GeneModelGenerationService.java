package org.jcvi.vigor.service;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.GeneModel;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.IDGenerator;
import org.jcvi.vigor.utils.NoteType;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.FormatVigorOutput;

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

	boolean isDebug = false;
    private static final Logger LOGGER = LogManager.getLogger(GeneModelGenerationService.class);

    public List<Model> generateGeneModel(List<Model> models, VigorForm form) throws ServiceException {
		List<Model> pseudoGenes = new ArrayList<Model>();
		isDebug = form.getConfiguration().get(ConfigurationParameters.Verbose).equals("true") ? true : false;
		List<Model> processedModels = determineGeneFeatures(models, pseudoGenes, form);
		// TODO process pseudogenes/partial genes, Not included in initial release

		processedModels.stream().forEach(model -> checkCoverage.evaluate(model, form));
		processedModels.stream().forEach(model -> {
			if (model.isPseudogene()) {
				pseudoGenes.add(model);
			}
		});
		processedModels.removeAll(pseudoGenes);
		processedModels.stream().forEach(model -> evaluateScores.evaluate(model, form));
		if(processedModels.size()<=0){
            LOGGER.error("No gene models found. Currently Vigor4 does not support annotating Partial or Pseudogenes ");
			return Collections.EMPTY_LIST;
        }
        List<Model> geneModels = filterGeneModels(processedModels);
		return geneModels;

	}
	private List<Model> sortModels(List<Model> models,String scoreType){
        Collections.sort(models, new Comparator<Model>() {
            @Override
            public int compare(Model m1, Model m2) {
                return Double.compare(m2.getScores().get(scoreType), m1.getScores().get(scoreType));
            }
        });
        return models;
    }
    private boolean checkExonOverlap (Model m1, Model m2){

         for(Exon M1exon : m1.getExons()){
             for(Exon M2exon : m2.getExons()){
                        if(M1exon.getRange().intersects(M2exon.getRange())) return true;
                    }
                }

        return false;
    }
    private boolean isUnoverlappedCandidateModel(List<Model> models , Model model) {
        for (Model model1 : models) {
            if (!checkExonOverlap(model1, model)) {
                return true;
            }
        }
        return false;
    }

  	public List<Model> filterGeneModels(List<Model> models){
       Map<String,List<Model>> groupedModels = models.stream().collect(groupingBy(model -> model.getAlignment().getViralProtein().getProteinID()));
       List<Model> filteredModels = new ArrayList<Model>();
       List<Model> geneModels = new ArrayList<Model>();
       Set<String> proteinIDs = groupedModels.keySet();
       for(String proteinID : proteinIDs){
           List<Model> similarModels = groupedModels.get(proteinID);
           similarModels = sortModels(similarModels,"totalScore");
           List<Model> unOverlappedModels=new ArrayList<Model>();
           Model highScoringModel = similarModels.get(0);
           unOverlappedModels.add(highScoringModel);
           similarModels.remove(highScoringModel);
           for(Model checkModel : similarModels) {
             if(isUnoverlappedCandidateModel(unOverlappedModels,checkModel)) unOverlappedModels.add(checkModel);
           }
           filteredModels.addAll(unOverlappedModels);
       }

       List<Model> candidateGenes = new ArrayList<Model>();
        if(isDebug) {
            FormatVigorOutput.printAllGeneModelsWithScores(filteredModels,"All Gene Models");
        }
       Map<String,List<Model>> genewiseModels = filteredModels.stream().collect(Collectors.groupingBy(x -> x.getGeneSymbol()));
       for(Map.Entry<String,List<Model>> entry : genewiseModels.entrySet()){
           List<Model> tempModels = entry.getValue();
           tempModels=sortModels(tempModels,"modelScore");
           List<Model> aGeneModels = new ArrayList<Model>();
           aGeneModels.add(tempModels.get(0));
           tempModels.remove(0);
           for(Model tempModel : tempModels){
               if (isUnoverlappedCandidateModel(aGeneModels, tempModel)) {
                   aGeneModels.add(tempModel);
               }
           }
           candidateGenes.addAll(aGeneModels);
       }
        if(isDebug) {
            FormatVigorOutput.printAllGeneModelsWithScores(candidateGenes,"Candidate Gene Models");
        }
        candidateGenes=sortModels(candidateGenes,"modelScore");
        Model geneModel = candidateGenes.get(0);
        candidateGenes.remove(geneModel);
        String genomeID = geneModel.getAlignment().getVirusGenome().getId();
        String[] genomeIDParts = genomeID.split(Pattern.quote("|"));
        String proteinIDOfGenome;
        if (genomeIDParts.length >= 2) {
            proteinIDOfGenome = genomeIDParts[0] + "_" + genomeIDParts[1];
        } else {
            proteinIDOfGenome = genomeIDParts[0];
        }

        IDGenerator idGenerator = IDGenerator.of(proteinIDOfGenome);
        geneModels.add(geneModel);
        List<String> sharedCDSList = geneModel.getAlignment().getViralProtein().getGeneAttributes().getStructuralSpecifications().getShared_cds();
        for(int j=candidateGenes.size()-1;j>=0;j--){
            Model candidateGene = candidateGenes.get(j);
            List<Exon> candidateExons = candidateGene.getExons();
            boolean overlap = false;
            for(Model model : geneModels){
                List<Exon> exons = model.getExons();
                for(Exon candidateExon : candidateExons){
                    for(Exon exon : exons){
                        if(candidateExon.getRange().intersects(exon.getRange())){
                            overlap=true;
                        }
                    }
                }
            }
            if(overlap && sharedCDSList!=null){
                for(String sharedCDS : sharedCDSList){
                    if(candidateGene.getGeneSymbol().equals(sharedCDS)){
                      geneModels.add(candidateGene);
                    }
                }
            }else if(!overlap){
               geneModels.add(candidateGene);
            }

        }
        geneModels = geneModels.stream()
                .sorted(Comparator.comparing(g -> g.getRange(), Range.Comparators.ARRIVAL))
                .collect(Collectors.toList());
        geneModels.stream().forEach(x->x.setGeneID(idGenerator.next()));

       return geneModels;
    }

	public List<Model> determineGeneFeatures(List<Model> models, List<Model> pseudoGenes, VigorForm form) throws ServiceException {

		List<Model> modelsWithMissingExonsDetermined=new ArrayList<Model>();
		List<Model> modelsAfterDeterminingStart =new ArrayList<Model>();
	 	List<Model> modelsAfterDeterminingViralTricks=new ArrayList<Model>();
	 	List<Model> modelsAfterAdjustingBounds = new ArrayList<Model>();
	 	List<Model> modelsAfterDeterminingStop = new ArrayList<Model>();

        if (isDebug) {
            FormatVigorOutput.printModels(models, "Candidate models after processing fragments");
        }
		/* Determine Start */
		for (Model model: models) {
		    if(!model.isPartial5p()) {
                List<Model> outputModels = determineStart.determine(model, form);
                outputModels.stream().forEach(model1 -> {
                    if (model.isPseudogene()) {
                        pseudoGenes.add(model1);
                    } else {
                        modelsAfterDeterminingStart.add(model1);
                    }

                });
            }else modelsAfterDeterminingStart.add(model);
		};

		if (isDebug) {
			FormatVigorOutput.printModels(modelsAfterDeterminingStart, "After Determining Start");
		}
		
		/*Adjust RNAEditing, Ribosomal Slippage and find StopCodonReadThrough*/
		for (Model model: modelsAfterDeterminingStart) {
			List<Model> outputModels = adjustViralTricks.determine(model, form);
			modelsAfterDeterminingViralTricks.addAll(outputModels);
		}

        if(isDebug){
            FormatVigorOutput.printModels(modelsAfterDeterminingViralTricks,"After determining viral tricks");
        }
		
		/*Adjust unedited Exon boundaries*/
		for (Model model: modelsAfterDeterminingViralTricks) {
			List<Model> outputModels = adjustUneditedExonBounds.determine(model, form);
			for(Model adjustedModel : outputModels){
				int exonsCount = adjustedModel.getExons().size();
				Model missingExonsDeterminedModel = determineMissingExons.determine(adjustedModel, form).get(0);
				int afterExonsCount = missingExonsDeterminedModel.getExons().size();
				if(afterExonsCount>exonsCount){
					List<Model> revisitedModels = adjustUneditedExonBounds.determine(missingExonsDeterminedModel, form);
					modelsWithMissingExonsDetermined.addAll(revisitedModels);
				}else{
					modelsWithMissingExonsDetermined.add(adjustedModel);
				}
			}
			modelsAfterAdjustingBounds.addAll(outputModels);
		}
		
		if(isDebug){
		    FormatVigorOutput.printModels(modelsWithMissingExonsDetermined,"After determining missing exons");
        }
		
		
		/* Determine Stop */
		for (Model model: modelsWithMissingExonsDetermined) {
		    if(!model.isPartial3p()) {
                List<Model> outputModels = determineStop.determine(model, form);
                outputModels.stream().forEach(m -> {
                    if (m.isPseudogene()) {
                        pseudoGenes.add(m);
                    } else {
                        modelsAfterDeterminingStop.add(m);
                    }
                });
            }else modelsAfterDeterminingStop.add(model);

		}
		if(isDebug) {
			FormatVigorOutput.printModels(modelsAfterDeterminingStop,"Models after determining stop");
		}
		return modelsAfterDeterminingStop;
	}
	
	
	
	}

