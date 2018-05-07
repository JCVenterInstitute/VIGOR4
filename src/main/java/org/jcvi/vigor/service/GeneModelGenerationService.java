package org.jcvi.vigor.service;
import java.text.Normalizer;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.service.exception.ServiceException;
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
		List<Model> partialGeneModels = new ArrayList<Model>();
		List<Model> pseudoGenes = new ArrayList<Model>();
		isDebug = form.isDebug();
		List<Model> processedModels = determineGeneFeatures(models, partialGeneModels, pseudoGenes, form);
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
        FormatVigorOutput.printSequenceFeatures(geneModels);
     	return geneModels;

	}
	public List<Model> filterGeneModels(List<Model> models){
       Map<String,List<Model>> groupedModels = models.stream().collect(groupingBy(model -> model.getAlignment().getViralProtein().getProteinID()));
       List<Model> filteredModels = new ArrayList<Model>();
       List<Model> geneModels = new ArrayList<Model>();
       Set<String> proteinIDs = groupedModels.keySet();
       for(String proteinID : proteinIDs){
           List<Model> similarModels = groupedModels.get(proteinID);
           similarModels.sort(new Comparator<Model>() {
               @Override
               public int compare(Model m1, Model m2) {
                   return Double.compare(m1.getScores().get("totalScore"), m2.getScores().get("totalScore"));
               }
           });
           filteredModels.add(similarModels.get(similarModels.size()-1));
       }
        filteredModels.sort(new Comparator<Model>() {
            @Override
            public int compare(Model m1, Model m2) {
                return Double.compare(m1.getScores().get("modelScore"), m2.getScores().get("modelScore"));
            }
        });
       if(isDebug) {
           FormatVigorOutput.printModels(filteredModels, "Candidate Models");
       }
       Model geneModel = filteredModels.get(filteredModels.size()-1);
       int x=1;
       String genomeID = geneModel.getAlignment().getVirusGenome().getId();
        String regex = Pattern.quote("|");
        String[] genomeIDParts = genomeID.split(regex);
        String geneID="";
        if(genomeIDParts.length>=2) {
            geneID = genomeIDParts[0] + "_" + genomeIDParts[1];
        }else{
            geneID=genomeIDParts[0];
        }
       geneModel.setGeneID(geneID+"."+x);
       geneModels.add(geneModel);
       filteredModels.remove(geneModel);
       List<String> sharedCDSList = geneModel.getAlignment().getViralProtein().getGeneAttributes().getStructuralSpecifications().getShared_cds();
       if(! (sharedCDSList ==null || filteredModels.isEmpty()) ){
           for(String sharedCDS : sharedCDSList){
               for(int j=filteredModels.size()-1;j >= 0; j--){
                 Model model = filteredModels.get(j);
                   if(model.getGeneSymbol().equals(sharedCDS)){
                       x++;
                       model.setGeneID(geneID+"."+x);
                       geneModels.add(model);
                       break;
                   }
               }
           }
       }

       return geneModels;
    }

	public List<Model> determineGeneFeatures(List<Model> models, List<Model> partialGeneModels, List<Model> pseudoGenes, VigorForm form) throws ServiceException {

		List<Model> modelsWithMissingExonsDetermined=new ArrayList<Model>();
		List<Model> modelsAfterDeterminingStart =new ArrayList<Model>();
	 	List<Model> modelsAfterDeterminingViralTricks=new ArrayList<Model>();
	 	List<Model> modelsAfterAdjustingBounds = new ArrayList<Model>();
	 	List<Model> modelsAfterDeterminingStop = new ArrayList<Model>();

	  				
		/* Determine Start */
		for (Model model: models) {
			if (!model.isPartial5p()) {
				List<Model> outputModels = determineStart.determine(model, form);
				outputModels.stream().forEach(model1 -> {
					if ((model.isPartial5p())) {
						partialGeneModels.add(model);
					} else if (model.isPseudogene()) {
						pseudoGenes.add(model1);
					} else {
						modelsAfterDeterminingStart.add(model1);
					}

				});
			} else {
				partialGeneModels.add(model);
			}
		};

		if (isDebug) {
			FormatVigorOutput.printModelsWithStart(modelsAfterDeterminingStart, "After Determining Start");
		}
		
		/*Adjust RNAEditing, Ribosomal Slippage and find StopCodonReadThrough*/
		for (Model model: modelsAfterDeterminingStart) {
			List<Model> outputModels = adjustViralTricks.determine(model, form);
			modelsAfterDeterminingViralTricks.addAll(outputModels);
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
			if(!model.isPartial3p()){
				List<Model> outputModels = determineStop.determine(model, form);
				outputModels.stream().forEach(m->{
					if(m.isPseudogene()) {
						pseudoGenes.add(m);
					}else if(m.isPartial3p()){
						partialGeneModels.add(m);
					}else{
						modelsAfterDeterminingStop.add(m);
					}
				});
			}else{
				partialGeneModels.add(model);
			}
		}
		if(isDebug) {
			FormatVigorOutput.printModels(modelsAfterDeterminingStop,"Models after determining stop");
		}
		return modelsAfterDeterminingStop;
	}
	
	
	
	}

