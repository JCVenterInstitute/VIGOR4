package org.jcvi.vigor.RegressionTest;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.utils.TBLParser.TBLFileParser;
import org.jcvi.vigor.utils.TBLParser.TBLModel;

public class GenerateVigor3Models {

    private final static Logger LOGGER = LogManager.getLogger(GenerateVigor3Models.class);

    public Map<String, List<Model>> generateModels (String outputDirectory, String outputPrefix) throws IOException {
        String tblFile = Paths.get(outputDirectory, outputPrefix + ".tbl").toString();
        String pepFile = Paths.get(outputDirectory, outputPrefix + ".pep").toString();
        return generateModels(tblFile, pepFile, null);
    }

    public Map<String, List<Model>> generateModels ( String TBLFilePath, String PEPFilePath, String fastaFilePath ) throws IOException {

        Map<String, List<Model>> vigor3Models = new HashMap<String, List<Model>>();
        if (fastaFilePath != null) {
            NucleotideFastaDataStore datastore = new NucleotideFastaFileDataStoreBuilder(new File(fastaFilePath)).build();
            LOGGER.debug("Number of records in the fasta file are : " + datastore.getNumberOfRecords());
        }
        TBLFileParser TBLParser = new TBLFileParser();
        List<TBLModel> TBLModels = TBLParser.getModels(TBLFilePath, PEPFilePath);
        LOGGER.debug("Total Number of models are :" + TBLModels.size());
        for (TBLModel tblModel : TBLModels) {
            Model model = new Model();
            List<Model> models = new ArrayList<Model>();
            Alignment alignment = new Alignment();
            ViralProtein viralProtein = new ViralProtein();
            VirusGenome virusGenome = new VirusGenome();
            model.setGeneSymbol(tblModel.getGene());
            model.setGeneID(tblModel.getGeneID());
            viralProtein.setProduct(tblModel.getProduct());
            viralProtein.setProteinID(tblModel.getViralProteinID());
            alignment.setViralProtein(viralProtein);
            virusGenome.setId(tblModel.getVirusGenomeID());
            alignment.setVirusGenome(virusGenome);
            model.setAlignment(alignment);
            model.setExons(tblModel.getExons());
            model.setPartial5p(tblModel.is5Partial());
            model.setPartial3p(tblModel.is3Partial());
            model.setPseudogene(tblModel.isPseudoGene());
            if (tblModel.getStopCodonReadThrough() != null)
                model.setReplaceStopCodonRange(tblModel.getStopCodonReadThrough());
            if (tblModel.isRiboSlippage())
                model.setRibosomalSlippageRange(Range.of(0)); //since in TBL there is no specific row for range.we just capture if this feature exists for model
            String virusGenomeID = model.getAlignment().getVirusGenome().getId();
            if (vigor3Models.containsKey(virusGenomeID)) {
                models = vigor3Models.get(virusGenomeID);
            }
            models.add(model);
            vigor3Models.put(virusGenomeID, models);
        }
        LOGGER.debug(vigor3Models.entrySet().size());
        vigor3Models.entrySet().forEach(entry -> {
            LOGGER.trace("genome ID: \"{}\" model count: {}",entry.getKey(), entry.getValue().size());
        });
        return vigor3Models;
    }
}
