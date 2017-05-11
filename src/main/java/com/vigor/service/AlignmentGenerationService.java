package com.vigor.service;

import com.vigor.component.*;
import com.vigor.forms.VigorForm;
import org.apache.commons.lang3.StringUtils;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.springframework.stereotype.Service;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 *
 * Created by snettem on 5/9/2017.
 */
@Service
public class AlignmentGenerationService {

    private Alignment alignment;
    private ViralProtein viralProtein;

    public void GenerateAlignment(VirusGenome virusGenome , AlignmentEvidence alignmentEvidence, VigorForm form){
        String alignmentTool = chooseAlignmentTool (alignmentEvidence);

        if(alignmentTool.equals ( "jillion" ))
        {
           //Call Jillion helper method to retrieve the Jillion alignment object
            // create another method to set the alignment attributes received from Jillion object

        }

           setViralProteinAttributes ();


    }

    public void setViralProteinAttributes() {
        viralProtein = new ViralProtein ();
        viralProtein.setDefline ( "db=\"flua_db\" gene=\"M2\" intron_size=650-750 product=\"matrix protein 2\" splice_form=\"e26i690e264\" spliced=Y length=97" );
        ProteinSequence sequence = new ProteinSequenceBuilder ( "MSLLTEVETPTRNGWECKCSDSSDPLAIAANIIGILHLILWILDRLFFKCIYRRLKYGLK\n" +
                "RGPSTEGVPESMREEYRQEQQSAVDVDDSHFVNIELE\n" ).build ();
        viralProtein.setSequence ( sequence );
        viralProtein.setProteinID ( "seg7prot2ADR3_V27A&S31N" );
        GeneAttributes geneAttributes = new GeneAttributes ();
        String defline = viralProtein.getDefline ();
        defline = StringUtils.normalizeSpace ( defline );
        List<String> tok = Arrays.asList ( defline.split ( "[ ]+(?=([^\"]*\"[^\"]*\")*[^\"]*$)" ) );
        Map<String, String> attributes = Pattern.compile ( "[ ]+(?=([^\"]*\"[^\"]*\")*[^\"]*$)" )
                                                .splitAsStream ( defline.trim () )
                                                .map ( s -> s.split ( "=", 2 ) )
                                                .collect ( Collectors.toMap ( a -> a[0], a -> a.length > 1 ? a[1] : "" ) );

        /*Set Splicing attributes*/
        if (attributes.containsKey ( "splice_form" ) || attributes.containsKey ( "ribosomal_slippage" )) {
            Splicing splicing = new Splicing ();
            splicing.setSpliced ( true );
            splicing.setSpliceform ( (attributes.get ( "splice_form" )).replaceAll ( "^\"|\"$", "" ) );
            if (attributes.containsKey ( "noncanonical_splicing" )) {
                if (!(attributes.get ( "noncanonical_splicing" ).equalsIgnoreCase ( "N" ))) {
                    splicing.setNonCanonical_spliceSites ( Pattern.compile ( "," ).splitAsStream ( attributes.get ( "noncanonical_splicing" ).trim () )
                                                                  .map ( s -> s.split ( "\\+", 2 ) )
                                                                  .collect ( Collectors.toMap ( (a -> a[0]), a -> a.length > 1 ? a[1] : "" ) ) );
                }
            }
        }

        /*Set RibosomalSplippage attributes*/
        if (attributes.containsKey ( "ribosomal_slippage" )) {
            Ribosomal_Slippage ribosomal_slippage = new Ribosomal_Slippage ();
            if (attributes.get ( "ribosomal_slippage" ).equalsIgnoreCase ( "Y" )) {
                ribosomal_slippage.set_ribosomal_slippage ( true );
                ribosomal_slippage.setSlippage_frameshift ( Integer.parseInt ( attributes.get ( "slippage_fraqmeshift" ).replaceAll ( "^\"|\"$", "" ) ) );
                ribosomal_slippage.setSlippage_offset ( Integer.parseInt ( attributes.get ( "slippage_offset" ).replaceAll ( "^\"|\"$", "" ) ) );
                ribosomal_slippage.setSlippage_motif ( attributes.get ( "slippage_motif" ).replaceAll ( "^\"|\"$", "" ) );
            }
        }

        /*Set StopTranslationException attributes*/
        if (attributes.containsKey ( "stop_codon_readthru" )) {
            String attribute = attributes.get ( "stop_codon_readthru" ).replaceAll ( "^\"|\"$", "" );
            StopTranslationException stopTranslationException = new StopTranslationException ();
            String stop_codon_readthru;
            if (attribute.equalsIgnoreCase ( "Y" )) {
                stop_codon_readthru = "X";
                stopTranslationException.setStopTranslationException ( true );
            } else if (attribute.matches ( "^(Y|y)(:)[a-zA-z]" )) {

                stopTranslationException.setStopTranslationException ( true );
                String[] temp = attribute.split ( ":" );
                stop_codon_readthru = temp[1];
            }

        }

        /*Set StartTranslationException attributes*/
        if (attributes.containsKey ( "alternate_startcodon" )) {
            String attribute = attributes.get ( "alternate_startcodon" ).replaceAll ( "^\"|\"$", "" );
            if(attribute.matches ( ".*" )) {
                StartTranslationException startTranslationException = new StartTranslationException ();
                startTranslationException.setAlternateStartCodons ( Arrays.asList ( attribute.split ( "," ) ) );
                startTranslationException.setStartTranslationException ( true );
            }
        }

        /*Set RNA_Editing attributes*/
        if (attributes.containsKey ( "rna_editing" )) {
            String attribute = attributes.get ( "rna_editing" ).replaceAll ( "^\"|\"$", "" );
            if (attribute.matches ( ".*/.*/.*/.*/" )) {
                RNA_Editing rna_editing = new RNA_Editing ();
                String[] temp = attribute.split ( "/" );
                rna_editing.setSize ( Integer.parseInt ( temp[0] ) );
                rna_editing.setRegExp ( temp[1] );
                rna_editing.set_RNA_editing ( true );
                rna_editing.setReplacementString ( temp[2] );
                rna_editing.setNote ( temp[3] );
            }
        }

        /* Set StructuralSpecificatoons */



    }

    public String chooseAlignmentTool(AlignmentEvidence alignmentEvidence){

        return "jillion";
    }


    public void DetermineGeneStructure(ViralProtein viralProtein){

    boolean is_ribosomal_slippage = viralProtein.getGeneAttributes ().getRibosomal_slippage ().is_ribosomal_slippage ();
    boolean is_spliced = viralProtein.getGeneAttributes ().getSplicing ().isSpliced ();
    GeneStructure geneStructure = new GeneStructure ();
    ArrayList<Exon> exons = new ArrayList<Exon> ();
       if( !(is_ribosomal_slippage) && !(is_spliced)) {
           Exon exon = new Exon ();
           Range range = Range.of ( 0, 3 * (viralProtein.getSequence ().getLength ()) );
           exon.setRange ( range );
           exon.setOrder ( 1 );
           exons.add ( exon );
           geneStructure.setExons ( exons );
       }
       else{
         String spliceform = viralProtein.getGeneAttributes ().getSplicing ().getSpliceform ();
         System.out.println("printing spliceform"+spliceform);







       }
    }

}
