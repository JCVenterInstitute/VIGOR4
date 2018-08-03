package org.jcvi.vigor.component;

import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import lombok.Data;

@Component
@Scope("prototype")
@Data
public class GeneAttributes {

    private Ribosomal_Slippage ribosomal_slippage = Ribosomal_Slippage.NO_SLIPPAGE;
    private Splicing splicing = Splicing.NO_SPLICING;
    private StartTranslationException startTranslationException = StartTranslationException.NO_EXCEPTION;
    private StopTranslationException stopTranslationException = StopTranslationException.NO_EXCEPTION;
    private RNA_Editing rna_editing = RNA_Editing.NO_EDITING;
    private StructuralSpecifications structuralSpecifications = new StructuralSpecifications();
}
