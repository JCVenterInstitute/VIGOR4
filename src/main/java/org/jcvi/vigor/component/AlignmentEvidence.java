package org.jcvi.vigor.component;

import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import lombok.Data;

@Component
@Scope("prototype")
@Data
public class AlignmentEvidence {

    private String reference_db;
    private String matpep_db;

    public AlignmentEvidence ( String ref_db ) {

        this.reference_db = ref_db;
    }

    public AlignmentEvidence ( String ref_db, String matpep_db ) {

        this.reference_db = ref_db;
        this.matpep_db = matpep_db;
    }

    public AlignmentEvidence () {

    }

    public AlignmentEvidence copy () {

        return new AlignmentEvidence(this.reference_db, this.matpep_db);
    }
}


