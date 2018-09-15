package org.jcvi.vigor.component;

import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import lombok.Data;

import java.io.File;

@Component
@Scope("prototype")
@Data
public class AlignmentEvidence {

    private String reference_db;
    private String matpep_db;
    private File results_directory;
    private File raw_alignment;

    public AlignmentEvidence ( String ref_db ) {

        this.reference_db = ref_db;
    }

    public AlignmentEvidence ( String ref_db, String matpep_db, File results_directory, File raw_alignment ) {

        this.reference_db = ref_db;
        this.matpep_db = matpep_db;
        this.results_directory = results_directory;
        this.raw_alignment = raw_alignment;
    }

    public AlignmentEvidence () {

    }

    public AlignmentEvidence copy () {

        return new AlignmentEvidence(this.reference_db, this.matpep_db, results_directory, raw_alignment);
    }
}


