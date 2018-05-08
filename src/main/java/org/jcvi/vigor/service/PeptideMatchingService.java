package org.jcvi.vigor.service;

import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.vigor.component.MaturePeptideMatch;
import org.jcvi.vigor.service.exception.ServiceException;

import java.io.File;
import java.util.List;

/**
 * TODO add filter argument
 */
public interface PeptideMatchingService {

    class Scores {
        final double identity;
        final double coverage;
        final double similarity;

        Scores(double identity, double coverage, double similarity) {
            this.identity = identity;
            this.coverage = coverage;
            this.similarity = similarity;
        }

        public static Scores of(double identity, double coverage, double similarity) {
            return new Scores(identity, coverage, similarity);
        }

    }

    List<MaturePeptideMatch> findPeptides(ProteinSequence protein, File peptideDatabase) throws ServiceException;
    List<MaturePeptideMatch> findPeptides(ProteinSequence protein, File peptideDatabase, Scores scores) throws ServiceException;
}
