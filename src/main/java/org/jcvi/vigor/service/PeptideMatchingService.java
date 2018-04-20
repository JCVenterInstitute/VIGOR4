package org.jcvi.vigor.service;

import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.vigor.component.MaturePeptideMatch;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.service.exception.ServiceException;

import java.io.File;
import java.util.List;

/**
 * TODO add filter argument
 */
public interface PeptideMatchingService {
    List<MaturePeptideMatch> findPeptides(ProteinSequence protein, File peptideDatabase) throws ServiceException;
}
