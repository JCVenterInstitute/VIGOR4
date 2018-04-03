package org.jcvi.vigor.service;

import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.vigor.component.MaturePeptide;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.service.exception.ServiceException;

import java.io.File;
import java.util.List;

/**
 * TODO return metadata as well as raw sequence.
 * TODO add filter argument
 */
public interface PeptideMatchingService {
    List<MaturePeptide> findPeptides(ViralProtein protein, File peptideDatabase) throws ServiceException;
}
