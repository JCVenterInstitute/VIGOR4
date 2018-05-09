package org.jcvi.vigor.service;

import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.service.exception.ServiceException;

import java.util.List;

public interface AlignmentService {

    List<Alignment> getAlignment(AlignmentEvidence alignmentEvidence, VirusGenome virusGenome, String referenceDB, String workspace) throws ServiceException;
}