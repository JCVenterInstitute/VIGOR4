package org.jcvi.vigor.service;

import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentTool;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.VigorConfiguration;

import java.util.List;

public interface AlignmentService {

    List<Alignment> getAlignment (VigorConfiguration config, VirusGenome virusGenome, String referenceDB, String workspace ) throws ServiceException;
    AlignmentTool getAlignmentTool();
}
