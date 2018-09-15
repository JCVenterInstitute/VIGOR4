package org.jcvi.vigor.service;

import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.VigorConfiguration;

public interface EvaluateModel {

    Model evaluate (Model model, VigorConfiguration configuration );
}
