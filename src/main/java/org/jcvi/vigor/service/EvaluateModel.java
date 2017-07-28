package org.jcvi.vigor.service;

import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;

public interface EvaluateModel {

	public Model determine(Model model,VigorForm form);
}
