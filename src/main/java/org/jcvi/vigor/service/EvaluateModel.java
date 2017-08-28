package org.jcvi.vigor.service;

import java.util.List;

import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;

public interface EvaluateModel {

	public List<Model> determine(Model model,VigorForm form);
}
