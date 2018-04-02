package org.jcvi.vigor.service;

import java.util.Map;

import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.VigorUtils;
import org.springframework.stereotype.Service;

@Service
public class EvaluateScores implements EvaluateModel {
	private float exonerateScoreFactor=1;
	private float startScoreFactor=1;
	private float splicingScoreFactor=1;
	private float stopScoreFactor=1;
	private float leakyStopScoreFactor=1;

	@Override
	public Model evaluate(Model model,VigorForm form) {
		
		Map<String,Double> scores = model.getScores();
		Map<String,String> params = form.getVigorParametersList();
		if (VigorUtils.is_Integer(params.get("exonerate_score_factor"))) {
			exonerateScoreFactor = Integer.parseInt(params.get("exonerate_score_factor"));
		}
		if (VigorUtils.is_Integer(params.get("start_score_factor"))) {
			startScoreFactor = Integer.parseInt(params.get("start_score_factor"));
		}
		if (VigorUtils.is_Integer(params.get("splicing_score_factor"))) {
			splicingScoreFactor = Integer.parseInt(params.get("splicing_score_factor"));
		}
		if (VigorUtils.is_Integer(params.get("stop_score_factor"))) {
			stopScoreFactor = Integer.parseInt(params.get("stop_score_factor"));
		}
		if (VigorUtils.is_Integer(params.get("leakystop_score_factor"))) {
			leakyStopScoreFactor = Integer.parseInt(params.get("leakystop_score_factor"));
		}
		double exonerateScore=0;
		double startCodonScore=0;
		double splicingScore=0;
		double leakyStopScore=0;
		double stopScore=0;
		double totalScore=0;
		if(scores.get("exonerateScore")!=null){
		exonerateScore = scores.get("exonerateScore")*exonerateScoreFactor;
		}
		if(scores.get("startCodonScore")!=null){
	    startCodonScore = scores.get("startCodonScore")*startScoreFactor;
		}
		if(scores.get("leakyStopScore") !=null){
	    leakyStopScore = scores.get("leakyStopScore")*leakyStopScoreFactor;
		}
		if(scores.get("spliceScore")!=null){
	    splicingScore = scores.get("spliceScore")*splicingScoreFactor;
		}
		if(scores.get("stopCodonScore")!=null){
	    stopScore = scores.get("stopCodonScore")*stopScoreFactor;
		}		
		scores.put("exonerateScore",exonerateScore);
		scores.put("startCodonScore",startCodonScore);
		scores.put("leakyStopScore",leakyStopScore);
		scores.put("spliceScore",splicingScore);
		scores.put("stopCodonScore",stopScore);
		totalScore = exonerateScore + startCodonScore + leakyStopScore + splicingScore + stopScore;
		scores.put("totalScore", totalScore);
		return model;
	}
	
}
