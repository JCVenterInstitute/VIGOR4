package org.jcvi.vigor.service;

import java.util.Map;

import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.springframework.stereotype.Service;

@Service
public class EvaluateScores implements EvaluateModel {


    @Override
    public Model evaluate ( Model model, VigorConfiguration defaultConfiguration ) {

        Map<String, Double> scores = model.getScores();
        VigorConfiguration configuration = model.getAlignment().getViralProtein().getConfiguration();
        double alignmentScoreFactor = configuration.getOrDefault(ConfigurationParameters.ScoreFactorAlignment, 1d);

        double startScoreFactor = configuration.getOrDefault(ConfigurationParameters.ScoreFactorStart, 1d);
        double splicingScoreFactor = configuration.getOrDefault(ConfigurationParameters.ScoreFactorSplicing, 1d);
        double stopScoreFactor = configuration.getOrDefault(ConfigurationParameters.ScoreFactorStop, 1d);
        double leakyStopScoreFactor = configuration.getOrDefault(ConfigurationParameters.ScoreFactorLeakyStop, 1d);

        double alignmentScore = scores.getOrDefault(Scores.ALIGNMENT_SCORE,0d) * alignmentScoreFactor;
        double startCodonScore = scores.getOrDefault(Scores.START_CODON_SCORE,0d) * startScoreFactor;;
        double splicingScore = scores.getOrDefault(Scores.SPLICE_SCORE,0d) * splicingScoreFactor;
        double leakyStopScore = scores.getOrDefault(Scores.LEAKY_STOP_SCORE,0d) * leakyStopScoreFactor;
        double stopScore = scores.getOrDefault(Scores.STOP_CODON_SCORE,0d) * stopScoreFactor;;
        double totalScore;
        scores.put(Scores.ALIGNMENT_SCORE, alignmentScore);
        scores.put(Scores.START_CODON_SCORE, startCodonScore);
        scores.put(Scores.LEAKY_STOP_SCORE, leakyStopScore);
        scores.put(Scores.SPLICE_SCORE, splicingScore);
        scores.put(Scores.STOP_CODON_SCORE, stopScore);
        totalScore = alignmentScore + startCodonScore + leakyStopScore + splicingScore + stopScore;
        scores.put(Scores.TOTAL_SCORE, totalScore);
        return model;
    }
}
