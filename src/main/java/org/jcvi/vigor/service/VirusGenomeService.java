package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.List;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.springframework.stereotype.Service;

import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.VigorUtils;

@Service
public class VirusGenomeService {

	/**
	 * 
	 * @param sequence
	 * @param form
	 * @return List<Range>: Only ranges having length greater or equal to
	 *         min_gap_length will be considered as a sequence gap
	 */

	public List<Range> findSequenceGapRanges(String minGapLenString, NucleotideSequence sequence) {
		List<Range> RangesOfNs = sequence.getRangesOfNs();
		long minGapLength = 20;
		if (VigorUtils.is_Integer(minGapLenString)) {
			minGapLength = Long.parseLong(minGapLenString);
		}
		List<Range> filteredRangesOfNs = new ArrayList<Range>();
		if (!RangesOfNs.isEmpty()) {
			Range previousRange = Range.of(0, 0);
			for (int i = 0; i < RangesOfNs.size(); i++) {
				if (RangesOfNs.get(i).getLength() >= minGapLength) {
					Range currentRange = RangesOfNs.get(i);
					if (previousRange.getBegin() != 0 && previousRange.getEnd() != 0) {
						if (previousRange.getEnd() <= currentRange.getBegin() + 6) {
							previousRange = Range.of(previousRange.getBegin(), currentRange.getEnd());
						}
					}

				}
				filteredRangesOfNs.add(previousRange);
			}
		}

		return RangesOfNs;

	}

}
