package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.utils.VigorUtils;

	@Service
	public class VirusGenomeService {

		/**
		 * 
		 * @param sequence
		 * @return List<Range>: Only ranges having length greater or equal to
		 *         min_gap_length will be considered as a sequence gap
		 */

		public static List<Range> findSequenceGapRanges(String minGapLenString, NucleotideSequence sequence) {
			List<Range> rangesOfNs = sequence.getRangesOfNs();
			long minGapLength = 20;
			if (VigorUtils.is_Integer(minGapLenString)) {
				minGapLength = Long.parseLong(minGapLenString);
			}
			List<Range> filteredRangesOfNs = new ArrayList<Range>();
			if (!rangesOfNs.isEmpty()) {
				Range previousRange = Range.of(0, 0);
				for (int i = 0; i < rangesOfNs.size(); i++) {
					if (rangesOfNs.get(i).getLength() >= minGapLength) {
						Range currentRange = rangesOfNs.get(i);
						if (previousRange.getBegin() != 0 && previousRange.getEnd() != 0) {
							if (previousRange.getEnd() <= currentRange.getBegin() + 6) {
								previousRange = Range.of(previousRange.getBegin(), currentRange.getEnd());
							}
						}

					}
					filteredRangesOfNs.add(previousRange);
				}
			}
		  return rangesOfNs;

		}

		public static Map<Frame,List<Long>> findInternalStops(NucleotideSequence NTSequence){
			Map<Frame,List<Long>> stops = IupacTranslationTables.STANDARD.findStops(NTSequence);
			stops = VigorFunctionalUtils.frameToSequenceFrame(stops);
			return stops;
		}

	}


