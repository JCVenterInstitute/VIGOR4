package com.vigor.blast;
import java.util.LinkedList;
import java.util.List;

import org.jcvi.jillion.experimental.align.blast.BlastHit;
import org.jcvi.jillion.experimental.align.blast.BlastVisitor;

public class BlastVistorImpl implements BlastVisitor {
	
	private String programName;
	private String programVersion;
	private String blastDb;
	private String queryId;
	private List<BlastHit> blastHitList;

	@Override
	public void visitEnd() {
		System.out.println("***************** in visit End *******************");
	}

	@Override
	public void visitInfo(String programName, String programVersion, String blastDb, String queryId) {
		System.out.println("***************** in visitInfo *******************");
		System.out.println("programName: "+programName);
		System.out.println("programVersion: "+programVersion);
		System.out.println("blastDb: "+blastDb);
		System.out.println("queryId: "+queryId);
		System.out.println("***************** in visitInfo *******************");
		this.programName = programName;
		this.programVersion = programVersion;
		this.blastDb = blastDb;
		this.queryId = queryId;
	}

	@Override
	public void visitHit(BlastHit hit) {
		System.out.println("***************** in visitHit *******************");
		System.out.println(hit.getBlastDbName());
		System.out.println(hit.getBlastProgramName());
		System.out.println(hit.getQueryId());
		System.out.println(hit.getSubjectDefinition());
		System.out.println(hit.getSubjectId());
		System.out.println(hit.getHsps());
		System.out.println(hit.getQueryLength());
		System.out.println(hit.getSubjectLength());
		System.out.println("***************** in visitHit *******************");
		if(blastHitList!=null){
			blastHitList.add(hit);
		}else{
			blastHitList = new LinkedList<BlastHit>();
			blastHitList.add(hit);
		}
	}

	public List<BlastHit> getBlastHitList() {
		return blastHitList;
	}

	public void setBlastHitList(List<BlastHit> blastHitList) {
		this.blastHitList = blastHitList;
	}

	public String getProgramName() {
		return programName;
	}

	public void setProgramName(String programName) {
		this.programName = programName;
	}

	public String getProgramVersion() {
		return programVersion;
	}

	public void setProgramVersion(String programVersion) {
		this.programVersion = programVersion;
	}

	public String getBlastDb() {
		return blastDb;
	}

	public void setBlastDb(String blastDb) {
		this.blastDb = blastDb;
	}

	public String getQueryId() {
		return queryId;
	}

	public void setQueryId(String queryId) {
		this.queryId = queryId;
	}
 
}
