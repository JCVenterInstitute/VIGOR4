package org.jcvi.vigor.exception;


public class VigorException extends Exception {
	private static final long serialVersionUID = -8729709786417672433L;

	public VigorException(String msg) {
		super(msg);
	}

	public VigorException(Throwable nestedException) {
		super(nestedException);
	}

	public VigorException(String message, Throwable nestedException) {
		super(message, nestedException);
	}
	
}
