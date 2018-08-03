package org.jcvi.vigor.exception;

public class VigorRuntimeException extends RuntimeException {

    public VigorRuntimeException() {
    }

    public VigorRuntimeException(String message) {
        super(message);
    }

    public VigorRuntimeException(String message, Throwable cause) {
        super(message, cause);
    }

    public VigorRuntimeException(Throwable cause) {
        super(cause);
    }

    public VigorRuntimeException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace) {
        super(message, cause, enableSuppression, writableStackTrace);
    }
}
