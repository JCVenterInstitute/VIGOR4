package org.jcvi.vigor.service.exception;

import org.jcvi.vigor.exception.VigorException;

public class UserFacingException extends VigorException {

    public UserFacingException(String msg) {
        super(msg);
    }

    public UserFacingException(Throwable nestedException) {
        super(nestedException);
    }

    public UserFacingException(String message, Throwable nestedException) {
        super(message, nestedException);
    }
}
