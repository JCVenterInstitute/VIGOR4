package org.jcvi.vigor.service.exception;

import org.jcvi.vigor.exception.VigorException;

public class ServiceException extends VigorException {

    public ServiceException ( String message ) {

        super(message);
    }

    public ServiceException ( String message, Throwable cause ) {

        super(message, cause);
    }

    public ServiceException ( Throwable cause ) {

        super(cause);
    }
}
