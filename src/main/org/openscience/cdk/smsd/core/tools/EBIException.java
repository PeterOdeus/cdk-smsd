

package org.openscience.cdk.smsd.core.tools;


public class EBIException extends Exception {

	private static final long serialVersionUID = 8371328769230823678L;

	/**
     * Constructs a new EBIException with the given message.
     *
     * @param message for the constructed exception
     */
    public EBIException(String message) {
        super( message );
    }

    /**
     * Constructs a new EBIException with the given message and the
     * Exception as cause.
     *
     * @param message for the constructed exception
     * @param cause   the Throwable that triggered this EBIException
     */
    public EBIException(String message, Throwable cause) {
        super(message, cause);
    }
}

