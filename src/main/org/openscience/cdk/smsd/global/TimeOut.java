/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.global;

/**
 *
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class TimeOut {

    private static TimeOut INSTANCE = null;
    private double time = -1;

    /**
     *
     * @return
     */
    public static synchronized TimeOut getInstance() {
        if (INSTANCE == null) {

            // it's ok, we can call this constructor
            INSTANCE = new TimeOut();
        }

        return INSTANCE;
    }

    protected TimeOut() {
    }

    public void setTimeOut(double timeout) {

        this.time = timeout;
    }

    public double getTimeOut() {

        return time;
    }
}
