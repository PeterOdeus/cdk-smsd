/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.core.tools;

import java.text.SimpleDateFormat;
import java.util.TimeZone;

/**
 * This is used in MCSPlus, RGraph and UniversialIsomorphismTester
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class EBITimeManager {

    /**
     * 
     * @param startTime
     * @return
     */
    private static EBITimeManager ref = null;
    private double startTime;
    private SimpleDateFormat dateFormat;

    public EBITimeManager() {

        dateFormat = new SimpleDateFormat("HH:mm:ss");

        dateFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
        startTime = System.currentTimeMillis();
    }

    /**
     *
     * @return
     */
    public double getElapsedTimeInHours() {
        double currentTime = System.currentTimeMillis();

        return (currentTime - startTime) / (60 * 60 * 1000);


    }

    /**
     * 
     * @return 
     */
    public double getElapsedTimeInMinutes() {

        //long diffSeconds = diff / 1000;
        //long diffMinutes = diff / (60 * 1000);
        //long diffHours = diff / (60 * 60 * 1000);
        //long diffDays = diff / (24 * 60 * 60 * 1000);

        double currentTime = System.currentTimeMillis();

        return (currentTime - startTime) / (60 * 1000);

    }

    /**
     *
     * @return
     */
    public double getElapsedTimeInSeconds() {



        double currentTime = System.currentTimeMillis();

        return ((currentTime - startTime) / 1000);

    }

    public double getElapsedTimeInMillisSeconds() {



        double currentTime = System.currentTimeMillis();

        return (currentTime - startTime);

    }
}
