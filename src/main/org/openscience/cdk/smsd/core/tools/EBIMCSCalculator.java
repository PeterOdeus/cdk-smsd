/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.core.tools;

import org.openscience.cdk.smsd.SubStructureFactory;

/**
 *
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class EBIMCSCalculator extends SubStructureFactory {



    /**
     *
     * @param subStructureMode
     * @param bondTypeMatch
     * @param removeHydrogen
     * @param stereoFilter
     * @param fragmentFilter
     * @param energyFilter
     * @throws Exception
     */
    public EBIMCSCalculator(boolean subStructureMode, boolean bondTypeMatch, boolean removeHydrogen, boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) throws Exception {
        super(subStructureMode, bondTypeMatch, removeHydrogen, stereoFilter, fragmentFilter, energyFilter);
    }
    /**
     *
     * @param algorithmType 0: default, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
     * @param subStructureMode
     * @param bondTypeMatch
     * @param removeHydrogen
     * @param stereoFilter
     * @param fragmentFilter
     * @param energyFilter
     * @throws Exception
     */
    public EBIMCSCalculator(int algorithmType, boolean subStructureMode, boolean bondTypeMatch, boolean removeHydrogen, boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) throws Exception {
        super(algorithmType,subStructureMode, bondTypeMatch, removeHydrogen, stereoFilter, fragmentFilter, energyFilter);
    }
}
