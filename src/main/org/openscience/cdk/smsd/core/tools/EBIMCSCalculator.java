/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.core.tools;

import org.openscience.cdk.smsd.SubStructureFactory;

/**
 *
 * Copyright (C) 2006-2009  Syed Asad Rahman {asad@ebi.ac.uk}
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/**
 * @cdk.module smsd
 */
public class EBIMCSCalculator extends SubStructureFactory {

    /**
    @param subStructureMode true for fast substructure search without
     * exhaustive MCS else false
     * @param bondTypeMatch   true will considered bond types for mapping else false
     * @param removeHydrogen true if Hydrogen are not to be mapped else false
     * @param stereoFilter   true if stereo match is considered else false
     * @param fragmentFilter true if fragement filter is switched on else false
     * @param energyFilter   true if bond energy filter is switched on else false
     * @throws Exception
     */
    public EBIMCSCalculator(boolean subStructureMode, boolean bondTypeMatch, boolean removeHydrogen, boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) throws Exception {
        super(subStructureMode, bondTypeMatch, removeHydrogen, stereoFilter, fragmentFilter, energyFilter);
    }

    /**
     *
     * @param algorithmType 0: default, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
     * @param subStructureMode true for fast substructure search without
     * exhaustive MCS else false
     * @param bondTypeMatch   true will considered bond types for mapping else false
     * @param removeHydrogen true if Hydrogen are not to be mapped else false
     * @param stereoFilter   true if stereo match is considered else false
     * @param fragmentFilter true if fragement filter is switched on else false
     * @param energyFilter   true if bond energy filter is switched on else false
     * @throws Exception
     */
    public EBIMCSCalculator(int algorithmType, boolean subStructureMode, boolean bondTypeMatch, boolean removeHydrogen, boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) throws Exception {
        super(algorithmType, subStructureMode, bondTypeMatch, removeHydrogen, stereoFilter, fragmentFilter, energyFilter);
    }
}
