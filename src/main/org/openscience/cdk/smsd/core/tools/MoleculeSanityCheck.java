/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.core.tools;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;

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
public class MoleculeSanityCheck {

    /**
     * 
     * @param molecule
     * @return
     */
    protected static IAtomContainer checkAndCleanMolecule(IAtomContainer molecule) {
        boolean isMarkush = false;
        for (IAtom atom : molecule.atoms()) {
            if (atom.getSymbol().equals("R")) {
                isMarkush = true;
                break;
            }
        }

        if (isMarkush) {
            System.err.println("Skipping Markush structure for sanity check");
        }

        // Check for salts and such
        if (!ConnectivityChecker.isConnected(molecule)) {
            // lets see if we have just two parts if so, we assume its a salt and just work
            // on the larger part. Ideally we should have a check to ensure that the smaller
            //  part is a metal/halogen etc.
            IMoleculeSet fragments = ConnectivityChecker.partitionIntoMolecules(molecule);
            if (fragments.getMoleculeCount() > 2) {
                System.err.println("More than 2 components. Skipped");
            } else {
                IMolecule frag1 = fragments.getMolecule(0);
                IMolecule frag2 = fragments.getMolecule(1);
                if (frag1.getAtomCount() > frag2.getAtomCount()) {
                    molecule = frag1;
                } else {
                    molecule = frag2;
                }
            }
        }
        fixAromaticity(molecule);
//
//        // do a aromaticity check
//        try {
//            CDKHueckelAromaticityDetector.detectAromaticity(molecule);
//        } catch (CDKException e) {
//            throw new CDKException("Error in aromaticity detection");
//        }

        return molecule;
    }

    public static void fixAromaticity(IAtomContainer m) {
        // need to find rings and aromaticity again since added H's

        IRingSet rs = null;
        try {
            AllRingsFinder arf = new AllRingsFinder();
            rs = arf.findAllRings(m);

            // SSSRFinder s = new SSSRFinder(m);
            // srs = s.findEssentialRings();

        } catch (Exception e) {
            e.printStackTrace();
        }

        try {
            // figure out which atoms are in aromatic rings:
            EBIAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(m);
            CDKHueckelAromaticityDetector.detectAromaticity(m);
            // figure out which rings are aromatic:
            RingSetManipulator.markAromaticRings(rs);
            // figure out which simple (non cycles) rings are aromatic:
            // HueckelAromaticityDetector.detectAromaticity(m, srs);
        } catch (Exception e) {
            e.printStackTrace();
        }


        // only atoms in 6 membered rings are aromatic
        // determine largest ring that each atom is a part of

        for (int i = 0; i <= m.getAtomCount() - 1; i++) {

            m.getAtom(i).setFlag(CDKConstants.ISAROMATIC, false);

            jloop:
            for (int j = 0; j <= rs.getAtomContainerCount() - 1; j++) {
                //logger.debug(i+"\t"+j);
                IRing r = (IRing) rs.getAtomContainer(j);
                if (!r.getFlag(CDKConstants.ISAROMATIC)) {
                    continue jloop;
                }

                boolean haveatom = r.contains(m.getAtom(i));

                //logger.debug("haveatom="+haveatom);

                if (haveatom && r.getAtomCount() == 6) {
                    m.getAtom(i).setFlag(CDKConstants.ISAROMATIC, true);
                }

            }

        }


    }
}
