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
 * MERCHANTABILITY or FITNESS FOR source PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.smsd.filters;

//~--- non-JDK imports --------------------------------------------------------
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.smsd.algorithm.mcsplus.MCSPlusHandler;
import org.openscience.cdk.smsd.algorithm.single.SingleMappingHandler;
import org.openscience.cdk.smsd.algorithm.vflib.VFlibMCSHandler;
import org.openscience.cdk.smsd.factory.MCSFactory;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.smsd.interfaces.IMCSBase;

/**
 * @cdk.module smsd
 */
public class FragmentMatcher implements IMCSBase {

    private MolHandler RMol;
    private MolHandler PMol;
    private IAtomContainerSet ReactantSet = DefaultChemObjectBuilder.getInstance().newAtomContainerSet();
    private IAtomContainerSet ProductSet = DefaultChemObjectBuilder.getInstance().newAtomContainerSet();
    private static List<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static List<TreeMap<Integer, Integer>> allMCS = null;
    private List<TreeMap<Integer, Integer>> GallMCS;
    private TreeMap<Integer, Integer> GfirstSolution;
    private List<Map<IAtom, IAtom>> GallAtomMCS;
    private Map<IAtom, IAtom> GfirstAtomMCS;
    private boolean removeHydrogen = false;

    //~--- constructors -------------------------------------------------------
    private void search() {
        int SolutionSize = 0;
        try {
            for (int i = 0; i < ReactantSet.getAtomContainerCount(); i++) {

                IAtomContainer source = ReactantSet.getAtomContainer(i);
                RMol = new MolHandler(source, false);
                for (int j = 0; j < ProductSet.getAtomContainerCount(); j++) {

                    IAtomContainer target = ProductSet.getAtomContainer(j);
                    PMol = new MolHandler(target, false);

                    builder();

                    if (SolutionSize < firstMCS.size()) {

                        GfirstSolution.clear();
                        GfirstAtomMCS.clear();
                        GallAtomMCS.clear();
                        GallMCS.clear();

                        GfirstSolution.putAll(firstMCS);
                        GallMCS.addAll(allMCS);
                        GfirstAtomMCS.putAll(atomsMCS);
                        GallAtomMCS.addAll(allAtomMCS);
                        SolutionSize = firstMCS.size();
                    } else if (SolutionSize == firstMCS.size()) {
                        GallMCS.addAll(allMCS);
                        GallAtomMCS.addAll(allAtomMCS);
                    }
                }


            }
        } catch (Exception ex) {
            Logger.getLogger(FragmentMatcher.class.getName()).log(Level.SEVERE, null, ex);

        }

    }

    private synchronized void builder() {

        try {

            int rBondCount = RMol.getMolecule().getBondCount();
            int pBondCount = PMol.getMolecule().getBondCount();

            int rAtomCount = RMol.getMolecule().getAtomCount();
            int pAtomCount = PMol.getMolecule().getAtomCount();

            if (rBondCount == 0 || rAtomCount == 1 || pBondCount == 0 || pAtomCount == 1) {
//                System.out.println("Single Mapping");
                singleMapping();
            } else {
                if (rBondCount >= 6 && rBondCount >= 6) {
                    vfLibMCS();
                    if (getFirstMapping() == null) {
                        System.gc();
                    }
//                    System.out.println("Mapped with vfLibMCS");
                } else {
                    mcsPlus();
//                    System.out.println("Mapped with mcsPlus");
                }
            }
            System.gc();
        } catch (Exception ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

//    private synchronized void CDKMCS() {
//        try {
//            CDKMCSHandler mcs = new CDKMCSHandler();
//            mcs.set(RMol, PMol);
//            mcs.searchMCS(removeHydrogen);
//            firstMCS = mcs.getFirstMapping();
//            allMCS = mcs.getAllMapping();
//            allAtomMCS = mcs.getAllAtomMapping();
//            atomsMCS = mcs.getFirstAtomMapping();
//
//        } catch (CDKException ex) {
//            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
//        }
//
//    }
    private synchronized void mcsPlus() {
        try {
            MCSPlusHandler mcs = new MCSPlusHandler();
            mcs.set(RMol, PMol);
            mcs.searchMCS(removeHydrogen);
            firstMCS = mcs.getFirstMapping();
            allMCS = mcs.getAllMapping();
            allAtomMCS = mcs.getAllAtomMapping();
            atomsMCS = mcs.getFirstAtomMapping();

        } catch (IOException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private void vfLibMCS() {

        try {
            VFlibMCSHandler mcs = new VFlibMCSHandler();
            mcs.set(RMol, PMol);
            mcs.searchMCS(removeHydrogen);


            firstMCS = mcs.getFirstMapping();
            allMCS = mcs.getAllMapping();
            allAtomMCS = mcs.getAllAtomMapping();
            atomsMCS = mcs.getFirstAtomMapping();

        } catch (IOException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    /**
     *
     * @param source
     * @param target
     * @param removeHydrogen
     */
    public FragmentMatcher(IAtomContainerSet source, IAtomContainerSet target, boolean removeHydrogen) {

        this.removeHydrogen = removeHydrogen;
        GallMCS = new ArrayList<TreeMap<Integer, Integer>>();
        GfirstSolution = new TreeMap<Integer, Integer>();
        GallAtomMCS = new ArrayList<Map<IAtom, IAtom>>();
        GfirstAtomMCS = new HashMap<IAtom, IAtom>();
        this.ReactantSet = source;
        this.ProductSet = target;
        search();

    }

    private void singleMapping() {
        try {

            SingleMappingHandler mcs = new SingleMappingHandler();
            mcs.set(RMol, PMol);
            mcs.searchMCS(removeHydrogen);

            firstMCS = mcs.getFirstMapping();
            allMCS = mcs.getAllMapping();
            allAtomMCS = mcs.getAllAtomMapping();
            atomsMCS = mcs.getFirstAtomMapping();


        } catch (IOException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     *
     * @return
     */
    @Override
    public List<Map<IAtom, IAtom>> getAllAtomMapping() {
        return GallAtomMCS;
    }

    /**
     *
     * @return
     */
    @Override
    public List<TreeMap<Integer, Integer>> getAllMapping() {
        return GallMCS;
    }

    /**
     *
     * @return
     */
    @Override
    public Map<IAtom, IAtom> getFirstAtomMapping() {
        return GfirstAtomMCS;
    }

    /**
     *
     * @return
     */
    @Override
    public TreeMap<Integer, Integer> getFirstMapping() {
        return GfirstSolution;
    }
}

