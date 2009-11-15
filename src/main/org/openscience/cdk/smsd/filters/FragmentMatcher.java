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
package org.openscience.cdk.smsd.filters;

//~--- non-JDK imports --------------------------------------------------------
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
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
    double tanimoto = 0;
    double euclidean = 0;
    int stereoScore = 0;
    private Vector<TreeMap<Integer, Integer>> GallMCS;
    private TreeMap<Integer, Integer> GfirstSolution;
    private Vector<Map<IAtom, IAtom>> GallAtomMCS;
    private Map<IAtom, IAtom> GfirstAtomMCS;
    double Gtanimoto;
    double Geuclidean;
    private boolean removeHydrogen = false;

    //~--- constructors -------------------------------------------------------
    private void search() {

//        System.out.println("In FragmentMatcher Builder->FragmentMatcher->search");

        //SimMatrix = new JMatrix(rowSize, colSize);
        int SolutionSize = 0;
        try {
            for (int i = 0; i < ReactantSet.getAtomContainerCount(); i++) {

                IAtomContainer A = ReactantSet.getAtomContainer(i);
                RMol = new MolHandler(A, false);
                for (int j = 0; j < ProductSet.getAtomContainerCount(); j++) {

                    IAtomContainer B = ProductSet.getAtomContainer(j);
                    PMol = new MolHandler(B, false);

                    Builder();

                    if (SolutionSize < firstMCS.size()) {

                        //System.out.println("First Solution: " + firstSolution);

                        GfirstSolution.clear();
                        GfirstAtomMCS.clear();
                        GallAtomMCS.clear();
                        GallMCS.clear();

                        GfirstSolution.putAll(firstMCS);
                        GallMCS.addAll(allMCS);
                        GfirstAtomMCS.putAll(atomsMCS);
                        GallAtomMCS.addAll(allAtomMCS);
                        Gtanimoto = tanimoto;
                        Geuclidean = euclidean;

                        //setStereoScore();

                        SolutionSize = firstMCS.size();

                        // System.out.println("Best Solution size: " + SolutionSize);


                    } else if (SolutionSize == firstMCS.size()) {

                        /*Uncomments this if you want tp join all the solutions*/
                        //joinSolutions();
                            /*If you uncomment the above method then comment these*/
                        GallMCS.addAll(allMCS);
                        GallAtomMCS.addAll(allAtomMCS);
                    }
                }


            }
        } catch (Exception ex) {
            Logger.getLogger(FragmentMatcher.class.getName()).log(Level.SEVERE, null, ex);

        }

    }

    private synchronized void Builder() {

        try {

            int rBondCount = RMol.getMolecule().getBondCount();
            int pBondCount = PMol.getMolecule().getBondCount();

            int rAtomCount = RMol.getMolecule().getAtomCount();
            int pAtomCount = PMol.getMolecule().getAtomCount();

//            int commonAtoms = checkCommonAtomCount(RMol.getMolecule(), PMol.getMolecule());

            /*This is importsnt because CDK fails to Generate makeAtomsMapOfBondsMap if bonds are less than 2*/

//            long startTime = System.currentTimeMillis();
//
            if (rBondCount == 0 || rAtomCount == 1 || pBondCount == 0 || pAtomCount == 1) {
//                System.out.println("Single Mapping");
                SingleMapping();
            } else {
                if (rBondCount >= 6 && rBondCount >= 6) {
                    VFLibMCS();
                    if (getFirstMapping() == null) {
                        System.gc();
                    }
//                    System.out.println("Mapped with VFLibMCS");
                } else {
                    MCSPlus();
//                    System.out.println("Mapped with MCSPlus");
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

    private synchronized void MCSPlus() {
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

    private void VFLibMCS() {

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
     * @param A
     * @param B
     * @param removeHydrogen
     */
    public FragmentMatcher(IAtomContainerSet A, IAtomContainerSet B, boolean removeHydrogen) {

        this.removeHydrogen = removeHydrogen;
        GallMCS = new Vector<TreeMap<Integer, Integer>>();
        GfirstSolution = new TreeMap<Integer, Integer>();
        GallAtomMCS = new Vector<Map<IAtom, IAtom>>();
        GfirstAtomMCS = new HashMap<IAtom, IAtom>();
        Gtanimoto = 0;

        this.ReactantSet = A;
        this.ProductSet = B;
        search();

    }

    private void SingleMapping() {
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
    public Vector<Map<IAtom, IAtom>> getAllAtomMapping() {
        return GallAtomMCS;
    }

    /**
     *
     * @return
     */
    @Override
    public Vector<TreeMap<Integer, Integer>> getAllMapping() {
        return GallMCS;
    }

    /**
     *
     * @return
     */
    @Override
    public Map<IAtom, IAtom> getFirstAtomMapping() {
        /*for(IAtom I: GfirstAtomMCS){
        System.out.println(" A " + I.getSymbol());
        }*/
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

