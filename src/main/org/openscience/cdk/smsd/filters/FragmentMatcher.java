
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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.smsd.algorithm.cdk.CDKMCSHandler;
import org.openscience.cdk.smsd.algorithm.mcsplus.MCSPlusHandler;
import org.openscience.cdk.smsd.algorithm.single.SingleMappingHandler;
import org.openscience.cdk.smsd.algorithm.vflib.VFlibMCSHandler;
import org.openscience.cdk.smsd.factory.MCSFactory;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.IFragment;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

/**
 *
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 */
/**
 * @cdk.module smsd
 */
public class FragmentMatcher implements IFragment {

    // private boolean[][] localFlagMatrix;
    private int rowSize = 0;
    private int colSize = 0;
    private boolean[][] FlagMatrix = null;
    private MolHandler RMol;
    private MolHandler PMol;
    //private EBIMatrix SimMatrix = null;
    private IAtomContainerSet ReactantSet = DefaultChemObjectBuilder.getInstance().newAtomContainerSet();
    private IAtomContainerSet ProductSet = DefaultChemObjectBuilder.getInstance().newAtomContainerSet();
    private static List<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static List<TreeMap<Integer, Integer>> allMCS = null;
    double tanimoto = 0;
    double euclidean = 0;
    int stereoScore = 0;
    List<TreeMap<Integer, Integer>> GallMCS;
    TreeMap<Integer, Integer> GfirstSolution;
    List<Map<IAtom, IAtom>> GallAtomMCS;
    Map<IAtom, IAtom> GfirstAtomMCS;
    double Gtanimoto;
    double Geuclidean;
    private boolean removeHydrogen = false;

    //~--- constructors -------------------------------------------------------
    private void search() {

        //System.out.println("In FragmentMatcher Builder->FragmentMatcher->search");

        //SimMatrix = new JMatrix(rowSize, colSize);
        int SolutionSize = 0;
        try {
            for (int i = 0; i < ReactantSet.getAtomContainerCount(); i++) {

                IAtomContainer A = ReactantSet.getAtomContainer(i);
                RMol = new MolHandler(A, false);
                for (int j = 0; j < ProductSet.getAtomContainerCount(); j++) {

                    IAtomContainer B = ProductSet.getAtomContainer(j);
                    PMol = new MolHandler(B, false);

                    //printMolecules(RMol.getMolecule(),PMol.getMolecule());

                    //System.out.println("Mol " + Afp + " Mol" + Bfp);
                    //System.out.println("R Mol AtomCount: " + RMol.getMolecule().getAtomCount() + " P: Atom Count " + PMol.getMolecule().getAtomCount());




                    Builder();

                    //System.out.println("Best Solution size: " + SolutionSize);
                    //System.out.println("Present Solution size: " + firstSolution.size());

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

    private int checkCommonAtomCount(IAtomContainer reactantMolecule, IAtomContainer productMolecule) {
        ArrayList<String> a = new ArrayList<String>();
        for (int i = 0; i < reactantMolecule.getAtomCount(); i++) {
            if (removeHydrogen && !reactantMolecule.getAtom(i).getSymbol().equals("H")) {
                a.add(reactantMolecule.getAtom(i).getSymbol());
            } else {
                a.add(reactantMolecule.getAtom(i).getSymbol());
            }
        }


        int common = 0;
        for (int i = 0; i < productMolecule.getAtomCount(); i++) {

            if (a.contains(productMolecule.getAtom(i).getSymbol())) {
                a.remove(productMolecule.getAtom(i).getSymbol());
                common++;
            }
        }
        return common - a.size();
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
            if ((rBondCount <= 1 && rAtomCount == 1) || (pBondCount <= 1 && pAtomCount == 1)) {
//                System.out.println("Single Mapping");
                SingleMapping();
            } else if (rBondCount >= 6 && pBondCount >= 6) {
//                System.err.println("MCSPlus");
                MCSPlus();
                if (getFirstMapping() == null) {
                    System.gc();
//                    System.err.println("Many d-edges, switching to VF-McGregor");
                    VFLibMCS();
//                    System.out.println("Mapped with VF-McGregor");
                } else {
//                    System.out.println("Mapped with MCSPlus");
                }
            } else {
//                System.out.println("MCSPlus");
                MCSPlus();

            }
//            else if (commonAtoms == rAtomCount ||
//                    commonAtoms == pAtomCount) {
//                System.out.println("VF-McGregor");
//                VFLibMCS();
//
//            }

//            else if (rAtomCount != pAtomCount && rBondCount != pBondCount) {
//
//                System.out.println("CDKMCS-MCSPlus-VF-McGregor");
//                CDKMCS();
//                if (getFirstMapping() == null) {
//                    //Reset the mapping
//                    System.gc();
//                    System.err.println("Time out, switching to MCSPlus");
//                    MCSPlus();
//                    if (getFirstMapping() == null) {
//                        System.gc();
//                        System.err.println("Many d-edges, switching to VF-McGregor");
//                        VFLibMCS();
//                        System.out.println("Mapped with VF-McGregor");
//                    } else {
//                        System.err.println("Mapped with MCSPlus");
//                    }
//                } else {
//                    System.out.println("Mapped with CDKMCS");
//                }
//
//            } else {
////                System.err.println("MCSPlus");
//                MCSPlus();
//                if (getFirstMapping() == null) {
//                    System.gc();
//                    System.err.println("Many d-edges, switching to VF-McGregor");
//                    VFLibMCS();
//                    System.out.println("Mapped with VF-McGregor");
//                } else {
//                    System.out.println("Mapped with MCSPlus");
//                }
//
//
//            }

//            System.out.println("MCS solution count:" + FragmentMatcher.allMCS.size());
//            System.out.println("solution Size:" + FragmentMatcher.allMCS.firstElement().size());
//            for (Map.Entry<IAtom, IAtom> map : allAtomMCS.firstElement().entrySet()) {
//                System.out.println(map.getKey().getSymbol() + ":" + map.getValue().getSymbol());
//                System.out.println(map.getKey().getID() + ":" + map.getValue().getID());
//            }
//            long endTime = System.currentTimeMillis();
            System.gc();
//            System.out.println("Calculation Time: " + (endTime - startTime) * 0.001 + " seconds");
//            System.out.println("!Done!\n");

        } catch (Exception ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        }



    }

    private synchronized void CDKMCS() {
        try {
            CDKMCSHandler mcs = new CDKMCSHandler();
            mcs.set(RMol, PMol);
            mcs.searchMCS(removeHydrogen);
            firstMCS = mcs.getFirstMapping();
            allMCS = mcs.getAllMapping();
            allAtomMCS = mcs.getAllAtomMapping();
            atomsMCS = mcs.getFirstAtomMapping();

        } catch (CDKException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

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
        // System.out.println("In FragmentMatcher Builder->FragmentMatcher");
        GallMCS = new Vector<TreeMap<Integer, Integer>>();
        GfirstSolution =
                new TreeMap<Integer, Integer>();
        GallAtomMCS =
                new Vector<Map<IAtom, IAtom>>();
        GfirstAtomMCS =
                new HashMap<IAtom, IAtom>();
        Gtanimoto = 0;

        rowSize = A.getAtomContainerCount();
        colSize = B.getAtomContainerCount();

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

//~--- get methods --------------------------------------------------------
    @Override
    public boolean getFlag() {

        // int rowSize = EdMap.size();
        // int colSize = PdMap.size();
        boolean flag = false;

        for (int i = 0; i <
                rowSize; i++) {
            for (int j = 0; j <
                    colSize; j++) {
                if (FlagMatrix[i][j]) {
                    flag = true;

                    break;
                }

            }
        }

        return flag;
    }

    /**
     * 
     * @return
     */
    @Override
    public boolean[][] getFlagMatrix() {
        return FlagMatrix;
    }

    private void initFlagMatrix() {
        for (int i = 0; i <
                rowSize; i++) {
            for (int j = 0; j <
                    colSize; j++) {
                this.FlagMatrix[i][j] = false;
            }

        }
    }

//    public void printMolecules(IAtomContainer Molecule1, IAtomContainer Molecule2) {
//
//        System.out.println("Molecule 1");
//        for (int i = 0; i <
//                Molecule1.getAtomCount(); i++) {
//
//            System.out.print(Molecule1.getAtom(i).getSymbol() + " ");
//        }
//
//        System.out.println();
//        System.out.println("Molecule 2");
//        for (int i = 0; i <
//                Molecule2.getAtomCount(); i++) {
//
//            System.out.print(Molecule2.getAtom(i).getSymbol() + " ");
//        }
//
//        System.out.println();
//
//    }
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
        //System.out.println("Solution: "+ GfirstSolution );
        return GfirstSolution;
    }
}
//~ Formatted by Jindent --- http://www.jindent.com

