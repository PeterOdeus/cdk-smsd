/* Copyright (C) 2006-2009  Syed Asad Rahman {asad@ebi.ac.uk}
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
 * MERCHANTABILITY or FITNESS FOR sourceAtom PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.smsd.algorithm.cdk;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smsd.helper.FinalMappings;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.smsd.interfaces.IMCS;

/**
 * @cdk.module smsd
 */
public class CDKMCSHandler implements IMCS {

//    //~--- fields -------------------------------------------------------------
    private IAtomContainer source;
    private IAtomContainer target;
    private boolean RonPFlag = false;
    private static List<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static List<TreeMap<Integer, Integer>> allMCS = null;

    //~--- constructors -------------------------------------------------------
    /** Creates a new instance of MappingHandler
     */
    public CDKMCSHandler() {

        allAtomMCS = new ArrayList<Map<IAtom, IAtom>>();
        atomsMCS = new HashMap<IAtom, IAtom>();
        firstMCS = new TreeMap<Integer, Integer>();
        allMCS = new ArrayList<TreeMap<Integer, Integer>>();
    }


    /**
     * @param source
     * @param target
     */
    @Override
    public void init(IAtomContainer source, IAtomContainer target) {

        IAtomContainer mol1 = source;
        IAtomContainer mol2 = target;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);

        init(Reactant, Product);

    }

    /**
     * @param source
     * @param target
     */
    public void init(IMolecule source, IMolecule target) throws CDKException {

        IMolecule mol1 = source;
        IMolecule mol2 = target;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);

        init(Reactant, Product);
    }

    /**
     * @param sourceMolFileName
     * @param targetMolFileName
     */
    @Override
    public void init(String sourceMolFileName, String targetMolFileName) {

        String mol1 = sourceMolFileName;
        String mol2 = targetMolFileName;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);
        init(Reactant, Product);


    }

    /**
     * Set the CDKMCS software
     * 
     * @param reactant
     * @param product
     */
    @Override
    public void init(MolHandler reactant, MolHandler product) {


        this.source = reactant.getMolecule();
        this.target = product.getMolecule();

    }
    /**
     * 
     * @return
     */
    @Override
    public int searchMCS(boolean removeHydrogen) throws CDKException {


        CDKRMapHandler rmap = new CDKRMapHandler();

        try {

            if ((source.getAtomCount() == target.getAtomCount()) && source.getBondCount() == target.getBondCount()) {
                RonPFlag = true;
                //System.err.print("True");

                rmap.calculateOverlapsAndReduceExactMatch(source, target);

            } else if (source.getAtomCount() > target.getAtomCount() && source.getBondCount() != target.getBondCount()) {
                RonPFlag = true;
                //System.err.print("True");

                rmap.calculateOverlapsAndReduce(source, target);

            } else {
                RonPFlag = false;
                rmap.calculateOverlapsAndReduce(target, source);


            }

            setAllMapping();
            setAllAtomMapping();
            setFirstMapping();
            setFirstAtomMapping();

        } catch (CDKException e) {
            rmap = null;
//            System.err.println("WARNING: CDKMCS: most probably time out error ");
        }

        return 0;
    }

    /**
     * 
     * @param mol
     * @param mcss
     * @return
     * @throws CDKException 
     */
    protected IMoleculeSet getUncommon(IAtomContainer mol, IAtomContainer mcss) throws CDKException {
        ArrayList<Integer> atomSerialsToDelete = new ArrayList<Integer>();

        List matches = CDKMCS.getSubgraphAtomsMaps(mol, mcss);
        List mapList = (List) matches.get(0);
        for (Object o : mapList) {
            CDKRMap rmap = (CDKRMap) o;
            atomSerialsToDelete.add(rmap.getId1());
        }

        // at this point we have the serial numbers of the bonds to delete
        // we should get the actual bonds rather than delete by serial numbers
        ArrayList<IAtom> atomsToDelete = new ArrayList<IAtom>();
        for (Integer serial : atomSerialsToDelete) {
            atomsToDelete.add(mol.getAtom(serial));
        }

        // now lets get rid of the bonds themselves
        for (IAtom atom : atomsToDelete) {
            mol.removeAtomAndConnectedElectronContainers(atom);
        }

        // now we probably have a set of disconnected components
        // so lets get a set of individual atom containers for
        // corresponding to each component
        return ConnectivityChecker.partitionIntoMolecules(mol);
    }

    //~--- get methods --------------------------------------------------------
    private final synchronized void setAllMapping() {

        //int count_final_sol = 1;
        //System.out.println("Output of the final FinalMappings: ");
        try {
            List<TreeMap<Integer, Integer>> sol = FinalMappings.getInstance().getFinalMapping();
            int counter = 0;
            for (TreeMap<Integer, Integer> final_solution : sol) {
                TreeMap<Integer, Integer> atomMappings = new TreeMap<Integer, Integer>();
                for (Map.Entry<Integer, Integer> Solutions : final_solution.entrySet()) {

                    int IIndex = Solutions.getKey().intValue();
                    int JIndex = Solutions.getValue().intValue();

                    if (RonPFlag) {
                        atomMappings.put(IIndex, JIndex);
                    } else {
                        atomMappings.put(JIndex, IIndex);
                    }
                }


                if (!allMCS.contains(atomMappings)) {
//                    System.out.println("atomMappings: " + atomMappings);
                    allMCS.add(counter++, atomMappings);
                }
            }

        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    private final synchronized void setAllAtomMapping() {
        List<TreeMap<Integer, Integer>> sol = allMCS;

        int counter = 0;
        for (TreeMap<Integer, Integer> final_solution : sol) {

            Map<IAtom, IAtom> atomMappings = new HashMap<IAtom, IAtom>();


            for (Map.Entry<Integer, Integer> Solutions : final_solution.entrySet()) {

                int IIndex = Solutions.getKey().intValue();
                int JIndex = Solutions.getValue().intValue();

                IAtom sourceAtom = null;
                IAtom targetAtom = null;

                sourceAtom = source.getAtom(IIndex);
                targetAtom = target.getAtom(JIndex);
                atomMappings.put(sourceAtom, targetAtom);


            }

            allAtomMCS.add(counter++, atomMappings);
        }



    }

    private final synchronized void setFirstMapping() {

        if (allMCS.size() > 0) {
            firstMCS = new TreeMap<Integer, Integer>(allMCS.get(0));
        }

    }

    private final synchronized void setFirstAtomMapping() {
        if (atomsMCS.size() > 0) {
            atomsMCS = new TreeMap<IAtom, IAtom>(allAtomMCS.get(0));
        }

    }

    @Override
    public List<TreeMap<Integer, Integer>> getAllMapping() {
        return allMCS;
    }

    @Override
    public TreeMap<Integer, Integer> getFirstMapping() {
        return firstMCS;
    }

    @Override
    public List<Map<IAtom, IAtom>> getAllAtomMapping() {
        return allAtomMCS;
    }

    @Override
    public Map<IAtom, IAtom> getFirstAtomMapping() {
        return atomsMCS;
    }
}
