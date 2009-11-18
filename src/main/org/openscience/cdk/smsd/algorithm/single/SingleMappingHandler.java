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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.smsd.algorithm.single;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.smsd.helper.FinalMappings;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.IMCS;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @cdk.module smsd
 */
public class SingleMappingHandler implements IMCS {

    private static List<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static List<TreeMap<Integer, Integer>> allMCS = null;
    private IAtomContainer ac1 = null;
    private IAtomContainer ac2 = null;

    public SingleMappingHandler() {


        allAtomMCS = new Vector<Map<IAtom, IAtom>>();
        atomsMCS = new HashMap<IAtom, IAtom>();
        firstMCS = new TreeMap<Integer, Integer>();
        allMCS = new Vector<TreeMap<Integer, Integer>>();

    }

    /** Creates a new instance of SearchCliques
     * @param Reactant
     * @param Product
     * @throws java.io.IOException
     *
     *
     */
    @Override
    public void set(MolHandler Reactant, MolHandler Product) throws IOException {

        this.ac1 = Reactant.getMolecule();
        this.ac2 = Product.getMolecule();


    }

    /** Creates a new instance of SearchCliques
     * @param ReactantMolFileName
     * @param ProductMolFileName
     * @throws java.io.IOException
     */
    @Override
    public void set(String ReactantMolFileName, String ProductMolFileName) throws IOException {


        String mol1 = ReactantMolFileName;
        String mol2 = ProductMolFileName;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);

        this.set(Reactant, Product);
        // System.exit(1);

    }

    /** Creates a new instance of SearchCliques
     * @param ReactantMol
     * @param ProductMol
     * @throws java.io.IOException
     */
    @Override
    public void set(IAtomContainer ReactantMol, IAtomContainer ProductMol) throws IOException {

        IAtomContainer mol1 = ReactantMol;
        IAtomContainer mol2 = ProductMol;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);

        this.set(Reactant, Product);
        // System.exit(1);

    }

    //Function is called by the main program and serves as a starting point for the comparision procedure.
    /**
     *
     * @return
     * @throws java.io.IOException
     */
    @Override
    public int searchMCS(boolean removeHydrogen) throws IOException, CDKException {


        SingleMapping SM = new SingleMapping();
        SM.getOverLaps(ac1, ac2, removeHydrogen);


        setAllMapping();
        setAllAtomMapping();
        setFirstMapping();
        setFirstAtomMapping();
        //setStereoScore();

        return 0;
    }

    public final void setAllMapping() {
//
//        int count_final_sol = 1;
        //System.out.println("Output of the final FinalMappings: ");
        try {

            List<TreeMap<Integer, Integer>> final_solution = FinalMappings.getInstance().getFinalMapping();
            int counter = 0;
            for (TreeMap<Integer, Integer> solution : final_solution) {
//                System.out.println("Numner of MCS solution: " + count_final_sol++);
                allMCS.add(counter++, solution);
            }

        } catch (Exception ex) {
            ex.printStackTrace();
        }


    }

    private final synchronized void setAllAtomMapping() {

        try {
            List<TreeMap<Integer, Integer>> final_solution = FinalMappings.getInstance().getFinalMapping();

            int counter = 0;
            for (TreeMap<Integer, Integer> solution : final_solution) {


                Map<IAtom, IAtom> atomMappings = new HashMap<IAtom, IAtom>();

//                System.out.println("Sol size " + solution.size());


                for (Map.Entry<Integer, Integer> map : solution.entrySet()) {

                    int IIndex = map.getKey();
                    int JIndex = map.getValue();
                    IAtom A = null;
                    IAtom B = null;


                    A = ac1.getAtom(IIndex);
                    B = ac2.getAtom(JIndex);

                    atomMappings.put(A, B);

                }

                allAtomMCS.add(counter++, atomMappings);
            }
        } catch (Exception I) {
            I.getCause();
        }



    }

    private final synchronized void setFirstMapping() {

        if (allMCS.size() > 0) {
            firstMCS = new TreeMap<Integer, Integer>(allMCS.get(0));
        }

    }

    private final synchronized void setFirstAtomMapping() {
        if (allAtomMCS.size() > 0) {
//            System.out.println("In MCS+handle First Atom MCS: " + allAtomMCS.get(0));
            atomsMCS = new HashMap<IAtom, IAtom>(allAtomMCS.get(0));
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
