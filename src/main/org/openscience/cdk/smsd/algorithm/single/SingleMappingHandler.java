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
package org.openscience.cdk.smsd.algorithm.single;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.smsd.helper.FinalMappings;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.IMCSAlgorithm;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;

/**
 * @cdk.module smsd
 */
public class SingleMappingHandler implements IMCSAlgorithm {

    private static List<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static List<TreeMap<Integer, Integer>> allMCS = null;
    private IAtomContainer source = null;
    private IAtomContainer target = null;
    private boolean removeHydrogen = false;

    /**
     * 
     * @param removeH true
     */
    public SingleMappingHandler(boolean removeH) {

        this.removeHydrogen = removeH;
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
    public void set(IAtomContainer source, IAtomContainer target) {

        IAtomContainer mol1 = source;
        IAtomContainer mol2 = target;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);

        set(Reactant, Product);

    }

    /**
     * @param source
     * @param target
     */
    public void set(IMolecule source, IMolecule target) throws CDKException {

        IMolecule mol1 = source;
        IMolecule mol2 = target;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);

        set(Reactant, Product);
    }

    /**
     * @param sourceMolFileName
     * @param targetMolFileName
     */
    @Override
    public void set(String sourceMolFileName, String targetMolFileName) {

        String mol1 = sourceMolFileName;
        String mol2 = targetMolFileName;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);
        set(Reactant, Product);


    }

    /**
     * @param source
     * @param target
     */
    @Override
    public void set(MolHandler source, MolHandler target) {
        this.source = source.getMolecule();
        this.target = target.getMolecule();
    }
    //Function is called by the main program and serves as a starting point for the comparision procedure.

    /**
     *
     * @return
     * @throws java.io.IOException
     */
    @Override
    public int searchMCS() throws IOException, CDKException {


        SingleMapping singleMapping = new SingleMapping();
        singleMapping.getOverLaps(source, target, removeHydrogen);


        setAllMapping();
        setAllAtomMapping();
        setFirstMapping();
        setFirstAtomMapping();
        //setStereoScore();

        return 0;
    }

    /**
     * Set the mappings
     */
    public final void setAllMapping() {
        try {

            List<TreeMap<Integer, Integer>> final_solution = FinalMappings.getInstance().getFinalMapping();
            int counter = 0;
            for (TreeMap<Integer, Integer> solution : final_solution) {
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
                for (Map.Entry<Integer, Integer> map : solution.entrySet()) {

                    int IIndex = map.getKey();
                    int JIndex = map.getValue();
                    IAtom sourceAtom = null;
                    IAtom targetAtom = null;


                    sourceAtom = source.getAtom(IIndex);
                    targetAtom = target.getAtom(JIndex);

                    atomMappings.put(sourceAtom, targetAtom);

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
