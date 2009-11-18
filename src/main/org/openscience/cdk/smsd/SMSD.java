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
package org.openscience.cdk.smsd;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.smsd.factory.MCSFactory;
import org.openscience.cdk.smsd.factory.SubGraphFactory;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.IMCSAlgorithm;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;

/**
 * @cdk.module smsd
 * <p>SMSD algorithm is described in {@cdk.cite SMSD2009}. Please cite {@cdk.cite SMSD2009}
 * if you are using SMSD module
 * </p>
 *
 * <p> <H2>Usage Example</H2>
 *          EBIMDLReader R1 = new EBIMDLReader(new FileInputStream(mol1));
            IMolecule A1 = R1.getMolecule();

            EBIMDLReader R2 = new EBIMDLReader(new FileInputStream(mol2));
            IMolecule A2 = R2.getMolecule();

            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(A1);
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(A2);

            A1 = (IMolecule) AtomContainerManipulator.removeHydrogens(A1);
            A2 = (IMolecule) AtomContainerManipulator.removeHydrogens(A2);

            CDKHueckelAromaticityDetector.detectAromaticity(A1);
            CDKHueckelAromaticityDetector.detectAromaticity(A2);

            boolean substructure = false;
            boolean bondSensitive = true;
            boolean removeHydrogen = true;
            boolean stereoMatch = true;
            boolean fragmentMinimization = true;
            boolean energyMinimization = true;

            SMSD comparison = new SMSD(substructure, bondSensitive,
                    removeHydrogen, stereoMatch, fragmentMinimization, energyMinimization);

            comparison.init(A1, A2);


            //Get modified Query and Target Molecules as Mappings will correspond to these molecules
            IAtomContainer Query = comparison.getReactantMolecule();
            IAtomContainer Target = comparison.getProductMolecule();


            for (Map.Entry<Integer, Integer> mappings : comparison.getFirstMapping().entrySet()) {
                //Get the mapped atom number in Query Molecule
                int queryMappingNumber = mappings.getKey();
                //Get the mapped atom number in Target Molecule
                int targetMappingNumber = mappings.getValue();

                //Get the mapped atom in Query Molecule
                IAtom queryAtom = Query.getAtom(queryMappingNumber);
                //Get the mapped atom in Target Molecule
                IAtom targetAtom = Target.getAtom(targetMappingNumber);
                //Print mapped atom numbers
                System.out.println(queryMappingNumber + " " +
                        (targetMappingNumber));
                //Print mapped atoms
                System.out.println(queryAtom.getSymbol() + " " +
                        targetAtom.getSymbol());
            }
            System.out.println("");

            System.out.println("");

            System.out.println("Stereo Match: " + comparison.getStereoScore(0));
            System.out.println("Stereo different: " + comparison.isStereoMisMatch());
            System.out.println("Fragment Size: " + comparison.getFragmentSize(0));
            System.out.println("Tanimoto Similarity Score: " + comparison.getTanimotoSimilarity());
            System.out.println("Tanimoto Euclidean Distance: " + comparison.getEuclideanDistance());
            System.out.println("");

 * </p>
 *
 */
public class SMSD implements IMCSAlgorithm {

    IMCSAlgorithm comparison;

    /**
     * 
     * @param subStructureMode true for fast substructure search without
     * exhaustive MCS else false
     * @param bondTypeFlag   true will considered bond types for mapping else false
     * @param removeHydrogen true if Hydrogen are not to be mapped else false
     * @param stereoFilter   true if stereo match is considered else false
     * @param fragmentFilter true if fragement filter is switched on else false
     * @param energyFilter   true if bond energy filter is switched on else false
     * @throws Exception
     */
    public SMSD(boolean subStructureMode, boolean bondTypeFlag, boolean removeHydrogen, boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) throws Exception {
        if (subStructureMode) {
            comparison = new SubGraphFactory(bondTypeFlag, removeHydrogen, stereoFilter, fragmentFilter, energyFilter);
        } else {
            comparison = new MCSFactory(bondTypeFlag, removeHydrogen, stereoFilter, fragmentFilter, energyFilter);
        }

        System.gc();

    }

    /**
     *
     * @param algorithmType 0: default, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
     * exhaustive MCS else false
     * @param bondTypeFlag   true will considered bond types for mapping else false
     * @param removeHydrogen true if Hydrogen are not to be mapped else false
     * @param stereoFilter   true if stereo match is considered else false
     * @param fragmentFilter true if fragement filter is switched on else false
     * @param energyFilter   true if bond energy filter is switched on else false
     * @throws Exception
     */
    public SMSD(int algorithmType, boolean bondTypeFlag, boolean removeHydrogen, boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) throws Exception {
        comparison = new MCSFactory(algorithmType, bondTypeFlag, removeHydrogen, stereoFilter, fragmentFilter, energyFilter);
        System.gc();

    }

    /**
     * 
     * @param Query
     * @param Target
     * @throws CDKException
     * 
     */
    @Override
    public void init(MolHandler Query, MolHandler Target) throws CDKException {

        if (Query.getMolecule().getAtomCount() > 0 && Target.getMolecule().getAtomCount() > 0) {
            comparison.init(Query, Target);
        } else {
            throw new CDKException("Each molecule should have atleast one atom to compare");
        }

        System.gc();

    }

    /**
     * 
     * @param Query
     * @param Target
     * @throws CDKException
     */
    @Override
    public synchronized void init(IMolecule Query, IMolecule Target) throws CDKException {

        if (Query.getAtomCount() > 0 && Target.getAtomCount() > 0) {
            comparison.init(Query, Target);
        } else {
            throw new CDKException("Each molecule should have atleast one atom to compare");
        }

        System.gc();
    }

    /**
     * 
     * @param Query
     * @param Target
     * @throws CDKException
     */
    @Override
    public synchronized void init(IAtomContainer Query, IAtomContainer Target) throws CDKException {

        if (Query.getAtomCount() > 0 && Target.getAtomCount() > 0) {
            comparison.init(Query, Target);
        } else {
            throw new CDKException("Each molecule should have atleast one atom to compare");
        }

        System.gc();

    }

    /**
     * 
     * @return All the possible mapped MCS atoms
     * List of mappings. Each map in the list has
     * atom-atom equivalence of the mappings
     * between query and target molecule i.e.
     * map.getKey() for the query and map.getValue()
     * for the target molecule
     */
    @Override
    public List<Map<IAtom, IAtom>> getAllAtomMapping() {
        return comparison.getAllAtomMapping();
    }

    /**
     *
     * @return All the possible mapped MCS atom index
     * List of mappings. Each map in the list has
     * atom-atom equivalence of the mappings
     * between query and target molecule i.e.
     * map.getKey() for the query and map.getValue()
     * for the target molecule
     */
    @Override
    public List<TreeMap<Integer, Integer>> getAllMapping() {
        return comparison.getAllMapping();
    }

    /**
     * 
     * @return One of the best MCS hits
     *
     */
    @Override
    public Map<IAtom, IAtom> getFirstAtomMapping() {
        return comparison.getFirstAtomMapping();
    }

    /**
     * 
     * @return One of the best MCS hits
     */
    @Override
    public TreeMap<Integer, Integer> getFirstMapping() {
        return comparison.getFirstMapping();
    }

    /**
     * 
     * @param Key position
     * @return total stereo matching 
     * score for the mapped part
     * 
     */
    @Override
    public Integer getStereoScore(int Key) {
        Integer solution = null;
        solution = comparison.getStereoScore(Key);
        if (solution == null) {
            solution = 0;
        }
        return solution;

    }

    /**
     * 
     * @param Key position
     * @return number of fragment(s) generated after removing the mapped parts
     * 
     */
    @Override
    public Integer getFragmentSize(int Key) {
        Integer solution = null;

        solution = comparison.getFragmentSize(Key);
        if (solution == null) {
            solution = 0;
        }
        return solution;
    }

    /**
     * 
     * @param Key position
     * @return total energy required to remove the mapped part
     *
     */
    @Override
    public Double getEnergyScore(int Key) {
        Double solution = 0.0;

        solution = comparison.getEnergyScore(Key);
        if (solution == null) {
            solution = 0.0;
        }
        return solution;
    }

    /**
     * 
     * @return Query Molecule AtomContainer
     */
    @Override
    public IAtomContainer getReactantMolecule() {
        return comparison.getReactantMolecule();
    }

    /**
     * 
     * @return Target Molecule AtomContainer
     */
    @Override
    public IAtomContainer getProductMolecule() {
        return comparison.getProductMolecule();
    }

    /**
     *
     * @return tanimoto similarity between query and target molecules
     * @throws java.io.IOException
     */
    @Override
    public double getTanimotoSimilarity() throws IOException {
        return comparison.getTanimotoSimilarity();
    }

    /**
     * 
     * @return true if mols have different stereo 
     * chemistry else flase if no stereo mismatch
     */
    @Override
    public boolean isStereoMisMatch() {
        return comparison.isStereoMisMatch();
    }

    /**
     * @return true if query is a subgraph of the target
     * else false
     */
    @Override
    public boolean isSubgraph() {
        return comparison.isSubgraph();

    }

    /**
     * @return Euclidean Distance between query and target molecule
     * @throws IOException
     */
    @Override
    public double getEuclideanDistance() throws IOException {
        return comparison.getEuclideanDistance();
    }
}
