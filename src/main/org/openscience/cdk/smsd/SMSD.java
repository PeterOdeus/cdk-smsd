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
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.smsd.factory.MCSFactory;
import org.openscience.cdk.smsd.factory.SubGraphFactory;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.IMCSAlgorithm;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;

/**
 *  This class implements the SMSD a multipurpose structure comparison tool.
 *  It allows users to find maximal common substructure (MCS), perform the
 *  mapping of a substructure in another structure, and the mapping of
 *  two isomorphic structures.
 *
 *  It also comes with various published algorithms and user is free to
 *  choose his favorite algorithm to perform MCS or substructure search.
 *  For example 0: SMSD algorithm, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
 *
 *  It also has a set of robust chemical filter (i.e. bond energy, fragment count,
 *  stereo & bond match) for sort out the reported MCS solution in a chemically
 *  relevant manner. Each comparison can be made with or without bond sensitive
 *  mode and with implicit or explicit hydrogens.
 *
 *  An example for MCS search:
 *  <pre>
 *  //Bond Sensitive is set true
 *  SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
 *  // acetic acid anhydride
 *  IAtomContainer A1 = sp.parseSmiles("CC");
 *  // acetic acid anhydride
 *  IAtomContainer A2 = sp.parseSmiles("CC(=O)OC(=O)C");
 *  SMSD comparison = new SMSD(3, true);
 *  // set molecules and remove hydrogens
 *  comparison.init(A1, A2, true);
 *  // set chemical filter true
 *  comparison.setChemFilters(true, true, true);
 *  System.out.println("Tanimoto coefficient:  " + comparison.getTanimotoSimilarity());
 *  System.out.println("A1 is a subgraph of A2:  " + comparison.isSubgraph());
 *  //Get Modified AtomContainer
 *  IAtomContainer Mol1 = comparison.getReactantMolecule();
 *  IAtomContainer Mol2 = comparison.getProductMolecule();
 *  // Print the mapping between molecules
 *  System.out.println(" Mappings: ");
 *  for (Map.Entry<Integer, Integer> mapping : comparison.getFirstMapping().entrySet()) {
 *      System.out.println((mapping.getKey() + 1) + " " + (mapping.getValue() + 1));
 *
 *      IAtom eAtom = Mol1.getAtom(mapping.getKey());
 *      IAtom pAtom = Mol2.getAtom(mapping.getValue());
 *      System.out.println(eAtom.getSymbol() + " " + pAtom.getSymbol());
 *  }
 *  System.out.println("");
 *  </pre>
 *
 *  @cdk.module smsd
 *  <p>SMSD algorithm is described in {@cdk.cite SMSD2009}.
 *  If you are using SMSD module</p>
 *  <p><font color="#FF0000">please cite {@cdk.cite SMSD2009}.</font>:
 *  </p>
 */

@TestClass("org.openscience.cdk.smsd.SMSDTest")
public class SMSD implements IMCSAlgorithm {

    IMCSAlgorithm comparison;

    /**
     * 
     * @param subStructureMode true for fast substructure search without
     * exhaustive MCS else false
     * @param bondTypeFlag   true will considered bond types for mapping else false
     */
    public SMSD(boolean subStructureMode, boolean bondTypeFlag) {
        if (subStructureMode) {
            comparison = new SubGraphFactory(bondTypeFlag);
        } else {
            comparison = new MCSFactory(bondTypeFlag);
        }
        System.gc();

    }

    /**
     *
     * @param algorithmType 0: default, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
     * exhaustive MCS else false
     * @param bondTypeFlag   true will considered bond types for mapping else false
     */
    public SMSD(int algorithmType, boolean bondTypeFlag) {
        comparison = new MCSFactory(algorithmType, bondTypeFlag);
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
    public void init(MolHandler Query, MolHandler Target, boolean removeHydrogen) throws CDKException {

        if (Query.getMolecule().getAtomCount() > 0 && Target.getMolecule().getAtomCount() > 0) {
            comparison.init(Query, Target, removeHydrogen);
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
    public synchronized void init(IMolecule Query, IMolecule Target, boolean removeHydrogen) throws CDKException {

        if (Query.getAtomCount() > 0 && Target.getAtomCount() > 0) {
            comparison.init(Query, Target, removeHydrogen);
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
    public synchronized void init(IAtomContainer Query, IAtomContainer Target, boolean removeHydrogen) throws CDKException {

        if (Query.getAtomCount() > 0 && Target.getAtomCount() > 0) {
            comparison.init(Query, Target, removeHydrogen);
        } else {
            throw new CDKException("Each molecule should have atleast one atom to compare");
        }

        System.gc();

    }

    public void init(String sourceMolFileName, String targetMolFileName, boolean removeHydrogen) throws CDKException {

        String mol1 = sourceMolFileName;
        String mol2 = targetMolFileName;

        MolHandler Query = new MolHandler(mol1, false);
        MolHandler Target = new MolHandler(mol2, false);

        init(Query, Target, removeHydrogen);
        System.gc();
    }

    /*
     * @param stereoFilter   true if stereo match is considered else false
     * @param fragmentFilter true if fragement filter is switched on else false
     * @param energyFilter   true if bond energy filter is switched on else false
     */
    @Override
    public void setChemFilters(boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) {
        comparison.setChemFilters(stereoFilter, fragmentFilter, energyFilter);
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
        return comparison.getStereoScore(Key);
    }

    /**
     * 
     * @param Key position
     * @return number of fragment(s) generated after removing the mapped parts
     * 
     */
    @Override
    public Integer getFragmentSize(int Key) {
        return comparison.getFragmentSize(Key);
    }

    /**
     * 
     * @param Key position
     * @return total energy required to remove the mapped part
     *
     */
    @Override
    public Double getEnergyScore(int Key) {
        return comparison.getEnergyScore(Key);
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
