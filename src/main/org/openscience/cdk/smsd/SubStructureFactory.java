/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.smsd.core.tools.EBIException;
import org.openscience.cdk.smsd.factory.MCSFactory;
import org.openscience.cdk.smsd.factory.SubGraphFactory;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.IMCSAlgorithm;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;

/**
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class SubStructureFactory implements IMCSAlgorithm {

    IMCSAlgorithm comparison;

    /**
     * 
     * @param subStructureMode
     * @param bondTypeFlag This flag if set to true will considered bond
     * types for mapping. Hence only similar bon types will be mapped.
     * @param removeHydrogen
     * @param stereoFilter
     * @param fragmentFilter
     * @param energyFilter
     * @throws Exception
     */
    public SubStructureFactory(boolean subStructureMode, boolean bondTypeFlag, boolean removeHydrogen, boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) throws Exception {

//        this.bondFlag = bondTypeFlag;

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
     * @param subStructureMode
     * @param bondTypeFlag This flag if set to true will considered bond
     * types for mapping. Hence only similar bon types will be mapped.
     * @param removeHydrogen
     * @param stereoFilter
     * @param fragmentFilter
     * @param energyFilter
     * @throws Exception
     */
    public SubStructureFactory(int algorithmType, boolean subStructureMode, boolean bondTypeFlag, boolean removeHydrogen, boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) throws Exception {

//        this.bondFlag = bondTypeFlag;

        if (subStructureMode) {

            comparison = new SubGraphFactory(bondTypeFlag, removeHydrogen, stereoFilter, fragmentFilter, energyFilter);

        } else {
            comparison = new MCSFactory(algorithmType, bondTypeFlag, removeHydrogen, stereoFilter, fragmentFilter, energyFilter);

        }

        System.gc();

    }

    /**
     * 
     * @param Reactant
     * @param Product
     * @throws EBIException 
     * 
     */
    @Override
    public void init(MolHandler Reactant, MolHandler Product) throws EBIException {

        if (Reactant.getMolecule().getAtomCount() > 0 && Product.getMolecule().getAtomCount() > 0) {

            comparison.init(Reactant, Product);

        } else {
            throw new EBIException("Each molecule should have atleast one atom to compare");
        }

        System.gc();

    }

    /**
     * 
     * @param Reactant
     * @param Product
     * @throws EBIException 
     */
    @Override
    public synchronized void init(IMolecule Reactant, IMolecule Product) throws EBIException {

        if (Reactant.getAtomCount() > 0 && Product.getAtomCount() > 0) {
//           
            comparison.init(Reactant, Product);

        } else {
            throw new EBIException("Each molecule should have atleast one atom to compare");
        }

        System.gc();
    }

    /**
     * 
     * @param Reactant
     * @param Product
     * @throws EBIException 
     */
    @Override
    public synchronized void init(IAtomContainer Reactant, IAtomContainer Product) throws EBIException {

        if (Reactant.getAtomCount() > 0 && Product.getAtomCount() > 0) {
//            System.out.println("SubStructureFactory");
//            printMolecules(Reactant, Product);
            comparison.init(Reactant, Product);

        } else {
            throw new EBIException("Each molecule should have atleast one atom to compare");
        }

        System.gc();

    }

    /**
     * 
     * @return All the possible MCS Atoms
     * List of mapped position
     * i.e if i=0 is reactant mapping index,
     * i+1 is product index
     */
    @Override
    public List<Map<IAtom, IAtom>> getAllAtomMapping() {
        return comparison.getAllAtomMapping();
    }

    /**
     * 
     * @return All the possible MCS mapping
     * List of mapped position
     * i.e if i=0 is reactant mapping index,
     * i+1 is product index
     */
    @Override
    public List<TreeMap<Integer, Integer>> getAllMapping() {
        return comparison.getAllMapping();
    }

    /**
     * 
     * @return MCS list of mapped Atoms
     * i.e if i=0 is reactant mapping index,
     * i+1 is product index
     */
    @Override
    public Map<IAtom, IAtom> getFirstAtomMapping() {
        return comparison.getFirstAtomMapping();
    }

    /**
     * 
     * @return MCS list of mapped positions
     * i.e if i=0 is reactant mapping index,
     * i+1 is product index
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

//        System.out.println("Key: " + Key);
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
    public Integer getFragmentSize(
            int Key) {
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
    public Double getEnergyScore(
            int Key) {
        Double solution = 0.0;

        solution = comparison.getEnergyScore(Key);
        if (solution == null) {
            solution = 0.0;
        }
        return solution;
    }

    /**
     * 
     * @return Reactant Molecule AtomContainer
     */
    @Override
    public IAtomContainer getReactantMolecule() {
        IAtomContainer solution;

        solution = comparison.getReactantMolecule();

        return solution;
    }

    /**
     * 
     * @return Product Molecule AtomContainer
     */
    @Override
    public IAtomContainer getProductMolecule() {
        IAtomContainer solution;

        solution = comparison.getProductMolecule();

        return solution;
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
     * 
     * @param Molecule1
     * @param Molecule2
     */
    private void printMolecules(IAtomContainer Molecule1, IAtomContainer Molecule2) {

        System.out.println("Molecule 1");

        for (int i = 0; i < Molecule1.getAtomCount(); i++) {

            System.out.print(Molecule1.getAtom(i).getSymbol() + " ");
        }

        System.out.println();
        System.out.println("Molecule 2");
        for (int i = 0; i < Molecule2.getAtomCount(); i++) {

            System.out.print(Molecule2.getAtom(i).getSymbol() + " ");
        }
        System.out.println();

    }

    /**
     * @return true if query is a subgraph of the target
     * else false
     */
    @Override
    public boolean isSubgraph() {
        return comparison.isSubgraph();

    }

    @Override
    public double getEuclideanDistance() throws IOException {
        return comparison.getEuclideanDistance();
    }
}
