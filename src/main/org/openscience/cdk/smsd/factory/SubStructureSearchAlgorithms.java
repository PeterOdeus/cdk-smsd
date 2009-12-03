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
 * You should have received rAtomCount copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.smsd.factory;

import org.openscience.cdk.smsd.interfaces.IMCSAlgorithm;
import org.openscience.cdk.smsd.interfaces.IMCS;
import org.openscience.cdk.smsd.algorithm.cdk.CDKMCSHandler;
import org.openscience.cdk.smsd.algorithm.mcsplus.MCSPlusHandler;
import org.openscience.cdk.smsd.algorithm.vflib.VFlibMCSHandler;
import org.openscience.cdk.smsd.algorithm.single.SingleMappingHandler;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.smsd.filters.ChemicalFilters;
import org.openscience.cdk.smsd.global.BondType;
import org.openscience.cdk.smsd.global.TimeOut;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smsd.algorithm.vflib.VFlibTurboHandler;

/**
 * @cdk.module smsd
 */
public class SubStructureSearchAlgorithms implements IMCS {

    private List<TreeMap<Integer, Integer>> allMCS = null;
    private TreeMap<Integer, Integer> firstSolution = null;
    private List<Map<IAtom, IAtom>> allAtomMCS = null;
    private Map<IAtom, IAtom> firstAtomMCS = null;
    private MolHandler RMol = null;
    private MolHandler PMol = null;
    private IAtomContainerSet RFrag = null;
    private IAtomContainerSet PFrag = null;
    private List<Double> stereoScore = null;
    private List<Integer> fragmentSize = null;
    private List<Double> bEnergies = null;
    private IMCSAlgorithm mcs = null;
    private Algorithm algorithmType;
    private boolean removeHydrogen = false;

    /**
     * 
     * @param algorithmType 0 default, 1 mcsPlus, 2 VFLib, 3 CDKMCS 4 SubStructure
     * @param bondTypeFlag
     */
    public SubStructureSearchAlgorithms(Algorithm algorithmType, boolean bondTypeFlag) {
        this.algorithmType = algorithmType;
        firstSolution = new TreeMap<Integer, Integer>();
        allMCS = new ArrayList<TreeMap<Integer, Integer>>();
        allAtomMCS = new ArrayList<Map<IAtom, IAtom>>();
        firstAtomMCS = new HashMap<IAtom, IAtom>();

        if (bondTypeFlag) {
            TimeOut tmo = TimeOut.getInstance();
            tmo.setTimeOut(0.10);
        } else {
            TimeOut tmo = TimeOut.getInstance();
            tmo.setTimeOut(0.15);
        }

        BondType bondType = BondType.getInstance();
        bondType.reset();
        bondType.setBondSensitiveFlag(bondTypeFlag);
    }

    private synchronized void mcsBuilder() {

        int rBondCount = RMol.getMolecule().getBondCount();
        int pBondCount = PMol.getMolecule().getBondCount();

        int rAtomCount = RMol.getMolecule().getAtomCount();
        int pAtomCount = PMol.getMolecule().getAtomCount();
        if (rBondCount == 0 || rAtomCount == 1 || pBondCount == 0 || pAtomCount == 1) {
            singleMapping();
        } else {

            switch (algorithmType) {
                case CDKMCS:
                    cdkMCS();
                    if (getFirstMapping() == null) {
                        mcs = null;
                        System.gc();
                    }
                    break;
                case DEFAULT:
                    if (BondType.getInstance().getBondSensitiveFlag()) {
                        cdkMCS();
                    } else {
                        mcsPlus();
                    }
                    if (getFirstMapping() == null) {
                        mcs = null;
                        System.gc();
                        vfLibMCS();
                    }
                    break;
                case MCSPlus:
                    mcsPlus();
                    if (getFirstMapping() == null) {
                        mcs = null;
                    }
                    break;
                case SubStructure:
                    if (rBondCount > 1 && pBondCount > 1) {
                        vfTurboHandler();
                        if (getFirstMapping() == null) {
                            mcs = null;
                            System.gc();
                        }
                    } else {
                        singleMapping();
                    }
                    System.gc();
                    break;
                case VFLibMCS:
                    if (rBondCount >= 6 && rBondCount >= 6) {
                        vfLibMCS();
                        if (getFirstMapping() == null) {
                            mcs = null;
                            System.gc();
                        }
                    } else {
                        mcsPlus();
                    }
                    break;

            }

        }
        System.gc();

    }

    private synchronized void fragmentBuilder() {

        FragmentMatcher fragmentMatcher = new FragmentMatcher(RFrag, PFrag, removeHydrogen);
        fragmentMatcher.searchMCS();

        firstSolution.clear();
        allMCS.clear();
        allAtomMCS.clear();
        firstAtomMCS.clear();
        firstSolution.putAll(fragmentMatcher.getFirstMapping());
        allMCS.addAll(fragmentMatcher.getAllMapping());

        firstAtomMCS.putAll(fragmentMatcher.getFirstAtomMapping());
        allAtomMCS.addAll(fragmentMatcher.getAllAtomMapping());

    }

    private synchronized void cdkMCS() {
        try {
            mcs = new CDKMCSHandler();

            mcs.set(RMol, PMol);
            mcs.searchMCS();

            firstSolution.clear();
            allMCS.clear();
            allAtomMCS.clear();
            firstAtomMCS.clear();

            firstSolution.putAll(mcs.getFirstMapping());
            allMCS.addAll(mcs.getAllMapping());

            firstAtomMCS.putAll(mcs.getFirstAtomMapping());
            allAtomMCS.addAll(mcs.getAllAtomMapping());


        } catch (IOException ex) {
            Logger.getLogger(SubStructureSearchAlgorithms.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(SubStructureSearchAlgorithms.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private synchronized void mcsPlus() {
        try {
            mcs = new MCSPlusHandler();

            mcs.set(RMol, PMol);
            mcs.searchMCS();

            firstSolution.clear();
            allMCS.clear();
            allAtomMCS.clear();
            firstAtomMCS.clear();

            firstSolution.putAll(mcs.getFirstMapping());
            allMCS.addAll(mcs.getAllMapping());

            firstAtomMCS.putAll(mcs.getFirstAtomMapping());
            allAtomMCS.addAll(mcs.getAllAtomMapping());


        } catch (IOException ex) {
            Logger.getLogger(SubStructureSearchAlgorithms.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(SubStructureSearchAlgorithms.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private void vfTurboHandler() {

        VFlibTurboHandler subGraphTurboSearch = new VFlibTurboHandler();
        subGraphTurboSearch.set(RMol, PMol);

        firstSolution.clear();
        allMCS.clear();
        allAtomMCS.clear();
        firstAtomMCS.clear();

        if (subGraphTurboSearch.isSubgraph()) {
            firstSolution.putAll(subGraphTurboSearch.getFirstMapping());
            allMCS.addAll(subGraphTurboSearch.getAllMapping());
            firstAtomMCS.putAll(subGraphTurboSearch.getFirstAtomMapping());
            allAtomMCS.addAll(subGraphTurboSearch.getAllAtomMapping());
        }
    }

    /**
     *
     * @param Reactant
     * @param Product
     * @param removeHydrogen
     *
     */
    private void init(MolHandler Reactant, MolHandler Product, boolean removeHydrogen) {
        this.removeHydrogen = removeHydrogen;
        this.RMol = new MolHandler(Reactant.getMolecule(), false, removeHydrogen);
        this.PMol = new MolHandler(Product.getMolecule(), false, removeHydrogen);

        if (RMol.getConnectedFlag() && PMol.getConnectedFlag()) {
            mcsBuilder();
        } else {
            this.RFrag = RMol.getFragmentedMolecule();
            this.PFrag = PMol.getFragmentedMolecule();
            fragmentBuilder();
        }

    }

    /**
     *
     * @param Reactant
     * @param Product
     */
    @Override
    public synchronized void init(IMolecule Reactant, IMolecule Product, boolean removeHydrogen) {
        this.removeHydrogen = removeHydrogen;
        this.RMol = new MolHandler(Reactant, false, removeHydrogen);
        this.PMol = new MolHandler(Product, false, removeHydrogen);
        init(RMol, PMol, removeHydrogen);
    }

    /**
     *
     * @param Reactant
     * @param Product
     */
    @Override
    public synchronized void init(IAtomContainer Reactant, IAtomContainer Product, boolean removeHydrogen) {
        this.removeHydrogen = removeHydrogen;
        this.RMol = new MolHandler(Reactant, false, removeHydrogen);
        this.PMol = new MolHandler(Product, false, removeHydrogen);
        init(RMol, PMol, removeHydrogen);
    }

    public void init(String sourceMolFileName, String targetMolFileName, boolean removeHydrogen) throws CDKException {
        String mol1 = sourceMolFileName;
        String mol2 = targetMolFileName;

        this.removeHydrogen = removeHydrogen;
        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);
        init(Reactant, Product, removeHydrogen);
    }

    public void setChemFilters(boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) {

        if (firstAtomMCS != null) {
            ChemicalFilters chemFilter = new ChemicalFilters(allMCS, allAtomMCS, firstSolution, firstAtomMCS, RMol, PMol);

            if (stereoFilter) {
                chemFilter.sortResultsByStereoAndBondMatch();
            }
            if (fragmentFilter) {
                chemFilter.sortResultsByFragments();
            }

            if (energyFilter) {
                try {
                    chemFilter.sortResultsByEnergies();
                } catch (CDKException ex) {
                    Logger.getLogger(SubStructureSearchAlgorithms.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

            this.stereoScore = chemFilter.getStereoMatches();
            this.fragmentSize = chemFilter.getSortedFragment();
            this.bEnergies = chemFilter.getSortedEnergy();
        }
    }

    @Override
    public synchronized Integer getFragmentSize(int Key) {
        return fragmentSize.get(Key) == null ? null : fragmentSize.get(Key);
    }

    @Override
    public Integer getStereoScore(int Key) {
        return stereoScore.get(Key) == null ? null : stereoScore.get(Key).intValue();
    }

    @Override
    public Double getEnergyScore(int Key) {
        return bEnergies.get(Key) == null ? null : bEnergies.get(Key);
    }

    @Override
    public synchronized TreeMap<Integer, Integer> getFirstMapping() {
        return firstSolution.isEmpty() ? null : firstSolution;
    }

    @Override
    public synchronized List<TreeMap<Integer, Integer>> getAllMapping() {
        return allMCS.isEmpty() ? null : allMCS;
    }

    @Override
    public synchronized Map<IAtom, IAtom> getFirstAtomMapping() {
        return firstAtomMCS.isEmpty() ? null : firstAtomMCS;
    }

    @Override
    public synchronized List<Map<IAtom, IAtom>> getAllAtomMapping() {
        return allAtomMCS.isEmpty() ? null : allAtomMCS;
    }

    @Override
    public IAtomContainer getReactantMolecule() {
        return RMol.getMolecule();
    }

    @Override
    public IAtomContainer getProductMolecule() {
        return PMol.getMolecule();
    }

    @Override
    public double getTanimotoSimilarity() throws IOException {
        int decimalPlaces = 4;
        int rAtomCount = 0;
        int pAtomCount = 0;
        if (!removeHydrogen) {
            rAtomCount = RMol.getMolecule().getAtomCount();
            pAtomCount = PMol.getMolecule().getAtomCount();
        } else {
            rAtomCount = RMol.getMolecule().getAtomCount() - getHCount(RMol.getMolecule());
            pAtomCount = PMol.getMolecule().getAtomCount() - getHCount(PMol.getMolecule());
        }
        double matchCount = getFirstMapping().size();
        double tanimoto = (matchCount) / (rAtomCount + pAtomCount - matchCount);

        BigDecimal tan = new BigDecimal(tanimoto);

        tan = tan.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
        tanimoto = tan.doubleValue();
        return tanimoto;
    }

    /**
     *
     * @return true if mols have different stereo
     * chemistry else flase if no stereo mismatch
     */
    @Override
    public boolean isStereoMisMatch() {
        boolean flag = false;
        IAtomContainer Reactant = RMol.getMolecule();
        IAtomContainer Product = PMol.getMolecule();
        int Score = 0;

        for (Map.Entry<IAtom, IAtom> mappingI : firstAtomMCS.entrySet()) {
            IAtom indexI = mappingI.getKey();
            IAtom indexJ = mappingI.getValue();
            for (Map.Entry<IAtom, IAtom> mappingJ : firstAtomMCS.entrySet()) {

                IAtom indexIPlus = mappingJ.getKey();
                IAtom indexJPlus = mappingJ.getValue();
                if (!indexI.equals(indexIPlus) && !indexJ.equals(indexJPlus)) {

                    IAtom sourceAtom1 = indexI;
                    IAtom sourceAtom2 = indexIPlus;

                    IBond RBond = Reactant.getBond(sourceAtom1, sourceAtom2);

                    IAtom targetAtom1 = indexJ;
                    IAtom targetAtom2 = indexJPlus;
                    IBond PBond = Product.getBond(targetAtom1, targetAtom2);

                    if ((RBond != null && PBond != null) && (RBond.getStereo() != PBond.getStereo())) {
                        Score++;
                    }
                }

            }
        }
        if (Score > 0) {
            flag = true;
        }
        return flag;
    }

    /**
     * @return true if ac1 is rAtomCount subgraph of ac2
     */
    @Override
    public boolean isSubgraph() {

        IAtomContainer Reactant = RMol.getMolecule();
        IAtomContainer Product = PMol.getMolecule();
        if (firstAtomMCS == null || firstAtomMCS.isEmpty()) {
            return false;
        }
        int score = getBondMatchScore(Reactant, Product);

        boolean flag = false;
        float size = firstSolution.size();
        int source = 0;
        int target = 0;
        if (!removeHydrogen) {
            source = Reactant.getAtomCount();
            target = Product.getAtomCount();
        } else {
            source = Reactant.getAtomCount() - getHCount(Reactant);
            target = Product.getAtomCount() - getHCount(Product);
        }
        if ((size == source && score / 2 == Reactant.getBondCount()) &&
                (target >= size && Product.getBondCount() >= score / 2)) {
            flag = true;
        }
        return flag;
    }

    @Override
    public double getEuclideanDistance() throws IOException {
        int decimalPlaces = 4;
        double source = 0;
        double target = 0;
        if (!removeHydrogen) {
            source = RMol.getMolecule().getAtomCount();
            target = PMol.getMolecule().getAtomCount();
        } else {
            source = RMol.getMolecule().getAtomCount() - getHCount(RMol.getMolecule());
            target = PMol.getMolecule().getAtomCount() - getHCount(PMol.getMolecule());
        }
        double common = getFirstMapping().size();

        double euclidean = Math.sqrt(source + target - 2 * common);

        BigDecimal dist = new BigDecimal(euclidean);

        dist = dist.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
        euclidean = dist.doubleValue();
        return euclidean;
    }

    private void vfLibMCS() {
        try {

            mcs = new VFlibMCSHandler();
            mcs.set(RMol, PMol);
            mcs.searchMCS();

            firstSolution.clear();
            allMCS.clear();
            allAtomMCS.clear();
            firstAtomMCS.clear();

            firstSolution.putAll(mcs.getFirstMapping());
            allMCS.addAll(mcs.getAllMapping());

            firstAtomMCS.putAll(mcs.getFirstAtomMapping());
            allAtomMCS.addAll(mcs.getAllAtomMapping());


        } catch (IOException ex) {
            Logger.getLogger(SubStructureSearchAlgorithms.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(SubStructureSearchAlgorithms.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private void singleMapping() {
        try {

            mcs = new SingleMappingHandler(removeHydrogen);
            mcs.set(RMol, PMol);
            mcs.searchMCS();

            firstSolution.clear();
            allMCS.clear();
            allAtomMCS.clear();
            firstAtomMCS.clear();

            firstSolution.putAll(mcs.getFirstMapping());
            allMCS.addAll(mcs.getAllMapping());

            firstAtomMCS.putAll(mcs.getFirstAtomMapping());
            allAtomMCS.addAll(mcs.getAllAtomMapping());


        } catch (IOException ex) {
            Logger.getLogger(SubStructureSearchAlgorithms.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(SubStructureSearchAlgorithms.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private int getHCount(IAtomContainer molecule) {
        int count = 0;
        for (IAtom atom : molecule.atoms()) {
            if (atom.getSymbol().equalsIgnoreCase("H")) {
                ++count;
            }
        }

        return count;
    }

    private int getBondMatchScore(IBond RBond, IBond PBond, BondType bondType) {
        int score = 0;
        if (RBond != null && PBond != null) {
            if (bondType.getBondSensitiveFlag()) {
                int RBondType = RBond.getOrder().ordinal();
                int PBondType = PBond.getOrder().ordinal();

                if (RBond.getFlag(CDKConstants.ISAROMATIC) == PBond.getFlag(CDKConstants.ISAROMATIC) && RBondType == PBondType) {
                    score++;
                } else if (RBond.getFlag(CDKConstants.ISAROMATIC) && PBond.getFlag(CDKConstants.ISAROMATIC)) {
                    score++;
                }
            } else {
                score++;
            }
        }
        return score;
    }

    private int getBondMatchScore(IAtomContainer Reactant, IAtomContainer Product) {
        int score = 0;

        BondType bondType = BondType.getInstance();
        for (Map.Entry<IAtom, IAtom> mappingI : firstAtomMCS.entrySet()) {
            IAtom indexI = mappingI.getKey();
            IAtom indexJ = mappingI.getValue();
            for (Map.Entry<IAtom, IAtom> mappingJ : firstAtomMCS.entrySet()) {

                IAtom indexIPlus = mappingJ.getKey();
                IAtom indexJPlus = mappingJ.getValue();
                if (!indexI.equals(indexIPlus) && !indexJ.equals(indexJPlus)) {

                    IAtom sourceAtom1 = indexI;
                    IAtom sourceAtom2 = indexIPlus;
                    IBond RBond = Reactant.getBond(sourceAtom1, sourceAtom2);

                    IAtom targetAtom1 = indexJ;
                    IAtom targetAtom2 = indexJPlus;
                    IBond PBond = Product.getBond(targetAtom1, targetAtom2);

                    score += getBondMatchScore(RBond, PBond, bondType);
                }
            }
        }
        return score;
    }
}

