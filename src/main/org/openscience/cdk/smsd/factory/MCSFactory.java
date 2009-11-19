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

import org.openscience.cdk.smsd.interfaces.IMCS;
import org.openscience.cdk.smsd.interfaces.IMCSAlgorithm;
import org.openscience.cdk.smsd.interfaces.IMCSBase;
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
import org.openscience.cdk.smsd.filters.FragmentMatcher;
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

/**
 * @cdk.module smsd
 */
public class MCSFactory implements IMCSAlgorithm {

    private List<TreeMap<Integer, Integer>> allMCS = null;
    private TreeMap<Integer, Integer> firstSolution = null;
    private List<Map<IAtom, IAtom>> allAtomMCS = null;
    private Map<IAtom, IAtom> firstAtomMCS = null;
    private MolHandler RMol = null;
    private MolHandler PMol = null;
    private IAtomContainerSet RFrag = null;
    private IAtomContainerSet PFrag = null;
    private List<Integer> StereoScore = null;
    private List<Integer> fragmentSize = null;
    private List<Double> bEnergies = null;
    private boolean removeHydrogen = false;
    private boolean stereoFilter = false;
    private boolean fragmentFilter = false;
    private boolean energyFilter = false;
    private IMCS mcs = null;
    private int algorithmType = 0;

    /**
     * 
     * @param algorithmType 0 default, 1 mcsPlus, 2 VFLib, 3 cdkMCS
     * @param bondTypeFlag
     * @param removeHydrogen
     * @param stereoFilter
     * @param fragmentFilter
     * @param energyFilter
     */
    public MCSFactory(
            int algorithmType,
            boolean bondTypeFlag,
            boolean removeHydrogen,
            boolean stereoFilter,
            boolean fragmentFilter,
            boolean energyFilter){

        this.removeHydrogen = removeHydrogen;
        this.stereoFilter = stereoFilter;
        this.fragmentFilter = fragmentFilter;
        this.energyFilter = energyFilter;
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

        BondType BT = BondType.getInstance();
        BT.reset();
        BT.setBondSensitiveFlag(bondTypeFlag);


    }

    public MCSFactory(
            boolean bondTypeFlag,
            boolean removeHydrogen,
            boolean stereoFilter,
            boolean fragmentFilter,
            boolean energyFilter) {

        this.removeHydrogen = removeHydrogen;
        this.stereoFilter = stereoFilter;
        this.fragmentFilter = fragmentFilter;
        this.energyFilter = energyFilter;


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

//            int commonAtoms = checkCommonAtomCount(RMol.getMolecule(), PMol.getMolecule());

            /*This is importsnt because CDK fails to Generate makeAtomsMapOfBondsMap if bonds are less than 2*/

//            long startTime = System.currentTimeMillis();
//
            if (rBondCount == 0 || rAtomCount == 1 || pBondCount == 0 || pAtomCount == 1) {
//                System.out.println("Single Mapping");
                singleMapping();
            } else if (algorithmType == 0) {
//                System.err.println("Default");
//                System.out.println("mcsBuilder");
//                printMolecules(RMol.getMolecule(), PMol.getMolecule());

                mcsPlus();
                if (getFirstMapping() == null) {
                    mcs = null;
                    System.gc();
//                    System.err.println("mcsPlus timeout ");
                    vfLibMCS();
//                    System.out.println("Mapped with VF-McGregor");
                }
            } else if (algorithmType == 1) {
//                System.err.println("mcsPlus");
//                System.out.println("mcsBuilder");
//                printMolecules(RMol.getMolecule(), PMol.getMolecule());
                mcsPlus();
                if (getFirstMapping() == null) {
                    mcs = null;
//                    System.out.println("Time-out occured mcsPlus");
//                    System.gc();
                }

            } else if (algorithmType == 2) {
//                System.err.println("vfLibMCS");
//                System.out.println("mcsBuilder");
//                printMolecules(RMol.getMolecule(), PMol.getMolecule());

                if (rBondCount >= 6 && rBondCount >= 6) {
                    vfLibMCS();
                    if (getFirstMapping() == null) {
                        mcs = null;
                        System.gc();
                    }
//                    System.out.println("Mapped with vfLibMCS");
                } else {
                    mcsPlus();
//                    System.out.println("Mapped with mcsPlus");
                }
            } else if (algorithmType == 3) {
//                 System.err.println("cdkMCS");
                cdkMCS();
                if (getFirstMapping() == null) {
//                    System.out.println("Time-out occured mcsPlus");
                    mcs = null;
                    System.gc();
                }

//              System.out.println("Mapped with cdkMCS");

            }
//
////
//            System.out.println("MCS solution count:" + this.allMCS.size());
//            System.out.println("solution Size:" + this.allMCS.firstElement().size());
//            for (Map.Entry<IAtom, IAtom> map : allAtomMCS.firstElement().entrySet()) {
//                System.out.println(map.getKey().getSymbol() + ":" + map.getValue().getSymbol());
//                System.out.println(map.getKey().getID() + ":" + map.getValue().getID());
//            }
//            long endTime = System.currentTimeMillis();
            System.gc();
//            System.out.println("Calculation Time: " + (endTime - startTime) * 0.001 + " seconds");
//            System.out.println("!Done!\n");


    }

    private synchronized void fragmentBuilder() {

        IMCSBase fragmentMatcher = new FragmentMatcher(RFrag, PFrag, removeHydrogen);


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
            mcs.searchMCS(removeHydrogen);

            firstSolution.clear();
            allMCS.clear();
            allAtomMCS.clear();
            firstAtomMCS.clear();


            firstSolution.putAll(mcs.getFirstMapping());
            allMCS.addAll(mcs.getAllMapping());

            firstAtomMCS.putAll(mcs.getFirstAtomMapping());
            allAtomMCS.addAll(mcs.getAllAtomMapping());


        } catch (IOException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private synchronized void mcsPlus() {
        try {
            mcs = new MCSPlusHandler();

            mcs.set(RMol, PMol);
            mcs.searchMCS(removeHydrogen);

            firstSolution.clear();
            allMCS.clear();
            allAtomMCS.clear();
            firstAtomMCS.clear();

            firstSolution.putAll(mcs.getFirstMapping());
            allMCS.addAll(mcs.getAllMapping());

            firstAtomMCS.putAll(mcs.getFirstAtomMapping());
            allAtomMCS.addAll(mcs.getAllAtomMapping());


        } catch (IOException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    /**
     *
     * @param Reactant
     * @param Product
     *
     */
    @Override
    public void init(MolHandler Reactant, MolHandler Product) {
        try {

            this.RMol = new MolHandler(Reactant.getMolecule(), false, removeHydrogen);
            this.PMol = new MolHandler(Product.getMolecule(), false, removeHydrogen);

            if (RMol.getConnectedFlag() && PMol.getConnectedFlag()) {
                mcsBuilder();
            } else {
                this.RFrag = RMol.getFragmentedMolecule();
                this.PFrag = PMol.getFragmentedMolecule();
                fragmentBuilder();
            }
            setChemFilters();
        } catch (CDKException ex) {
            Logger.getLogger(SubGraphFactory.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     *
     * @param Reactant
     * @param Product
     */
    @Override
    public synchronized void init(IMolecule Reactant, IMolecule Product) {
        this.RMol = new MolHandler(Reactant, false, removeHydrogen);
        this.PMol = new MolHandler(Product, false, removeHydrogen);
        init(RMol, PMol);
    }

    /**
     *
     * @param Reactant
     * @param Product
     */
    @Override
    public synchronized void init(IAtomContainer Reactant, IAtomContainer Product) {
        this.RMol = new MolHandler(Reactant, false, removeHydrogen);
        this.PMol = new MolHandler(Product, false, removeHydrogen);
        init(RMol, PMol);
    }

    public synchronized void setChemFilters() throws CDKException {
        if (firstAtomMCS != null) {
            ChemicalFilters chemFilter = new ChemicalFilters(allMCS, allAtomMCS, firstSolution, firstAtomMCS, RMol, PMol);

            if (stereoFilter) {
                chemFilter.sortResultsByStereoAndBondMatch();
            }
            if (fragmentFilter) {
                chemFilter.sortResultsByFragments();
            }

            if (energyFilter) {
                chemFilter.sortResultsByEnergies();
            }

            this.StereoScore = chemFilter.getStereoMatches();
            this.fragmentSize = chemFilter.getSortedFragment();
            this.bEnergies = chemFilter.getSortedEnergy();
        }
    }

    @Override
    public synchronized Integer getFragmentSize(int Key) {
        Integer Value = null;
        // System.out.println("Key" + Key);
        if (fragmentSize != null && firstSolution.size() > 0 &&
                fragmentSize.size() > Key && Key >= 0) {
            Value = fragmentSize.get(Key);
        }
        return Value;

    }

    @Override
    public Integer getStereoScore(int Key) {

        Integer Value = null;
        if (StereoScore != null && StereoScore.size() > Key && Key >= 0) {
            Value = StereoScore.get(Key);
        }

        return Value;
    }

    @Override
    public Double getEnergyScore(int Key) {

        Double Value = null;
        if (bEnergies != null && bEnergies.size() > Key &&
                Key >= 0) {
            Value = bEnergies.get(Key);
        }
        return Value;
    }

    /**
     *
     * @return
     */
    @Override
    public synchronized TreeMap<Integer, Integer> getFirstMapping() {
        if (firstSolution.size() > 0) {
            return firstSolution;
        } else {
            return null;
        }

    }

    /**
     *
     * @return
     */
    @Override
    public synchronized List<TreeMap<Integer, Integer>> getAllMapping() {
        if (allMCS.size() > 0) {
            return allMCS;
        } else {
            return null;
        }

    }

    /**
     *
     * @return
     */
    @Override
    public synchronized Map<IAtom, IAtom> getFirstAtomMapping() {
        if (firstSolution.size() > 0) {
            return firstAtomMCS;
        } else {

            return null;
        }

    }

    /**
     *
     * @return
     */
    @Override
    public synchronized List<Map<IAtom, IAtom>> getAllAtomMapping() {
        if (allMCS.size() > 0) {
            return allAtomMCS;
        } else {
            try {
                throw new CDKException("No Solutions!");
            } catch (CDKException ex) {
                Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
            }
            return null;
        }

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

    @Override
    /**
     *
     * @return true if mols have different stereo
     * chemistry else flase if no stereo mismatch
     */
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
                if (indexI.equals(indexIPlus) && indexJ.equals(indexJPlus)) {

                    IAtom sourceAtom1 = indexI;
                    IAtom sourceAtom2 = indexIPlus;

                    IBond RBond = Reactant.getBond(sourceAtom1, sourceAtom2);

                    if (RBond != null) {

                        IAtom targetAtom1 = indexJ;
                        IAtom targetAtom2 = indexJPlus;
                        IBond PBond = Product.getBond(targetAtom1, targetAtom2);

                        if ((PBond != null) && (RBond.getStereo() != PBond.getStereo())) {
                            Score++;
                        }
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
     * @return true if ac1 is sourceAtomCount subgraph of ac2
     */
    @Override
    public boolean isSubgraph() {
        if (firstAtomMCS == null || firstAtomMCS.isEmpty()) {

            return false;
        }
        BondType bondType = BondType.getInstance();
        int score = 0;
        for (Map.Entry<IAtom, IAtom> mappingI : firstAtomMCS.entrySet()) {
            IAtom indexI = mappingI.getKey();
            IAtom indexJ = mappingI.getValue();
            for (Map.Entry<IAtom, IAtom> mappingJ : firstAtomMCS.entrySet()) {

                IAtom indexIPlus = mappingJ.getKey();
                IAtom indexJPlus = mappingJ.getValue();
                if (!indexI.equals(indexIPlus) && !indexJ.equals(indexJPlus)) {

                    IAtom sourceAtom1 = indexI;
                    IAtom sourceAtom2 = indexIPlus;

                    IBond RBond = RMol.getMolecule().getBond(sourceAtom1, sourceAtom2);

                    if (RBond != null) {

                        IAtom targetAtom1 = indexJ;
                        IAtom targetAtom2 = indexJPlus;
                        IBond PBond = PMol.getMolecule().getBond(targetAtom1, targetAtom2);

                        if (PBond != null) {

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
                    }
                }
            }
        }
        boolean flag = false;
        float size = firstSolution.size();
        int sourceAtomCount = 0;
        int targetAtomCount = 0;
        if (!removeHydrogen) {
            sourceAtomCount = RMol.getMolecule().getAtomCount();
            targetAtomCount = PMol.getMolecule().getAtomCount();
        } else {
            sourceAtomCount = RMol.getMolecule().getAtomCount() - getHCount(RMol.getMolecule());
            targetAtomCount = PMol.getMolecule().getAtomCount() - getHCount(PMol.getMolecule());
        }
        if ((size == sourceAtomCount && score / 2 == RMol.getMolecule().getBondCount()) &&
                (targetAtomCount >= size && PMol.getMolecule().getBondCount() >= score / 2)) {
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
            mcs.searchMCS(removeHydrogen);

            firstSolution.clear();
            allMCS.clear();
            allAtomMCS.clear();
            firstAtomMCS.clear();

            firstSolution.putAll(mcs.getFirstMapping());
            allMCS.addAll(mcs.getAllMapping());

            firstAtomMCS.putAll(mcs.getFirstAtomMapping());
            allAtomMCS.addAll(mcs.getAllAtomMapping());


        } catch (IOException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private void singleMapping() {
        try {

            mcs = new SingleMappingHandler();
            mcs.set(RMol, PMol);
            mcs.searchMCS(removeHydrogen);

            firstSolution.clear();
            allMCS.clear();
            allAtomMCS.clear();
            firstAtomMCS.clear();

            firstSolution.putAll(mcs.getFirstMapping());
            allMCS.addAll(mcs.getAllMapping());

            firstAtomMCS.putAll(mcs.getFirstAtomMapping());
            allAtomMCS.addAll(mcs.getAllAtomMapping());


        } catch (IOException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
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
}

