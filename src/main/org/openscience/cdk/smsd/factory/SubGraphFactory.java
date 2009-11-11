/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.factory;

import java.io.IOException;
import java.math.BigDecimal;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.smsd.algorithm.single.SingleMappingHandler;
import org.openscience.cdk.smsd.algorithm.vflib.VFlibTurboHandler;
import org.openscience.cdk.smsd.core.tools.EBIException;
import org.openscience.cdk.smsd.filters.ChemicalFilters;
import org.openscience.cdk.smsd.global.BondType;
import org.openscience.cdk.smsd.global.TimeOut;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.IMCSAlgorithm;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;

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
/**
 * @cdk.module smsd
 */
public class SubGraphFactory implements IMCSAlgorithm {

    private List<TreeMap<Integer, Integer>> allMCS = null;
    private TreeMap<Integer, Integer> firstSolution = null;
    private List<Map<IAtom, IAtom>> allAtomMCS = null;
    private Map<IAtom, IAtom> firstAtomMCS = null;
    private double tanimoto = -1;
    private double euclidean = -1;
    private MolHandler RMol = null;
    private MolHandler PMol = null;
    private List<Integer> StereoScore;
    private List<Integer> fragmentSize;
    private List<Double> bEnergies;
    private boolean removeHydrogen = false;
    private boolean stereoFilter = false;
    private boolean fragmentFilter = false;
    private boolean energyFilter = false;
    private boolean subGraphFlag = false;

    /**
     * 
     * @param BondTypeFlag 
     * @param removeHydrogen
     * @param stereoFilter
     * @param fragmentFilter
     * @param energyFilter
     * @throws java.lang.Exception
     */
    public SubGraphFactory(boolean BondTypeFlag, boolean removeHydrogen, boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) throws Exception {

        this.removeHydrogen = removeHydrogen;
        this.stereoFilter = stereoFilter;
        this.fragmentFilter = fragmentFilter;
        this.energyFilter = energyFilter;

        firstSolution = new TreeMap<Integer, Integer>();
        allMCS = new Vector<TreeMap<Integer, Integer>>();
        allAtomMCS = new Vector<Map<IAtom, IAtom>>();
        firstAtomMCS = new HashMap<IAtom, IAtom>();

        TimeOut tmo = TimeOut.getInstance();
        tmo.setTimeOut(0.30);

        BondType BT = BondType.getInstance();
        BT.setBondSensitiveFlag(BondTypeFlag);


    }

    private synchronized void MCSBuilder() {

        try {


            int rBondCount = RMol.getMolecule().getBondCount();
            int pBondCount = PMol.getMolecule().getBondCount();


            /*This is importsnt because CDK fails to Generate makeAtomsMapOfBondsMap if bonds are less than 2*/

            if (rBondCount > 1 && pBondCount > 1) {
                //System.out.println("\nVF-MCS\n");
                VFLibMCS();

            } else {
                SingleMapping();
            }

            System.gc();

        } catch (Exception ex) {
            Logger.getLogger(SubGraphFactory.class.getName()).log(Level.SEVERE, null, ex);
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


            MCSBuilder();
            setChemFilters();
        } catch (EBIException ex) {
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
        //System.out.println(" Container Size " + Reactant.getAtomCount());
        //System.out.println(" Container Size " + Product.getAtomCount());
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

        //System.out.println(" Container Size " + Reactant.getBondCount());
        //System.out.println(" Container Size " + Product.getBondCount());
        this.RMol = new MolHandler(Reactant, false, removeHydrogen);
        this.PMol = new MolHandler(Product, false, removeHydrogen);

        init(RMol, PMol);

    }

    public synchronized void setChemFilters() throws EBIException {
        if (firstAtomMCS != null) {
            ChemicalFilters CF = new ChemicalFilters(allMCS, allAtomMCS, firstSolution, firstAtomMCS, RMol, PMol);

            if (stereoFilter) {
                CF.sortResultsByStereoAndBondMatch();
            }
            if (fragmentFilter) {
                CF.sortResultsByFragments();
            }

            if (energyFilter) {
                CF.sortResultsByEnergies();
            }

            this.StereoScore = CF.getStereoMatches();
            this.fragmentSize = CF.getSortedFragment();
            this.bEnergies = CF.getSortedEnergy();
        }
    }

    @Override
    public synchronized Integer getFragmentSize(int Key) {
        Integer Value = null;
        // System.out.println("Key" + Key);
        if (firstSolution.size() > 0 && fragmentSize.size() > Key && Key >= 0) {
            Value = fragmentSize.get(Key);
        }

        return Value;

    }

    @Override
    public Integer getStereoScore(int Key) {

        Integer Value = null;
//        System.out.println(StereoScore.size() + " :Key " + Key);
        if (StereoScore.size() > Key && Key >= 0) {
            Value = StereoScore.get(Key);
//             System.out.println("\nScore " + Value);
        }
        return Value;
    }

    @Override
    public Double getEnergyScore(int Key) {

        Double Value = null;
//        System.out.println(StereoScore.size() + " :Key " + Key);
        if (bEnergies.size() > Key && Key >= 0) {
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
            //System.out.println("Total Sol= " + allMCS.size());
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
            //System.out.println("firstSolution: " + firstAtomMCS);
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


        int a = 0;
        int b = 0;
        if (!removeHydrogen) {
            a = RMol.getMolecule().getAtomCount();
            b = PMol.getMolecule().getAtomCount();
        } else {
            a = RMol.getMolecule().getAtomCount() - HCount(RMol.getMolecule());
            b = PMol.getMolecule().getAtomCount() - HCount(PMol.getMolecule());
        }
        double c = getFirstMapping().size();

        tanimoto = (c) / (a + b - c);

//        System.out.println("c" + c);
//        System.out.println("a" + a);
//        System.out.println("b" + b);
//        System.out.println("a+b-c" + (a + b - c));
//


        BigDecimal tan = new BigDecimal(tanimoto);

        tan = tan.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
        //System.out.println("Tanimoto Coefficient: " + (c)/(a+b-c));
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
                if (indexI.equals(indexIPlus) && indexJ.equals(indexJPlus)) {

                    IAtom R1 = indexI;
                    IAtom R2 = indexIPlus;

                    IBond RBond = Reactant.getBond(R1, R2);

                    if (RBond != null) {

                        IAtom P1 = indexJ;
                        IAtom P2 = indexJPlus;
                        IBond PBond = Product.getBond(P1, P2);

                        if (PBond != null) {

                            if (RBond.getStereo() != PBond.getStereo()) {
                                Score++;
                            }
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
     * @return true if ac1 is a subgraph of ac2
     */
    @Override
    public boolean isSubgraph() {
        if (firstAtomMCS == null || firstAtomMCS.isEmpty()) {

            return false;
        }
        BondType BT = BondType.getInstance();
        int score = 0;
        for (Map.Entry<IAtom, IAtom> mappingI : firstAtomMCS.entrySet()) {
            IAtom indexI = mappingI.getKey();
            IAtom indexJ = mappingI.getValue();
            for (Map.Entry<IAtom, IAtom> mappingJ : firstAtomMCS.entrySet()) {

                IAtom indexIPlus = mappingJ.getKey();
                IAtom indexJPlus = mappingJ.getValue();
                if (!indexI.equals(indexIPlus) && !indexJ.equals(indexJPlus)) {

                    IAtom R1 = indexI;
                    IAtom R2 = indexIPlus;

                    IBond RBond = RMol.getMolecule().getBond(R1, R2);

                    if (RBond != null) {

                        IAtom P1 = indexJ;
                        IAtom P2 = indexJPlus;
                        IBond PBond = PMol.getMolecule().getBond(P1, P2);

                        if (PBond != null) {

                            if (BT.getBondSensitiveFlag()) {
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
        int a = 0;
        int b = 0;
        if (!removeHydrogen) {
            a = RMol.getMolecule().getAtomCount();
            b = PMol.getMolecule().getAtomCount();
        } else {
            a = RMol.getMolecule().getAtomCount() - HCount(RMol.getMolecule());
            b = PMol.getMolecule().getAtomCount() - HCount(PMol.getMolecule());
        }
        if (size == a && score / 2 == RMol.getMolecule().getBondCount()) {

            if (b >= size && PMol.getMolecule().getBondCount() >= score / 2) {
                flag = true;
            }

        }
        return flag;
    }

    @Override
    public double getEuclideanDistance() throws IOException {
        int decimalPlaces = 4;
        double a = 0;
        double b = 0;
        if (!removeHydrogen) {
            a = RMol.getMolecule().getAtomCount();
            b = PMol.getMolecule().getAtomCount();
        } else {
            a = RMol.getMolecule().getAtomCount() - HCount(RMol.getMolecule());
            b = PMol.getMolecule().getAtomCount() - HCount(PMol.getMolecule());
        }
        double c = getFirstMapping().size();

        euclidean = Math.sqrt(a + b - 2 * c);

        BigDecimal dist = new BigDecimal(euclidean);

        dist = dist.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
        //System.out.println("euclidean Coefficient: " + (c)/(a+b-c));
        euclidean = dist.doubleValue();
        return euclidean;
    }

    private void VFLibMCS() {

        try {
            VFlibTurboHandler mcs = new VFlibTurboHandler();
            mcs.set(RMol, PMol, removeHydrogen);
            this.subGraphFlag = mcs.isSubgraph();

            firstSolution.clear();
            allMCS.clear();
            allAtomMCS.clear();
            firstAtomMCS.clear();

            if (subGraphFlag) {
                firstSolution.putAll(mcs.getFirstMapping());
                allMCS.addAll(mcs.getAllMapping());


                firstAtomMCS.putAll(mcs.getFirstAtomMapping());
                allAtomMCS.addAll(mcs.getAllAtomMapping());
            }
//
//
//            System.out.println("First Atom MCS: " + firstSolution);
//            System.out.println("First Atom AllMCS: " + allMCS);


        } catch (IOException ex) {
            Logger.getLogger(SubGraphFactory.class.getName()).log(Level.SEVERE, null, ex);
        } catch (EBIException ex) {
            Logger.getLogger(SubGraphFactory.class.getName()).log(Level.SEVERE, null, ex);
        }


    }

    private void SingleMapping() {
        try {

            SingleMappingHandler mcs = new SingleMappingHandler();

            mcs.set(RMol, PMol);

            mcs.search_MCS(removeHydrogen);

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
        } catch (EBIException ex) {
            Logger.getLogger(MCSFactory.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private int HCount(IAtomContainer molecule) {
        int count = 0;
        for (IAtom atom : molecule.atoms()) {
            if (atom.getSymbol().equalsIgnoreCase("H")) {
                ++count;
            }
        }

        return count;
    }
}
