/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.filters;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.smsd.tools.BondEnergies;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;

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
 * You should have received eAtom copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/**
 * @cdk.module smsd
 */
public class ChemicalFilters {

    private List<TreeMap<Integer, Integer>> allMCS = null;
    private TreeMap<Integer, Integer> firstSolution = null;
    private List<Map<IAtom, IAtom>> allAtomMCS = null;
    private Map<IAtom, IAtom> firstAtomMCS = null;
    private List<Integer> StereoScore = null;
    private List<Integer> fragmentSize = null;
    private List<Double> bEnergies = null;
    private MolHandler RMol = null;
    private MolHandler PMol = null;

    /**
     *
     * @param allMCS
     * @param allAtomMCS
     * @param firstSolution
     * @param firstAtomMCS
     * @param RMol
     * @param PMol
     */
    public ChemicalFilters(List<TreeMap<Integer, Integer>> allMCS,
            List<Map<IAtom, IAtom>> allAtomMCS,
            TreeMap<Integer, Integer> firstSolution,
            Map<IAtom, IAtom> firstAtomMCS,
            MolHandler RMol,
            MolHandler PMol) {
        this.allAtomMCS = allAtomMCS;
        this.allMCS = allMCS;
        this.firstAtomMCS = firstAtomMCS;
        this.firstSolution = firstSolution;
        this.PMol = PMol;
        this.RMol = RMol;

        StereoScore = new ArrayList<Integer>();
        fragmentSize = new ArrayList<Integer>();
        bEnergies = new ArrayList<Double>();

    }

    /**
     *
     */
    public synchronized void sortResultsByStereoAndBondMatch() {

//        System.out.println("\nSort By ResultsByStereoAndBondMatch");

        Map<Integer, TreeMap<Integer, Integer>> allStereoMCS = new HashMap<Integer, TreeMap<Integer, Integer>>();
        Map<Integer, Map<IAtom, IAtom>> allStereoAtomMCS = new HashMap<Integer, Map<IAtom, IAtom>>();
        Map<Integer, Double> stereoScoreMap = new HashMap<Integer, Double>();
        Map<Integer, Boolean> stereoSelectionScore = new HashMap<Integer, Boolean>();
        Map<Integer, Integer> fragmentScoreMap = new TreeMap<Integer, Integer>();
        Map<Integer, Double> EnergySelectionMap = new TreeMap<Integer, Double>();

        Integer Index = 0;

        for (Map<IAtom, IAtom> atomsMCS : allAtomMCS) {
            allStereoAtomMCS.put(Index, atomsMCS);
            stereoScoreMap.put(Index, 0.0);
            stereoSelectionScore.put(Index, false);
            Index++;
        }
        Index = 0;

        for (TreeMap<Integer, Integer> MCS : allMCS) {
            allStereoMCS.put(Index, MCS);
//            System.out.println(" Stereo Key " + Index + " Value: " + allStereoMCS.get(Index));
            Index++;
        }

        Index = 0;
        for (Double score : bEnergies) {
            EnergySelectionMap.put(Index, score);
            Index++;
        }
        Index = 0;
        for (Integer score : fragmentSize) {
            fragmentScoreMap.put(Index, score);
            Index++;
        }


        IAtomContainer Reactant = RMol.getMolecule();
        IAtomContainer Product = PMol.getMolecule();
        try {
            CDKHueckelAromaticityDetector.detectAromaticity(Reactant);
            CDKHueckelAromaticityDetector.detectAromaticity(Product);
        } catch (CDKException ex) {
            Logger.getLogger(ChemicalFilters.class.getName()).log(Level.SEVERE, null, ex);
        }
        boolean stereoMatchFlag = false;

        for (Integer Key : allStereoMCS.keySet()) {
            double score = 0.0;
//            System.out.println("\nStart score " + score);
            TreeMap<Integer, Integer> atomsMCS = allStereoMCS.get(Key);
            Map<IAtom, IAtom> atomMapMCS = allStereoAtomMCS.get(Key);

            score = getAtomScore(score, atomMapMCS, Reactant, Product);
            Map<IBond, IBond> bondMaps = makeBondMapsOfAtomMaps(RMol.getMolecule(), PMol.getMolecule(), atomsMCS);
            IAtomContainer subgraphRContainer = getMappedFragment(RMol.getMolecule(), atomMapMCS, 1);
            IAtomContainer subgraphPContainer = getMappedFragment(PMol.getMolecule(), atomMapMCS, 2);

            double score1 = getBondScore(score, bondMaps);

            if ((score1 - score) > 0) {
                stereoMatchFlag = true;
                score = score1;
            }
//            System.out.println("\nStart score1 " + score);
            score = getRingMatchScore(score, subgraphRContainer, subgraphPContainer);

            stereoScoreMap.put(Key, score);
        }

        boolean flag = false;
        if (stereoMatchFlag) {

            Map<Integer, Integer> sortedMap = sortMapByValueInDecendingOrder(stereoScoreMap);

            flag = true;
//            System.out.println(" sortedMap `Size: " + sortedMap.size());

            firstSolution.clear();
            allMCS.clear();
            allAtomMCS.clear();
            firstAtomMCS.clear();
            StereoScore.clear();
            fragmentSize.clear();
            bEnergies.clear();

            int higest_stereoScore = 0;
            for (Integer I : sortedMap.keySet()) {
                higest_stereoScore = sortedMap.get(I);
                break;
            }

            /*Put back the sorted solutions*/

            int counter = 0;
            for (Integer I : sortedMap.keySet()) {


                if (higest_stereoScore == sortedMap.get(I)) {
                    allAtomMCS.add(counter, allStereoAtomMCS.get(I));
                    allMCS.add(counter, allStereoMCS.get(I));
                    StereoScore.add(counter, sortedMap.get(I));
                    fragmentSize.add(counter, fragmentScoreMap.get(I));
                    bEnergies.add(counter, EnergySelectionMap.get(I));
                    counter++;
                }
//                System.out.println("Sorted Map Key " + I + " Sorted Value: " + sortedMap.get(I));
//                System.out.println("allStereoMCS Key " + I + " Sorted Value: " + allStereoMCS.get(I));

            }
            if (flag) {
                firstSolution.putAll(allMCS.get(0));
                firstAtomMCS.putAll(allAtomMCS.get(0));
            }
        }

    }

    private Map<IBond, IBond> makeBondMapsOfAtomMaps(IAtomContainer ac1, IAtomContainer ac2, TreeMap<Integer, Integer> mappings) {

        HashMap<IBond, IBond> maps = new HashMap<IBond, IBond>();

        for (IAtom atoms : ac1.atoms()) {

            int ac1AtomNumber = ac1.getAtomNumber(atoms);

            if (mappings.containsKey(ac1AtomNumber)) {

                int ac2AtomNumber = mappings.get(ac1AtomNumber);

                List<IAtom> connectedAtoms = ac1.getConnectedAtomsList(atoms);

                for (IAtom cAtoms : connectedAtoms) {
                    int ac1ConnectedAtomNumber = ac1.getAtomNumber(cAtoms);

                    if (mappings.containsKey(ac1ConnectedAtomNumber)) {
                        {

                            int ac2ConnectedAtomNumber = mappings.get(ac1ConnectedAtomNumber);

                            IBond ac1Bond = ac1.getBond(atoms, cAtoms);

                            IBond ac2Bond = ac2.getBond(ac2.getAtom(ac2AtomNumber), ac2.getAtom(ac2ConnectedAtomNumber));

                            if (ac2Bond == null) {
                                ac2Bond = ac2.getBond(ac2.getAtom(ac2ConnectedAtomNumber), ac2.getAtom(ac2AtomNumber));

                            }

                            if (ac1Bond != null && ac2Bond != null) {
                                maps.put(ac1Bond, ac2Bond);
                            }

                        }


                    }
                }
            }
        }

//        System.out.println("Mol Map size:" + maps.size());
        return maps;

    }

    private synchronized int getMappedMoleculeFragmentSize(Map<IAtom, IAtom> MCSAtomSolution) {

//      System.out.println("Mol Size Eorg: " + RMol.getMolecule().getAtomCount() + " , Mol Size Porg: " + PMol.getMolecule().getAtomCount());

        IAtomContainer Educt = DefaultChemObjectBuilder.getInstance().newMolecule(RMol.getMolecule());
        IAtomContainer Product = DefaultChemObjectBuilder.getInstance().newMolecule(PMol.getMolecule());


        if (MCSAtomSolution != null) {
            for (Map.Entry<IAtom, IAtom> map : MCSAtomSolution.entrySet()) {

                IAtom atomE = map.getKey();
                IAtom atomP = map.getValue();
                Educt.removeAtomAndConnectedElectronContainers(atomE);
                Product.removeAtomAndConnectedElectronContainers(atomP);

            }
        }
        boolean EductFragmentFlag = true;
        boolean ProductFragmentFlag = true;

        IAtomContainerSet EductFragmentMolSet = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
        IAtomContainerSet ProductFragmentMolSet = DefaultChemObjectBuilder.getInstance().newMoleculeSet();



        int countEFrag = 0;

        if (Educt.getAtomCount() > 0) {
            EductFragmentFlag = ConnectivityChecker.isConnected(Educt);
            if (!EductFragmentFlag) {
                EductFragmentMolSet.add(ConnectivityChecker.partitionIntoMolecules(Educt));

            } else {

                EductFragmentMolSet.addAtomContainer(Educt);
            }

            countEFrag = EductFragmentMolSet.getAtomContainerCount();
        }
        int countPFrag = 0;

        if (Product.getAtomCount() > 0) {
            ProductFragmentFlag = ConnectivityChecker.isConnected(Product);

            if (!ProductFragmentFlag) {
//
                ProductFragmentMolSet.add(ConnectivityChecker.partitionIntoMolecules(Product));

            } else {

                ProductFragmentMolSet.addAtomContainer(Product);
            }
            countPFrag = ProductFragmentMolSet.getAtomContainerCount();
        }

        int FragmentSize = countEFrag + countPFrag;

        return FragmentSize;
    }

    public synchronized void sortResultsByFragments() {

//        System.out.println("\nSort By Fragment");
        Map<Integer, TreeMap<Integer, Integer>> allFragmentMCS = new TreeMap<Integer, TreeMap<Integer, Integer>>();
        Map<Integer, Integer> StereoScoreCopy = new TreeMap<Integer, Integer>();
        Map<Integer, Double> EnergySelectionCopy = new TreeMap<Integer, Double>();

        Map<Integer, Map<IAtom, IAtom>> allFragmentAtomMCS = new TreeMap<Integer, Map<IAtom, IAtom>>();
        Map<Integer, Integer> fragmentScoreMap = new TreeMap<Integer, Integer>();
        Map<Integer, Integer> fragmentSelectionScore = new TreeMap<Integer, Integer>();


        Integer Index = 0;

        for (Map<IAtom, IAtom> atomsMCS : allAtomMCS) {
            allFragmentAtomMCS.put(Index, new HashMap<IAtom, IAtom>(atomsMCS));
            fragmentScoreMap.put(Index, -1);
            fragmentSelectionScore.put(Index, -1);
            Index++;
        }
        Index = 0;

        for (TreeMap<Integer, Integer> MCS : allMCS) {
            allFragmentMCS.put(Index, new TreeMap<Integer, Integer>(MCS));
            Index++;
        }


        Index = 0;
        for (Integer MCS : StereoScore) {
            StereoScoreCopy.put(Index, MCS);
            Index++;
        }

        Index = 0;
        for (Double score : bEnergies) {
            EnergySelectionCopy.put(Index, score);
            Index++;
        }

        Integer _minFragmentScore = 9999;
        for (Integer Key : allFragmentAtomMCS.keySet()) {


            Map<IAtom, IAtom> mcsAtom = allFragmentAtomMCS.get(Key);
            int FragmentCount = getMappedMoleculeFragmentSize(mcsAtom);
//            System.out.println("FragmentCount " + FragmentCount);

            fragmentScoreMap.put(Key, FragmentCount);
            fragmentSelectionScore.put(Key, FragmentCount);
            if (_minFragmentScore > FragmentCount) {
                _minFragmentScore = FragmentCount;
            }
        }

        for (Integer I : fragmentScoreMap.keySet()) {
            if (fragmentScoreMap.get(I) > _minFragmentScore) {
                fragmentSelectionScore.remove(I);
                StereoScoreCopy.remove(I);
                allFragmentAtomMCS.remove(I);
                allFragmentMCS.remove(I);
            }

        }

//        System.out.println("allFragmentAtomMCS: " + allFragmentMCS.size());

        boolean flag = false;

        int counter = 0;

        for (Integer I : fragmentSelectionScore.keySet()) {

            if (_minFragmentScore > 0) {

                if (!flag) {
                    flag = true;

                    firstSolution.clear();
                    allMCS.clear();
                    allAtomMCS.clear();
                    firstAtomMCS.clear();
                    StereoScore.clear();
                    fragmentSize.clear();
                    bEnergies.clear();

                }
                allAtomMCS.add(counter, allFragmentAtomMCS.get(I));
                allMCS.add(counter, allFragmentMCS.get(I));
                fragmentSize.add(counter, fragmentSelectionScore.get(I));
                StereoScore.add(counter, StereoScoreCopy.get(I));
                fragmentSize.add(counter, fragmentScoreMap.get(I));
                bEnergies.add(counter, EnergySelectionCopy.get(I));
                counter++;
//                System.out.println("Fragment Key " + I + " Size: " + fragmentScoreMap.get(I));
//                System.out.println("Fragment MCS " + allFragmentMCS.get(I) + " Fragment Value: " + fragmentSelectionScore.get(I));

            }

        }

        if (flag) {
            firstSolution.putAll(allMCS.get(0));
            firstAtomMCS.putAll(allAtomMCS.get(0));
        }

    }

    private synchronized Double getMappedMoleculeEnergies(TreeMap<Integer, Integer> MCSAtomSolution) throws CDKException {

//        System.out.println("\nSort By Energies");
        BondEnergies bondEnergy = BondEnergies.getInstance();
        double totalBondEnergy = -9999.0;

        IAtomContainer Educt = DefaultChemObjectBuilder.getInstance().newMolecule(RMol.getMolecule());
        IAtomContainer Product = DefaultChemObjectBuilder.getInstance().newMolecule(PMol.getMolecule());


        Iterator<IAtom> eIterator = Educt.atoms().iterator();
        Iterator<IAtom> pIterator = Product.atoms().iterator();

        while (eIterator.hasNext()) {
            eIterator.next().setFlag(0, false);
        }

        while (pIterator.hasNext()) {
            pIterator.next().setFlag(0, false);
        }


        if (MCSAtomSolution != null) {
            for (Map.Entry<Integer, Integer> map : MCSAtomSolution.entrySet()) {
                int ENum = map.getKey();
                int PNum = map.getValue();

                IAtom eAtom = Educt.getAtom(ENum);
                IAtom pAtom = Product.getAtom(PNum);

                eAtom.setFlag(0, true);
                pAtom.setFlag(0, true);
            }
        }

        if (MCSAtomSolution != null) {

            Double eEnergy = 0.0;
            for (int i = 0; i < Educt.getBondCount(); i++) {
                IBond bond = Educt.getBond(i);
                if ((bond.getAtom(0).getFlag(0) == true && bond.getAtom(1).getFlag(0) == false) || (bond.getAtom(0).getFlag(0) == false && bond.getAtom(1).getFlag(0) == true)) {
                    Integer val = bondEnergy.getEnergies(bond.getAtom(0), bond.getAtom(1), bond.getOrder());
                    if (val != null) {
                        eEnergy += val;
                    }
                }

            }
            Double pEnergy = 0.0;
            for (int j = 0; j < Product.getBondCount(); j++) {
                IBond bond = Product.getBond(j);
                if ((bond.getAtom(0).getFlag(0) == true && bond.getAtom(1).getFlag(0) == false) || (bond.getAtom(0).getFlag(0) == false && bond.getAtom(1).getFlag(0) == true)) {
                    Integer val = bondEnergy.getEnergies(bond.getAtom(0), bond.getAtom(1), bond.getOrder());
                    if (val != null) {
                        eEnergy += val;
                    }
                }
            }
            totalBondEnergy = eEnergy + pEnergy;
        }
        return totalBondEnergy;
    }

    public synchronized void sortResultsByEnergies() throws CDKException {

//        System.out.println("\nSort By Energies");
        Map<Integer, TreeMap<Integer, Integer>> allEnergyMCS = new TreeMap<Integer, TreeMap<Integer, Integer>>();
        Map<Integer, Map<IAtom, IAtom>> allAtomEnergyMCS = new TreeMap<Integer, Map<IAtom, IAtom>>();

        Map<Integer, Integer> StereoScoreMap = new TreeMap<Integer, Integer>();
        Map<Integer, Integer> fragmentScoreMap = new TreeMap<Integer, Integer>();
        Map<Integer, Double> EnergySelectionMap = new TreeMap<Integer, Double>();

        Integer Index = 0;

        for (Map<IAtom, IAtom> atomsMCS : allAtomMCS) {
            allAtomEnergyMCS.put(Index, new HashMap<IAtom, IAtom>(atomsMCS));
            EnergySelectionMap.put(Index, -1.0);
            Index++;
        }

        Index = 0;

        for (TreeMap<Integer, Integer> MCS : allMCS) {
            allEnergyMCS.put(Index, new TreeMap<Integer, Integer>(MCS));
            Index++;
        }
        Index = 0;
        for (Integer score : StereoScore) {
            StereoScoreMap.put(Index, score);
            Index++;
        }
        Index = 0;
        for (Integer score : fragmentSize) {
            fragmentScoreMap.put(Index, score);
            Index++;
        }

        for (Integer Key : allEnergyMCS.keySet()) {
            TreeMap<Integer, Integer> mcsAtom = allEnergyMCS.get(Key);
            Double Energies = getMappedMoleculeEnergies(mcsAtom);
            EnergySelectionMap.put(Key, Energies);
        }
        Map<Integer, Double> sortedEnergyMap = sortMapByValueInAccendingOrder(EnergySelectionMap);
        boolean flag = false;

        int counter = 0;
        for (Integer I : sortedEnergyMap.keySet()) {

            if (!flag) {
                flag = true;

                firstSolution.clear();
                allMCS.clear();
                allAtomMCS.clear();
                firstAtomMCS.clear();
                StereoScore.clear();
                fragmentSize.clear();
                bEnergies.clear();

            }

            allAtomMCS.add(counter, allAtomEnergyMCS.get(I));
            allMCS.add(counter, allEnergyMCS.get(I));
            fragmentSize.add(counter, fragmentScoreMap.get(I));
            StereoScore.add(counter, StereoScoreMap.get(I));
            bEnergies.add(counter, sortedEnergyMap.get(I));
            counter++;
//            System.out.println("Energy Key " + I + " Size: " + fragmentScoreMap.get(I));
//            System.out.println("Energy " + allEnergyMCS.get(I) + " Sorted Energy Value: " + sortedEnergyMap.get(I));

        }

        if (flag) {
            firstSolution.putAll(allMCS.get(0));
            firstAtomMCS.putAll(allAtomMCS.get(0));
        }

    }

    static Map<Integer, Double> sortMapByValueInAccendingOrder(Map<Integer, Double> map) {
        List<Map.Entry<Integer, Double>> list = new LinkedList<Map.Entry<Integer, Double>>(map.entrySet());
        // Sort the list using an annonymous inner class implementing Comparator for the compare method
        java.util.Collections.sort(list, new Comparator<Map.Entry<Integer, Double>>() {

            @Override
            public int compare(Map.Entry<Integer, Double> entry, Map.Entry<Integer, Double> entry1) {
                // Return 0 for eAtom match, -1 for less than and +1 for more then (Aceending Order Sort)
                return (entry.getValue().equals(entry1.getValue()) ? 0 : (entry.getValue() > entry1.getValue() ? 1 : -1));
            }
        });
        // logger.info(list);
        Map<Integer, Double> result = new LinkedHashMap<Integer, Double>();
        for (Iterator<Map.Entry<Integer, Double>> it = list.iterator(); it.hasNext();) {
            Map.Entry<Integer, Double> entry = it.next();
            result.put(entry.getKey(), entry.getValue());
        }
        return result;
    }

    static Map<Integer, Integer> sortMapByValueInDecendingOrder(Map<Integer, Double> map) {
        List<Map.Entry<Integer, Double>> list = new LinkedList<Map.Entry<Integer, Double>>(map.entrySet());
        // Sort the list using an annonymous inner class implementing Comparator for the compare method
        java.util.Collections.sort(list, new Comparator<Map.Entry<Integer, Double>>() {

            @Override
            public int compare(Map.Entry<Integer, Double> entry, Map.Entry<Integer, Double> entry1) {
                // Return 0 for eAtom match, -1 for less than and +1 for more then (Decending Order Sort)
                return (entry.getValue().equals(entry1.getValue()) ? 0 : (entry.getValue() < entry1.getValue() ? 1 : -1));
            }
        });
        // logger.info(list);
        Map<Integer, Integer> result = new LinkedHashMap<Integer, Integer>();
        for (Iterator<Map.Entry<Integer, Double>> it = list.iterator(); it.hasNext();) {
            Map.Entry<Integer, Double> entry = it.next();
            result.put(entry.getKey(), entry.getValue().intValue());
        }
        return result;
    }

    /**
     * 
     * @return
     */
    public List<Double> getSortedEnergy() {
        return bEnergies;
    }

    /**
     * 
     * @return
     */
    public List<Integer> getSortedFragment() {
        return fragmentSize;
    }

    /**
     * 
     * @return
     */
    public List<Integer> getStereoMatches() {

        return StereoScore;
    }

    private IAtomContainer getMappedFragment(IAtomContainer molecule, Map<IAtom, IAtom> atomsMCS, int i) {
        IAtomContainer subgraphContainer = DefaultChemObjectBuilder.getInstance().newAtomContainer(molecule);
        atomsMCS.keySet();
        if (i == 1) {
            for (IAtom atoms : molecule.atoms()) {
                if (!atomsMCS.containsKey(atoms)) {
                    subgraphContainer.removeAtomAndConnectedElectronContainers(atoms);
                }
            }
        } else if (i == 2) {
            for (IAtom atoms : molecule.atoms()) {
                if (!atomsMCS.containsValue(atoms)) {
                    subgraphContainer.removeAtomAndConnectedElectronContainers(atoms);
                }
            }
        } else {

            System.out.println("1: Reactant, 2: Product " + i + "is invalid option");
        }

        return subgraphContainer;

    }

    private double getAtomScore(double score, Map<IAtom, IAtom> atomMapMCS, IAtomContainer Reactant, IAtomContainer Product) {
        for (Map.Entry<IAtom, IAtom> mappings : atomMapMCS.entrySet()) {
            IAtom rAtom = mappings.getKey();
            IAtom pAtom = mappings.getValue();

            int rHCount = 0;
            int pHCount = 0;
            double rBO = Reactant.getBondOrderSum(rAtom);
            double pBO = Product.getBondOrderSum(pAtom);

            if (rAtom.getHydrogenCount() != null) {
                rHCount = rAtom.getHydrogenCount();
            }
            if (pAtom.getHydrogenCount() != null) {
                pHCount = pAtom.getHydrogenCount();
            }

            int HScore = Math.abs(rHCount - pHCount);
            double BOScore = Math.abs(rBO - pBO);

            if (rHCount != pHCount) {
                score = score - HScore;
            } else {
                score = score + HScore;
            }

            if (rBO != pBO) {
                score = score - BOScore;
            } else {
                score = score + BOScore;
            }
        }

        return score;

    }

    private double getBondScore(double score, Map<IBond, IBond> bondMaps) {
        for (Map.Entry<IBond, IBond> matchedBonds : bondMaps.entrySet()) {

            IBond RBond = matchedBonds.getKey();
            IBond PBond = matchedBonds.getValue();
            int RBondType = RBond.getOrder().ordinal();
            int PBondType = PBond.getOrder().ordinal();


            if (RBond.getFlag(CDKConstants.ISAROMATIC) == PBond.getFlag(CDKConstants.ISAROMATIC) && RBondType == PBondType) {
                score += 2;
            } else if (RBond.getFlag(CDKConstants.ISAROMATIC) && PBond.getFlag(CDKConstants.ISAROMATIC)) {
                score += 4;
            }
            if (RBond.getStereo() != PBond.getStereo()) {
                score -= 2;
            } else {
                score += 2;
            }
            if (RBondType != PBondType) {
                score = score - Math.abs(RBondType - PBondType);
            } else {
                score = score + Math.abs(RBondType - PBondType);
            }
            if (RBond.getAtom(0).getFormalCharge() == PBond.getAtom(0).getFormalCharge()) {
                score = score + Math.abs(RBond.getAtom(0).getFormalCharge() - PBond.getAtom(0).getFormalCharge());

            } else {
                score = score - Math.abs(RBond.getAtom(0).getFormalCharge() - PBond.getAtom(0).getFormalCharge());

            }
            if (RBond.getAtom(1).getFormalCharge() == PBond.getAtom(1).getFormalCharge()) {
                score = score + Math.abs(RBond.getAtom(1).getFormalCharge() - PBond.getAtom(1).getFormalCharge());

            } else {
                score = score - Math.abs(RBond.getAtom(1).getFormalCharge() - PBond.getAtom(1).getFormalCharge());
            }
        }
        return score;
    }

    private double getRingMatchScore(double score, IAtomContainer subgraphRContainer, IAtomContainer subgraphPContainer) {

        SSSRFinder ringFinderR = new SSSRFinder(subgraphRContainer);
        IRingSet rRings = ringFinderR.findRelevantRings();
        SSSRFinder ringFinderP = new SSSRFinder(subgraphPContainer);
        IRingSet pRings = ringFinderP.findRelevantRings();

        int rLength = RingSetManipulator.getAtomCount(rRings);
        int pLength = RingSetManipulator.getAtomCount(pRings);

        for (IAtomContainer ac : RingSetManipulator.getAllAtomContainers(rRings)) {
            boolean flag = true;
            for (IAtom a : ac.atoms()) {
                for (Map<IAtom, IAtom> aMCS : allAtomMCS) {

                    if (!aMCS.containsKey(a)) {
                        flag = false;
                        break;
                    }
                }
            }

            if (flag) {
                score += 10;
            }
        }

        for (IAtomContainer ac : RingSetManipulator.getAllAtomContainers(pRings)) {
            boolean flag = true;
            for (IAtom a : ac.atoms()) {
                for (Map<IAtom, IAtom> aMCS : allAtomMCS) {

                    if (!aMCS.containsValue(a)) {
                        flag = false;
                        break;
                    }
                }
            }

            if (flag) {
                score += 10;
            }
        }

        if (rLength > 0) {

            if (rLength == pLength) {
                score += rLength * 2;
            }

            if (rLength > pLength) {
                score += (rLength - pLength) * 2;
            }

        }
        return score;
    }
}
