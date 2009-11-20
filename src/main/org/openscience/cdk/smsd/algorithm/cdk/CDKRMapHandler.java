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
 * You should have received sourceAtom copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.smsd.algorithm.cdk;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Stack;
import java.util.TreeMap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.smsd.helper.FinalMappings;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.manipulator.BondManipulator;

/**
 * This algorithm derives from the algorithm described in
 * [Tonnelier, C. and Jauffret, Ph. and Hanser, Th. and Jauffret, Ph. and Kaufmann, G.,
 * Machine Learning of generic reactions:
 * 3. An efficient algorithm for maximal common substructure determination,
 * Tetrahedron Comput. Methodol., 1990, 3:351-358] and modified in the thesis of
 * T. Hanser [Unknown BibTeXML type: HAN93].
 *
 * @cdk.module smsd
 */
public class CDKRMapHandler {

//    boolean RonPFlag = false;
    private List<TreeMap<Integer, Integer>> _mapping;
    static IAtomContainer source;
    static IAtomContainer target;
    private boolean timeoutFlag = false;

    /**
     * This function calculates all the possible combinations of MCS
     * @param Molecule1
     * @param Molecule2
     * @throws CDKException
     */
    public void calculateOverlapsAndReduce(IAtomContainer Molecule1, IAtomContainer Molecule2) throws CDKException {

        source = Molecule1;
        target = Molecule2;

        _mapping = new ArrayList<TreeMap<Integer, Integer>>();

        //System.out.println("Searching: ");
        //List overlaps = UniversalIsomorphismTesterBondTypeInSensitive.getSubgraphAtomsMap(source, target);



        if ((source.getAtomCount() == 1) || (target.getAtomCount() == 1)) {

            // System.out.println("Searched Single Molecule Mapping Case: ");

            List overlaps = CDKMCS.checkSingleAtomCases(source, target);
            //timeoutFlag=graphContainer.getTimeOutFlag();


            int nAtomsMatched = overlaps.size();
            nAtomsMatched = (nAtomsMatched > 0) ? 1 : 0;
            if (nAtomsMatched > 0) {
                /*UnComment this to get one Unique Mapping*/
                //List reducedList = removeRedundantMappingsForSingleAtomCase(overlaps);
                //int counter = 0;
                identifyMatchedPartsforSingleAtoms(overlaps, source, target);

            }

        } else {

            /*Time out set by Asad equalt 2 min*/
            //UniversalIsomorphismTesterBondTypeInSensitive.timeout=12000000;

            //source (Mol1), target (Mol2), RGraph1, RGraph2, true, true( search all MCS)

            List overlaps = CDKMCS.search(source, target, new BitSet(), new BitSet(), true, true);


            //timeoutFlag=graphContainer.getTimeOutFlag();
            //System.out.println("Searched: ");

            List reducedList = removeSubGraph(overlaps);


            Stack<List> allMaxOverlaps = getAllMaximum(reducedList);

            //Stack<List> allMaxOverlaps = getAllMaximum(overlaps);


            while (!allMaxOverlaps.empty()) {
//                System.out.println("source: " + source.getAtomCount() + ", target: " + target.getAtomCount() + ", overl: " + allMaxOverlaps.peek().size());
                @SuppressWarnings("unchecked")
                List maxOverlapsAtoms = makeAtomsMapOfBondsMap(allMaxOverlaps.peek(), source, target);
//                System.out.println("size of maxOverlaps: " + maxOverlapsAtoms.size());
                identifyMatchedParts(maxOverlapsAtoms, source, target);
//                identifyMatchedParts(allMaxOverlaps.peek(), source, target);
                allMaxOverlaps.pop();
            }
        }


        FinalMappings.getInstance().set(_mapping);
//
//        System.out.println("Mapping count: " + _mapping.size());
//        System.out.println("Mapping Size: " + _mapping.firstElement());


    }

    /**
     * This function calculates only one solution (exact) because we are looking at the
     * molecules which are exactly same in terms of the bonds and atoms determined by the
     * Fingerprint
     * @param Molecule1
     * @param Molecule2
     * @throws CDKException
     */
    public void calculateOverlapsAndReduceExactMatch(IAtomContainer Molecule1, IAtomContainer Molecule2) throws CDKException {

        source = Molecule1;
        target = Molecule2;

        _mapping = new ArrayList<TreeMap<Integer, Integer>>();

        //System.out.println("Searching: ");
        //List overlaps = UniversalIsomorphismTesterBondTypeInSensitive.getSubgraphAtomsMap(source, target);



        if ((source.getAtomCount() == 1) || (target.getAtomCount() == 1)) {

            // System.out.println("Searched Single Molecule Mapping Case: ");

            List overlaps = CDKMCS.checkSingleAtomCases(source, target);
            //timeoutFlag=graphContainer.getTimeOutFlag();


            int nAtomsMatched = overlaps.size();
            nAtomsMatched = (nAtomsMatched > 0) ? 1 : 0;
            if (nAtomsMatched > 0) {
                /*UnComment this to get one Unique Mapping*/
                //List reducedList = removeRedundantMappingsForSingleAtomCase(overlaps);
                //int counter = 0;
                identifyMatchedPartsforSingleAtoms(overlaps, source, target);

            }

        } else {

            /*Time out set by Asad equalt 2 min*/
            //UniversalIsomorphismTesterBondTypeInSensitive.timeout=12000000;

            //source (Mol1), target (Mol2), RGraph1, RGraph2, true, true( search all MCS)

            List overlaps = CDKMCS.search(source, target, new BitSet(), new BitSet(), true, true);
//            List overlaps = graphContainer.search(source, target, new BitSet(), new BitSet(), true, false);
            //timeoutFlag=graphContainer.getTimeOutFlag();
            //System.out.println("Searched: ");

            List reducedList = removeSubGraph(overlaps);


            Stack<List> allMaxOverlaps = getAllMaximum(reducedList);

            //Stack<List> allMaxOverlaps = getAllMaximum(overlaps);


            while (!allMaxOverlaps.empty()) {
//                System.out.println("source: " + source.getAtomCount() + ", target: " + target.getAtomCount() + ", overl: " + allMaxOverlaps.peek().size());
                @SuppressWarnings("unchecked")
                List maxOverlapsAtoms = makeAtomsMapOfBondsMap(allMaxOverlaps.peek(), source, target);
//                System.out.println("size of maxOverlaps: " + maxOverlapsAtoms.size());
                identifyMatchedParts(maxOverlapsAtoms, source, target);
//                identifyMatchedParts(allMaxOverlaps.peek(), source, target);

                allMaxOverlaps.pop();
            }
        }
        FinalMappings.getInstance().set(_mapping);
//
//        System.out.println("Mapping count: " + _mapping.size());
//        System.out.println("Mapping Size: " + _mapping.firstElement());

    }

    /**
     *
     * @param overlaps
     * @return
     */
    protected List removeSubGraph(List overlaps) {
        @SuppressWarnings("unchecked")
        List reducedList = new ArrayList(overlaps);

        for (int i = 0; i < overlaps.size(); i++) {
            //System.out.println("i: " + i + ", overlaps.size(): " + overlaps.size());
            //System.out.println("overlaps i: " + overlaps.get(i).getClass().getName());
            List gi = (List) overlaps.get(i);


            for (int j = i + 1; j < overlaps.size(); j++) {
                //System.out.println("j: " + j + ", overlaps.size(): " + overlaps.size());
                //System.out.println("overlaps j: " + overlaps.get(j).getClass().getName());

                List gj = (List) overlaps.get(j);

                // Gi included in Gj or Gj included in Gi then
                // reduce the irrelevant solution
                if (gi.size() != gj.size()) {
                    if (isSubgraph(gj, gi)) {
                        reducedList.remove(gi);
                    } else if (isSubgraph(gi, gj)) {
                        reducedList.remove(gj);
                    }
                }

            }
        }
        return reducedList;
    }

    /**
     *
     * @param overlaps
     * @return
     */
    @SuppressWarnings("unchecked")
    protected List removeRedundantMappingsForSingleAtomCase(List overlaps) {
        @SuppressWarnings("unchecked")

        List reducedList = new ArrayList();

        reducedList.add(overlaps.get(0));
        //reducedList.add(overlaps.get(1));

        return reducedList;

    }

    /**
     *  This makes sourceAtom map of matching atoms out of sourceAtom map of matching bonds as produced by the get(Subgraph|Ismorphism)Map methods.
     *
     * @param  l   The list produced by the getMap method.
     * @param  g1  first molecule. Must not be an IQueryAtomContainer.
     * @param  g2  second molecule. May be an IQueryAtomContainer.
     * @return     The mapping found projected on g1. This is sourceAtom List of CDKRMap objects containing Ids of matching atoms.
     */
    public static List<CDKRMap> makeAtomsMapOfBondsMap(List<CDKRMap> l, IAtomContainer g1, IAtomContainer g2) {
        if (l == null) {
            return (l);
        }
        List<CDKRMap> result = new ArrayList<CDKRMap>();
        for (int i = 0; i < l.size(); i++) {
            IBond bond1 = g1.getBond(l.get(i).getId1());
            IBond bond2 = g2.getBond(l.get(i).getId2());
            IAtom[] atom1 = BondManipulator.getAtomArray(bond1);
            IAtom[] atom2 = BondManipulator.getAtomArray(bond2);
            for (int j = 0; j < 2; j++) {
                List<IBond> bondsConnectedToAtom1j = g1.getConnectedBondsList(atom1[j]);
                for (int k = 0; k < bondsConnectedToAtom1j.size(); k++) {
                    if (bondsConnectedToAtom1j.get(k) != bond1) {
                        IBond testBond = bondsConnectedToAtom1j.get(k);
                        for (int m = 0; m < l.size(); m++) {
                            IBond testBond2;
                            if ((l.get(m)).getId1() == g1.getBondNumber(testBond)) {
                                testBond2 = g2.getBond((l.get(m)).getId2());
                                for (int n = 0; n < 2; n++) {
                                    List<IBond> bondsToTest = g2.getConnectedBondsList(atom2[n]);
                                    if (bondsToTest.contains(testBond2)) {
                                        CDKRMap map;
                                        if (j == n) {
                                            map = new CDKRMap(g1.getAtomNumber(atom1[0]), g2.getAtomNumber(atom2[0]));
                                        } else {
                                            map = new CDKRMap(g1.getAtomNumber(atom1[1]), g2.getAtomNumber(atom2[0]));
                                        }
                                        if (!result.contains(map)) {
                                            result.add(map);
                                        }
                                        CDKRMap map2;
                                        if (j == n) {
                                            map2 = new CDKRMap(g1.getAtomNumber(atom1[1]), g2.getAtomNumber(atom2[1]));
                                        } else {
                                            map2 = new CDKRMap(g1.getAtomNumber(atom1[0]), g2.getAtomNumber(atom2[1]));
                                        }
                                        if (!result.contains(map2)) {
                                            result.add(map2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return result;
    }

//    /**
//     *  This makes sourceAtom map of matching atoms out of sourceAtom map of matching bonds as produced by the get(Subgraph|Ismorphism)Map methods.
//     *
//     * @param  l   The list produced by the getMap method.
//     * @param  g1  first molecule. Must not be an IQueryAtomContainer.
//     * @param  g2  second molecule. May be an IQueryAtomContainer.
//     * @return     The mapping found projected on g1. This is sourceAtom List of CDKRMap objects containing Ids of matching atoms.
//     */
//    @SuppressWarnings("unchecked")
//    protected List makeAtomsMapOfBondsMap(List l, IAtomContainer g1, IAtomContainer g2) {
//        if (l == null) {
//            return (l);
//        }
//        List result = new ArrayList();
//        for (int i = 0; i < l.size(); i++) {
//            IBond bond1 = g1.getBond(((CDKRMap) l.get(i)).getId1());
//            IBond bond2 = g2.getBond(((CDKRMap) l.get(i)).getId2());
//            IAtom[] atom1 = BondManipulator.getAtomArray(bond1);
//            IAtom[] atom2 = BondManipulator.getAtomArray(bond2);
//            for (int j = 0; j < 2; j++) {
//                List bondsConnectedToAtom1j = g1.getConnectedBondsList(atom1[j]);
//                for (int k = 0; k < bondsConnectedToAtom1j.size(); k++) {
//                    if (bondsConnectedToAtom1j.get(k) != bond1) {
//                        IBond testBond = (IBond) bondsConnectedToAtom1j.get(k);
//                        for (int m = 0; m < l.size(); m++) {
//                            IBond testBond2;
//                            if (((CDKRMap) l.get(m)).getId1() == g1.getBondNumber(testBond)) {
//                                testBond2 = g2.getBond(((CDKRMap) l.get(m)).getId2());
//                                for (int n = 0; n < 2; n++) {
//                                    List bondsToTest = g2.getConnectedBondsList(atom2[n]);
//                                    if (bondsToTest.contains(testBond2)) {
//                                        CDKRMap map;
//                                        if (j == n) {
//                                            map = new CDKRMap(g1.getAtomNumber(atom1[0]), g2.getAtomNumber(atom2[0]));
//                                        } else {
//                                            map = new CDKRMap(g1.getAtomNumber(atom1[1]), g2.getAtomNumber(atom2[0]));
//                                        }
//                                        if (!result.contains(map)) {
//                                            //System.out.println("Added Solution");
//                                            result.add(map);
//                                        }
//                                        CDKRMap map2;
//                                        if (j == n) {
//                                            map2 = new CDKRMap(g1.getAtomNumber(atom1[1]), g2.getAtomNumber(atom2[1]));
//                                        } else {
//                                            map2 = new CDKRMap(g1.getAtomNumber(atom1[0]), g2.getAtomNumber(atom2[1]));
//                                        }
//                                        if (!result.contains(map2)) {
//                                            result.add(map2);
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        return (result);
//    }
    /**
     *
     * @param overlaps
     * @return
     */
    protected List getMaximum(List overlaps) {


        ArrayList list = null;
        int count = 0;
        for (Object o : overlaps) {
            ArrayList a = (ArrayList) o;
            if (a.size() > count) {
                list = a;
                count = a.size();
            }

        }
        return list;
    }

    /**
     *
     * @param overlaps
     * @return
     */
    protected Stack<List> getAllMaximum(List overlaps) {

        Stack<List> allMaximumMappings = null;

        int count = -1;

        for (Object o : overlaps) {

            ArrayList a = (ArrayList) o;

            //System.out.println("O size" + sourceAtom.size());

            if (a.size() > count) {
                @SuppressWarnings("unchecked")
                List list = new ArrayList(a);
                count = a.size();

                //System.out.println("List size" + list.size());

                //Collection threadSafeList = Collections.synchronizedCollection( list );
                allMaximumMappings = new Stack<List>();
                //allMaximumMappings.clear();
                allMaximumMappings.push(list);
            } else if (a.size() == count) {
                @SuppressWarnings("unchecked")
                List list = new ArrayList(a);
                count = a.size();
                allMaximumMappings.push(list);
            }

        }
        return allMaximumMappings;
    }

    /**
     *
     * @param list
     * @param source
     * @param target
     */
    protected void identifyMatchedParts(List list, IAtomContainer source, IAtomContainer target) {

        List<IAtom> array1 = new ArrayList<IAtom>();
        List<IAtom> array2 = new ArrayList<IAtom>();

        TreeMap<Integer, Integer> atomNumbersFromContainer = new TreeMap<Integer, Integer>();


        /*We have serial numbers of the bonds/Atoms to delete
         * Now we will collect the actual bond/Atoms rather than
         * serial number for deletion. RonP flag check whether Reactant is
         * mapped on Product or Vise Versa
         */



        for (Object o : list) {
            //System.err.print("Map " + o.getClass());
            CDKRMap rmap = (CDKRMap) o;
            IAtom sourceAtom = source.getAtom(rmap.getId1());
            IAtom targetAtom = target.getAtom(rmap.getId2());


            array1.add(sourceAtom);
            array2.add(targetAtom);


//            System.err.print("Map " + sourceAtom.getSymbol());
//            System.err.print("Map " + targetAtom.getSymbol());

            int IndexI = source.getAtomNumber(sourceAtom);
            int IndexJ = target.getAtomNumber(targetAtom);

            atomNumbersFromContainer.put(IndexI, IndexJ);



            /*Added the Mapping Numbers to the FinalMapping*
             */
            _mapping.add(atomNumbersFromContainer);



        }

    }

    /**
     *
     * @param list
     * @param source
     * @param target
     */
    protected void identifyMatchedPartsforSingleAtoms(List list,
            IAtomContainer source,
            IAtomContainer target) {

        List<IAtom> array1 = new ArrayList<IAtom>();
        List<IAtom> array2 = new ArrayList<IAtom>();



        /*We have serial numbers of the bonds/Atoms to delete
         * Now we will collect the actual bond/Atoms rather than
         * serial number for deletion. RonP flag check whether Reactant is
         * mapped on Product or Vise Versa
         */

        TreeMap<Integer, Integer> atomNumbersFromContainer = new TreeMap<Integer, Integer>();



        for (Object o : list) {
            //System.err.print("Map " + o.getClass());
            CDKRMap rmap = (CDKRMap) o;

            IAtom a = source.getAtom(rmap.getId1());
            IAtom b = target.getAtom(rmap.getId2());


            array1.add(a);
            array2.add(b);


//            System.err.print("Map " + sourceAtom.getSymbol());
//            System.err.print("Map " + targetAtom.getSymbol());

            int IndexI = source.getAtomNumber(a);
            int IndexJ = target.getAtomNumber(b);


            atomNumbersFromContainer.put(IndexI, IndexJ);



            /*Added the Mapping Numbers to the FinalMapping*
             */
            _mapping.add(atomNumbersFromContainer);


        }



        //Asad: Uncomment this part of the code to remove the matched parts from the Molecules.

        /*for (IAtom atom : array1) {
        source.removeAtomAndConnectedElectronContainers(atom);
        }
        for (IAtom atom : array2) {
        target.removeAtomAndConnectedElectronContainers(atom);
        }*/
    }

    /**
     *
     * @param rmaps1
     * @param rmaps2
     * @return
     */
    protected boolean isSubgraph(List rmaps1, List rmaps2) {
        //System.out.println("Entering isSubgraph.");
        ArrayList rmaps2clone = (ArrayList) ((ArrayList) rmaps2).clone();
        for (Object o : rmaps1) {
            CDKRMap rmap1 = (CDKRMap) o;
            boolean found = false;
            for (int i = 0; i <
                    rmaps2clone.size(); ++i) {
                CDKRMap rmap2 = (CDKRMap) rmaps2clone.get(i);
                if (isSameRMap(rmap1, rmap2)) {
                    rmaps2clone.remove(i);
                    found =
                            true;
                    break;
                }

            }
            if (!found) {
                return false;
            }

        }
        return true;
    }

    /**
     *
     * @param sourceRMap sourceAtom
     * @param targetRMap targetAtom
     * @return
     */
    protected boolean isSameRMap(CDKRMap sourceRMap, CDKRMap targetRMap) {
        if (sourceRMap.getId1() == targetRMap.getId1() && sourceRMap.getId2() == targetRMap.getId2()) {
            return true;
        }

        return false;
    }

    public boolean getTimeOutFlag() {

        boolean flag = false;

        if (timeoutFlag) {
            flag = true;
        }

        return flag;
    }
}
