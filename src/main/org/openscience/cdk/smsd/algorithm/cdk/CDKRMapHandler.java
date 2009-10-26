/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.algorithm.cdk;


import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Stack;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.smsd.core.tools.EBIException;
import org.openscience.cdk.smsd.helper.FinalMappings;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.manipulator.BondManipulator;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 * This algorithm derives from the algorithm described in
 * [Tonnelier, C. and Jauffret, Ph. and Hanser, Th. and Jauffret, Ph. and Kaufmann, G.,
 * Machine Learning of generic reactions:
 * 3. An efficient algorithm for maximal common substructure determination,
 * Tetrahedron Comput. Methodol., 1990, 3:351-358] and modified in the thesis of
 * T. Hanser [Unknown BibTeXML type: HAN93].
 *
 */
public class CDKRMapHandler {

//    boolean RonPFlag = false;
    private Vector<TreeMap<Integer, Integer>> _mapping;
    static IAtomContainer ac1;
    static IAtomContainer ac2;
    private boolean timeoutFlag = false;

    /**
     * This function calculates all the possible combinations of MCS
     * @param Molecule1
     * @param Molecule2
     * @throws EBIException
     */
    public void calculateOverlapsAndReduce(IAtomContainer Molecule1, IAtomContainer Molecule2) throws EBIException {

        ac1 = Molecule1;
        ac2 = Molecule2;

        _mapping = new Vector<TreeMap<Integer, Integer>>();

        //System.out.println("Searching: ");
        //List overlaps = UniversalIsomorphismTesterBondTypeInSensitive.getSubgraphAtomsMap(ac1, ac2);



        if ((ac1.getAtomCount() == 1) || (ac2.getAtomCount() == 1)) {

            // System.out.println("Searched Single Molecule Mapping Case: ");

            List overlaps = CDKMCS.checkSingleAtomCases(ac1, ac2);
            //timeoutFlag=CDKMCS.getTimeOutFlag();


            int nAtomsMatched = overlaps.size();
            nAtomsMatched = (nAtomsMatched > 0) ? 1 : 0;
            if (nAtomsMatched > 0) {
                /*UnComment this to get one Unique Mapping*/
                //List reducedList = removeRedundantMappingsForSingleAtomCase(overlaps);
                //int counter = 0;
                identifyMatchedPartsforSingleAtoms(overlaps, ac1, ac2);

            }

        } else {

            /*Time out set by Asad equalt 2 min*/
            //UniversalIsomorphismTesterBondTypeInSensitive.timeout=12000000;

            //ac1 (Mol1), ac2 (Mol2), RGraph1, RGraph2, true, true( search all MCS)

            List overlaps = CDKMCS.search(ac1, ac2, new BitSet(), new BitSet(), true, true);


            //timeoutFlag=CDKMCS.getTimeOutFlag();
            //System.out.println("Searched: ");

            List reducedList = removeSubGraph(overlaps);


            Stack<List> allMaxOverlaps = getAllMaximum(reducedList);

            //Stack<List> allMaxOverlaps = getAllMaximum(overlaps);


            while (!allMaxOverlaps.empty()) {
//                System.out.println("ac1: " + ac1.getAtomCount() + ", ac2: " + ac2.getAtomCount() + ", overl: " + allMaxOverlaps.peek().size());
                @SuppressWarnings("unchecked")
                List maxOverlapsAtoms = makeAtomsMapOfBondsMap(allMaxOverlaps.peek(), ac1, ac2);
//                System.out.println("size of maxOverlaps: " + maxOverlapsAtoms.size());
                identifyMatchedParts(maxOverlapsAtoms, ac1, ac2);
//                identifyMatchedParts(allMaxOverlaps.peek(), ac1, ac2);
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
     * @throws EBIException
     */
    public void calculateOverlapsAndReduceExactMatch(IAtomContainer Molecule1, IAtomContainer Molecule2) throws EBIException {

        ac1 = Molecule1;
        ac2 = Molecule2;

        _mapping = new Vector<TreeMap<Integer, Integer>>();

        //System.out.println("Searching: ");
        //List overlaps = UniversalIsomorphismTesterBondTypeInSensitive.getSubgraphAtomsMap(ac1, ac2);



        if ((ac1.getAtomCount() == 1) || (ac2.getAtomCount() == 1)) {

            // System.out.println("Searched Single Molecule Mapping Case: ");

            List overlaps = CDKMCS.checkSingleAtomCases(ac1, ac2);
            //timeoutFlag=CDKMCS.getTimeOutFlag();


            int nAtomsMatched = overlaps.size();
            nAtomsMatched = (nAtomsMatched > 0) ? 1 : 0;
            if (nAtomsMatched > 0) {
                /*UnComment this to get one Unique Mapping*/
                //List reducedList = removeRedundantMappingsForSingleAtomCase(overlaps);
                //int counter = 0;
                identifyMatchedPartsforSingleAtoms(overlaps, ac1, ac2);

            }

        } else {

            /*Time out set by Asad equalt 2 min*/
            //UniversalIsomorphismTesterBondTypeInSensitive.timeout=12000000;

            //ac1 (Mol1), ac2 (Mol2), RGraph1, RGraph2, true, true( search all MCS)

            List overlaps = CDKMCS.search(ac1, ac2, new BitSet(), new BitSet(), true, true);
//            List overlaps = CDKMCS.search(ac1, ac2, new BitSet(), new BitSet(), true, false);
            //timeoutFlag=CDKMCS.getTimeOutFlag();
            //System.out.println("Searched: ");

            List reducedList = removeSubGraph(overlaps);


            Stack<List> allMaxOverlaps = getAllMaximum(reducedList);

            //Stack<List> allMaxOverlaps = getAllMaximum(overlaps);


            while (!allMaxOverlaps.empty()) {
//                System.out.println("ac1: " + ac1.getAtomCount() + ", ac2: " + ac2.getAtomCount() + ", overl: " + allMaxOverlaps.peek().size());
                @SuppressWarnings("unchecked")
                List maxOverlapsAtoms = makeAtomsMapOfBondsMap(allMaxOverlaps.peek(), ac1, ac2);
//                System.out.println("size of maxOverlaps: " + maxOverlapsAtoms.size());
                identifyMatchedParts(maxOverlapsAtoms, ac1, ac2);
//                identifyMatchedParts(allMaxOverlaps.peek(), ac1, ac2);

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
        List reducedList = new Vector(overlaps);

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

        List reducedList = new Vector();

        reducedList.add(overlaps.get(0));
        //reducedList.add(overlaps.get(1));

        return reducedList;

    }

    /**
     *  This makes a map of matching atoms out of a map of matching bonds as produced by the get(Subgraph|Ismorphism)Map methods.
     *
     * @param  l   The list produced by the getMap method.
     * @param  g1  first molecule. Must not be an IQueryAtomContainer.
     * @param  g2  second molecule. May be an IQueryAtomContainer.
     * @return     The mapping found projected on g1. This is a List of CDKRMap objects containing Ids of matching atoms.
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
//     *  This makes a map of matching atoms out of a map of matching bonds as produced by the get(Subgraph|Ismorphism)Map methods.
//     *
//     * @param  l   The list produced by the getMap method.
//     * @param  g1  first molecule. Must not be an IQueryAtomContainer.
//     * @param  g2  second molecule. May be an IQueryAtomContainer.
//     * @return     The mapping found projected on g1. This is a List of CDKRMap objects containing Ids of matching atoms.
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

            //System.out.println("O size" + a.size());

            if (a.size() > count) {
                @SuppressWarnings("unchecked")
                List list = new Vector(a);
                count = a.size();

                //System.out.println("List size" + list.size());

                //Collection threadSafeList = Collections.synchronizedCollection( list );
                allMaximumMappings = new Stack<List>();
                //allMaximumMappings.clear();
                allMaximumMappings.push(list);
            } else if (a.size() == count) {
                @SuppressWarnings("unchecked")
                List list = new Vector(a);
                count = a.size();
                allMaximumMappings.push(list);
            }

        }
        return allMaximumMappings;
    }

    /**
     *
     * @param list
     * @param ac1
     * @param ac2
     */
    protected void identifyMatchedParts(List list, IAtomContainer ac1, IAtomContainer ac2) {

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
            IAtom a = ac1.getAtom(rmap.getId1());
            IAtom b = ac2.getAtom(rmap.getId2());


            array1.add(a);
            array2.add(b);


//            System.err.print("Map " + a.getSymbol());
//            System.err.print("Map " + b.getSymbol());

            int IndexI = ac1.getAtomNumber(a);
            int IndexJ = ac2.getAtomNumber(b);

            atomNumbersFromContainer.put(IndexI, IndexJ);



            /*Added the Mapping Numbers to the FinalMapping*
             */
            _mapping.add(atomNumbersFromContainer);



        }

    }

    /**
     *
     * @param list
     * @param ac1
     * @param ac2
     */
    protected void identifyMatchedPartsforSingleAtoms(List list,
            IAtomContainer ac1,
            IAtomContainer ac2) {

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

            IAtom a = ac1.getAtom(rmap.getId1());
            IAtom b = ac2.getAtom(rmap.getId2());


            array1.add(a);
            array2.add(b);


//            System.err.print("Map " + a.getSymbol());
//            System.err.print("Map " + b.getSymbol());

            int IndexI = ac1.getAtomNumber(a);
            int IndexJ = ac2.getAtomNumber(b);


            atomNumbersFromContainer.put(IndexI, IndexJ);



            /*Added the Mapping Numbers to the FinalMapping*
             */
            _mapping.add(atomNumbersFromContainer);


        }



        //Asad: Uncomment this part of the code to remove the matched parts from the Molecules.

        /*for (IAtom atom : array1) {
        ac1.removeAtomAndConnectedElectronContainers(atom);
        }
        for (IAtom atom : array2) {
        ac2.removeAtomAndConnectedElectronContainers(atom);
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
     * @param a
     * @param b
     * @return
     */
    protected boolean isSameRMap(CDKRMap a, CDKRMap b) {
        if (a.getId1() == b.getId1() && a.getId2() == b.getId2()) {
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
