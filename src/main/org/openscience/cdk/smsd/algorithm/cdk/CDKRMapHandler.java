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

/**
 * This algorithm derives from the algorithm described in
 * [Tonnelier, C. and Jauffret, Ph. and Hanser,
 * Th. and Jauffret, Ph. and Kaufmann, G.,
 * Machine Learning of generic reactions:
 * 3. An efficient algorithm for maximal common substructure determination,
 * Tetrahedron Comput. Methodol., 1990, 3:351-358] and modified in the thesis of
 * T. Hanser [Unknown BibTeXML type: HAN93].
 *
 * @cdk.module smsd
 */
public class CDKRMapHandler {

    private List<TreeMap<Integer, Integer>> mapping;
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

        mapping = new ArrayList<TreeMap<Integer, Integer>>();


        if ((source.getAtomCount() == 1) || (target.getAtomCount() == 1)) {
            List<CDKRMap> overlaps = CDKMCS.checkSingleAtomCases(source, target);
            int nAtomsMatched = overlaps.size();
            nAtomsMatched = (nAtomsMatched > 0) ? 1 : 0;
            if (nAtomsMatched > 0) {
                /*UnComment this to get one Unique Mapping*/
                //List reducedList = removeRedundantMappingsForSingleAtomCase(overlaps);
                //int counter = 0;
                identifyMatchedPartsforSingleAtoms(overlaps, source, target);

            }

        } else {
            List<List<CDKRMap>> overlaps = CDKMCS.search(source, target, new BitSet(), new BitSet(), true, true);

            List<List<CDKRMap>> reducedList = removeSubGraph(overlaps);
            Stack<List<CDKRMap>> allMaxOverlaps = getAllMaximum(reducedList);
            while (!allMaxOverlaps.empty()) {
//                System.out.println("source: " + source.getAtomCount() + ", target: " + target.getAtomCount() + ", overl: " + allMaxOverlaps.peek().size());
                List<CDKRMap> maxOverlapsAtoms = CDKMCS.makeAtomsMapOfBondsMap(allMaxOverlaps.peek(), source, target);
//                System.out.println("size of maxOverlaps: " + maxOverlapsAtoms.size());
                identifyMatchedParts(maxOverlapsAtoms, source, target);
//                identifyMatchedParts(allMaxOverlaps.peek(), source, target);
                allMaxOverlaps.pop();
            }
        }

        FinalMappings.getInstance().set(mapping);

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

        mapping = new ArrayList<TreeMap<Integer, Integer>>();

        //System.out.println("Searching: ");
        //List overlaps = UniversalIsomorphismTesterBondTypeInSensitive.getSubgraphAtomsMap(source, target);

        if ((source.getAtomCount() == 1) || (target.getAtomCount() == 1)) {

            List<CDKRMap> overlaps = CDKMCS.checkSingleAtomCases(source, target);
            int nAtomsMatched = overlaps.size();
            nAtomsMatched = (nAtomsMatched > 0) ? 1 : 0;
            if (nAtomsMatched > 0) {
                identifyMatchedPartsforSingleAtoms(overlaps, source, target);
            }

        } else {

            List<List<CDKRMap>> overlaps = CDKMCS.search(source, target, new BitSet(), new BitSet(), true, true);

            List<List<CDKRMap>> reducedList = removeSubGraph(overlaps);
            Stack<List<CDKRMap>> allMaxOverlaps = getAllMaximum(reducedList);

            while (!allMaxOverlaps.empty()) {
                List<CDKRMap> maxOverlapsAtoms = CDKMCS.makeAtomsMapOfBondsMap(allMaxOverlaps.peek(), source, target);
                identifyMatchedParts(maxOverlapsAtoms, source, target);
                allMaxOverlaps.pop();
            }
        }
        FinalMappings.getInstance().set(mapping);
    }

    /**
     *
     * @param overlaps
     * @return
     */
    protected List<List<CDKRMap>> removeSubGraph(List<List<CDKRMap>> overlaps) {

        List<List<CDKRMap>> reducedList = new ArrayList<List<CDKRMap>>(overlaps);

        for (int i = 0; i < overlaps.size(); i++) {
            List<CDKRMap> graphI = overlaps.get(i);

            for (int j = i + 1; j < overlaps.size(); j++) {

                List<CDKRMap> graphJ = overlaps.get(j);

                // Gi included in Gj or Gj included in Gi then
                // reduce the irrelevant solution
                if (graphI.size() != graphJ.size()) {
                    if (isSubgraph(graphJ, graphI)) {
                        reducedList.remove(graphI);
                    } else if (isSubgraph(graphI, graphJ)) {
                        reducedList.remove(graphJ);
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
    protected List<CDKRMap> removeRedundantMappingsForSingleAtomCase(List<CDKRMap> overlaps) {
        List<CDKRMap> reducedList = new ArrayList<CDKRMap>();
        reducedList.add(overlaps.get(0));
        //reducedList.add(overlaps.get(1));
        return reducedList;

    }

    /**
     *
     * @param overlaps
     * @return
     */
    protected List<CDKRMap> getMaximum(List<List<CDKRMap>> overlaps) {


        List<CDKRMap> list = null;
        int count = 0;
        for (List<CDKRMap> arrayList : overlaps) {
            if (arrayList.size() > count) {
                list = arrayList;
                count = arrayList.size();
            }
        }
        return list;
    }

    /**
     *
     * @param overlaps
     * @return
     */
    protected Stack<List<CDKRMap>> getAllMaximum(List<List<CDKRMap>> overlaps) {

        Stack<List<CDKRMap>> allMaximumMappings = null;

        int count = -1;

        for (List<CDKRMap> arrayList : overlaps) {

            //System.out.println("O size" + sourceAtom.size());

            if (arrayList.size() > count) {

                List<CDKRMap> list = new ArrayList<CDKRMap>(arrayList);
                count = arrayList.size();

                //System.out.println("List size" + list.size());

                //Collection threadSafeList = Collections.synchronizedCollection( list );
                allMaximumMappings = new Stack<List<CDKRMap>>();
                //allMaximumMappings.clear();
                allMaximumMappings.push(list);
            } else if (arrayList.size() == count) {

                List<CDKRMap> list = new ArrayList<CDKRMap>(arrayList);
                count = arrayList.size();
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
    protected void identifyMatchedParts(List<CDKRMap> list, IAtomContainer source, IAtomContainer target) {

        List<IAtom> array1 = new ArrayList<IAtom>();
        List<IAtom> array2 = new ArrayList<IAtom>();

        TreeMap<Integer, Integer> atomNumbersFromContainer = new TreeMap<Integer, Integer>();


        /*We have serial numbers of the bonds/Atoms to delete
         * Now we will collect the actual bond/Atoms rather than
         * serial number for deletion. RonP flag check whether Reactant is
         * mapped on Product or Vise Versa
         */



        for (CDKRMap rmap : list) {
            IAtom sourceAtom = source.getAtom(rmap.getId1());
            IAtom targetAtom = target.getAtom(rmap.getId2());


            array1.add(sourceAtom);
            array2.add(targetAtom);

            int IndexI = source.getAtomNumber(sourceAtom);
            int IndexJ = target.getAtomNumber(targetAtom);

            atomNumbersFromContainer.put(IndexI, IndexJ);

            /*Added the Mapping Numbers to the FinalMapping*/
            mapping.add(atomNumbersFromContainer);
        }
    }

    /**
     *
     * @param list
     * @param source
     * @param target
     */
    protected void identifyMatchedPartsforSingleAtoms(List<CDKRMap> list,
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



        for (CDKRMap rmap : list) {

            IAtom sAtom = source.getAtom(rmap.getId1());
            IAtom tAtom = target.getAtom(rmap.getId2());


            array1.add(sAtom);
            array2.add(tAtom);

            int IndexI = source.getAtomNumber(sAtom);
            int IndexJ = target.getAtomNumber(tAtom);


            atomNumbersFromContainer.put(IndexI, IndexJ);

            /*Added the Mapping Numbers to the FinalMapping*/
            mapping.add(atomNumbersFromContainer);


        }
    }

    /**
     *
     * @param rmaps1
     * @param rmaps2
     * @return
     */
    protected boolean isSubgraph(List<CDKRMap> rmaps1, List<CDKRMap> rmaps2) {
        //System.out.println("Entering isSubgraph.");
        List<CDKRMap> rmaps2clone = new ArrayList<CDKRMap>(rmaps2);
        for (CDKRMap rmap1 : rmaps1) {
            boolean found = false;
            for (int i = 0; i < rmaps2clone.size(); ++i) {
                CDKRMap rmap2 = rmaps2clone.get(i);
                if (isSameRMap(rmap1, rmap2)) {
                    rmaps2clone.remove(i);
                    found = true;
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

    /**
     *
     * @return time out flag
     */
    public boolean getTimeOutFlag() {

        boolean flag = false;

        if (timeoutFlag) {
            flag = true;
        }

        return flag;
    }
}
