/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.core.tools;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

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

public class SubGraphMapper {

    /**
     *
     * @param mainMolecule
     * @param matchedCore
     * @return
     */
    public static IAtomContainer getNeedle(IAtomContainer mainMolecule, List<Integer> matchedCore) {
        IAtomContainer needle = DefaultChemObjectBuilder.getInstance().newAtomContainer();

        assert matchedCore != null;
//        System.out.println("Number of matched atoms = " + matchedCore.size());
        for (int i = 0; i < matchedCore.size(); i++) {
            for (int j = 0; j < matchedCore.size(); j++) {

                int iIndex = matchedCore.get(i);
                int jIndex = matchedCore.get(j);

                if (i != j && iIndex != jIndex) {


                    IAtom atom0 = mainMolecule.getAtom(iIndex);
                    IAtom atom1 = mainMolecule.getAtom(jIndex);

                    atom0.setID(String.valueOf(iIndex));
                    atom1.setID(String.valueOf(jIndex));

                    IBond bond = mainMolecule.getBond(atom0, atom1);
                    if (bond != null) {
                        needle.addBond(bond);

                    }
                }
            }
        }

        return needle;
    }

    /**
     * 
     * @param matchedCore
     * @return
     */
    public static List<List<Integer>> getNeedleList(List<Integer> matchedCore) {
        List<List<Integer>> matchedList = new Vector<List<Integer>>();
        List<Integer> rList = new LinkedList<Integer>();
        List<Integer> pList = new LinkedList<Integer>();

        for (int i = 0; i < matchedCore.size(); i += 2) {

            rList.add(matchedCore.get(i));
            pList.add(matchedCore.get(i + 1));

        }

        matchedList.add(0, rList);
        matchedList.add(1, pList);
        return matchedList;
    }

    /**
     * 
     * @param matchedCore
     * @return
     */
    public static List<List<Integer>> getNeedleList(Map<Integer, Integer> matchedCore) {
        List<List<Integer>> matchedList = new Vector<List<Integer>>();
        List<Integer> rList = new LinkedList<Integer>();
        List<Integer> pList = new LinkedList<Integer>();

        for (Map.Entry<Integer, Integer> M : matchedCore.entrySet()) {

            rList.add(M.getKey());
            pList.add(M.getValue());

        }

        matchedList.add(0, rList);
        matchedList.add(1, pList);
        return matchedList;
    }

    /**
     * 
     * @param ac1
     * @param ac2
     * @param mappings
     * @return
     */
    public static Map<IBond, IBond> makeBondMapsOfAtomMaps(IAtomContainer ac1, IAtomContainer ac2, TreeMap<Integer, Integer> mappings) {

        HashMap<IBond, IBond> maps = new HashMap<IBond, IBond>();

        for (Map.Entry<Integer, Integer> mapS : mappings.entrySet()) {
            int indexI = mapS.getKey();
            int indexJ = mapS.getValue();

            for (Map.Entry<Integer, Integer> mapD : mappings.entrySet()) {

                int indexIPlus = mapD.getKey();
                int indexJPlus = mapD.getValue();

                if (indexI != indexIPlus && indexJ != indexJPlus) {

                    IAtom atomI0 = ac1.getAtom(indexI);
                    IAtom atomI1 = ac1.getAtom(indexIPlus);

                    atomI0.setID(String.valueOf(indexI));
                    atomI1.setID(String.valueOf(indexIPlus));

                    IBond ac1Bond = ac1.getBond(atomI0, atomI1);

                    if (ac1Bond != null) {

                        IAtom atomJ0 = ac1.getAtom(indexJ);
                        IAtom atomJ1 = ac1.getAtom(indexJPlus);

                        atomJ0.setID(String.valueOf(indexJ));
                        atomJ1.setID(String.valueOf(indexJPlus));

                        IBond ac2Bond = ac1.getBond(atomJ0, atomJ1);

                        if (ac2Bond != null) {

                            maps.put(ac1Bond, ac2Bond);
                        }

                    }

                }

            }
        }

//        System.out.println("bond Map size:" + maps.size());

        return maps;

    }
}
