/* Copyright (C) 2006-2009  Syed Asad Rahman <asad@ebi.ac.uk>
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
package org.openscience.cdk.smsd.tools;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond.Order;

/**
 * @cdk.module smsd
 */
@TestClass("org.openscience.cdk.smsd.tools.BondEnergiesTest")
public class BondEnergies {

    private static TreeMap<Integer, List<Object>> bondEngergies = null;
    private static BondEnergies _instance = null;
    private Integer key = null;

    /**
     * 
     * @return
     * @throws CDKException 
     */
    @TestMethod("testGetInstance")
    public synchronized static BondEnergies getInstance()
            throws CDKException {
        if (null == _instance) {
            _instance = new BondEnergies();
        }

        return _instance;
    }

    protected BondEnergies() {
        key = 1;
        // =========Hydrogen Block==============

        bondEngergies = new TreeMap<Integer, List<Object>>();
        List<Object> element = new ArrayList<Object>();
        element.add("H");
        element.add("H");
        element.add(Order.SINGLE);
        element.add(432);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("H");
        element.add("B");
        element.add(Order.SINGLE);
        element.add(389);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("H");
        element.add("C");
        element.add(Order.SINGLE);
        element.add(411);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("H");
        element.add("Si");
        element.add(Order.SINGLE);
        element.add(318);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("H");
        element.add("Ge");
        element.add(Order.SINGLE);
        element.add(288);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("H");
        element.add("Sn");
        element.add(Order.SINGLE);
        element.add(251);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("H");
        element.add("N");
        element.add(Order.SINGLE);
        element.add(386);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("H");
        element.add("P");
        element.add(Order.SINGLE);
        element.add(322);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("H");
        element.add("As");
        element.add(Order.SINGLE);
        element.add(247);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("H");
        element.add("O");
        element.add(Order.SINGLE);
        element.add(459);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("H");
        element.add("S");
        element.add(Order.SINGLE);
        element.add(363);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("H");
        element.add("Se");
        element.add(Order.SINGLE);
        element.add(276);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("H");
        element.add("Te");
        element.add(Order.SINGLE);
        element.add(238);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("H");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(565);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("H");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(428);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("H");
        element.add("Br");
        element.add(Order.SINGLE);
        element.add(362);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("H");
        element.add("I");
        element.add(Order.SINGLE);
        element.add(295);
        bondEngergies.put(key++, element);


        //==================Group 13=================

        element = new ArrayList<Object>();
        element.add("B");
        element.add("B");
        element.add(Order.SINGLE);
        element.add(293);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("B");
        element.add("O");
        element.add(Order.SINGLE);
        element.add(536);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("B");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(613);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("B");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(456);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("B");
        element.add("Br");
        element.add(Order.SINGLE);
        element.add(377);
        bondEngergies.put(key++, element);


        //===================Group 14 Part 1=================

        element = new ArrayList<Object>();
        element.add("C");
        element.add("C");
        element.add(Order.SINGLE);
        element.add(346);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("C");
        element.add(Order.DOUBLE);
        element.add(602);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("C");
        element.add(Order.TRIPLE);
        element.add(835);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("Si");
        element.add(Order.SINGLE);
        element.add(318);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("Ge");
        element.add(Order.SINGLE);
        element.add(238);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("Sn");
        element.add(Order.SINGLE);
        element.add(130);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("N");
        element.add(Order.SINGLE);
        element.add(305);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("N");
        element.add(Order.DOUBLE);
        element.add(615);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("N");
        element.add(Order.TRIPLE);
        element.add(887);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("P");
        element.add(Order.SINGLE);
        element.add(264);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("O");
        element.add(Order.SINGLE);
        element.add(358);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("C");
        element.add("O");
        element.add(Order.DOUBLE);
        element.add(799);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("C");
        element.add("O");
        element.add(Order.TRIPLE);
        element.add(1072);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("B");
        element.add(Order.SINGLE);
        element.add(356);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("C");
        element.add("S");
        element.add(Order.SINGLE);
        element.add(272);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("C");
        element.add("S");
        element.add(Order.DOUBLE);
        element.add(573);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(485);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(327);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("Br");
        element.add(Order.SINGLE);
        element.add(285);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("C");
        element.add("I");
        element.add(Order.SINGLE);
        element.add(213);
        bondEngergies.put(key++, element);

        //===================Group 14 Part 2=================

        element = new ArrayList<Object>();
        element.add("Si");
        element.add("Si");
        element.add(Order.SINGLE);
        element.add(222);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Si");
        element.add("N");
        element.add(Order.SINGLE);
        element.add(355);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Si");
        element.add("O");
        element.add(Order.SINGLE);
        element.add(452);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Si");
        element.add("S");
        element.add(Order.SINGLE);
        element.add(293);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Si");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(222);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Si");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(565);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Si");
        element.add("Br");
        element.add(Order.SINGLE);
        element.add(310);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Si");
        element.add("I");
        element.add(Order.SINGLE);
        element.add(234);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Ge");
        element.add("Ge");
        element.add(Order.SINGLE);
        element.add(188);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Ge");
        element.add("N");
        element.add(Order.SINGLE);
        element.add(257);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Ge");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(470);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Ge");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(349);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Ge");
        element.add("Br");
        element.add(Order.SINGLE);
        element.add(276);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Ge");
        element.add("I");
        element.add(Order.SINGLE);
        element.add(212);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Sn");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(414);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Sn");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(323);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("Sn");
        element.add("Br");
        element.add(Order.SINGLE);
        element.add(273);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Sn");
        element.add("I");
        element.add(Order.SINGLE);
        element.add(205);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Pb");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(331);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Pb");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(243);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("Pb");
        element.add("Br");
        element.add(Order.SINGLE);
        element.add(201);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Pb");
        element.add("I");
        element.add(Order.SINGLE);
        element.add(142);
        bondEngergies.put(key++, element);


        element = new ArrayList<Object>();
        element.add("N");
        element.add("N");
        element.add(Order.SINGLE);
        element.add(167);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("N");
        element.add("N");
        element.add(Order.DOUBLE);
        element.add(418);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("N");
        element.add("N");
        element.add(Order.TRIPLE);
        element.add(942);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("N");
        element.add("O");
        element.add(Order.SINGLE);
        element.add(201);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("N");
        element.add("O");
        element.add(Order.DOUBLE);
        element.add(607);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("N");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(283);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("N");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(313);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("P");
        element.add("P");
        element.add(Order.SINGLE);
        element.add(201);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("P");
        element.add("O");
        element.add(Order.SINGLE);
        element.add(335);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("P");
        element.add("O");
        element.add(Order.DOUBLE);
        element.add(544);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("P");
        element.add("S");
        element.add(Order.DOUBLE);
        element.add(335);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("P");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(490);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("P");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(326);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("P");
        element.add("Br");
        element.add(Order.SINGLE);
        element.add(264);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("P");
        element.add("I");
        element.add(Order.SINGLE);
        element.add(184);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("As");
        element.add("As");
        element.add(Order.SINGLE);
        element.add(146);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("As");
        element.add("O");
        element.add(Order.SINGLE);
        element.add(301);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("As");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(484);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("As");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(322);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("As");
        element.add("Br");
        element.add(Order.SINGLE);
        element.add(458);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("As");
        element.add("I");
        element.add(Order.SINGLE);
        element.add(200);
        bondEngergies.put(key++, element);


        element = new ArrayList<Object>();
        element.add("Sb");
        element.add("Sb");
        element.add(Order.SINGLE);
        element.add(121);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Sb");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(440);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Sb");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(248);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Sb");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(315);
        bondEngergies.put(key++, element);



        //===================Group 16=================

        element = new ArrayList<Object>();
        element.add("O");
        element.add("O");
        element.add(Order.SINGLE);
        element.add(142);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("O");
        element.add("O");
        element.add(Order.DOUBLE);
        element.add(494);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("O");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(190);
        bondEngergies.put(key++, element);


        element = new ArrayList<Object>();
        element.add("S");
        element.add("O");
        element.add(Order.SINGLE);
        element.add(365);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("S");
        element.add("O");
        element.add(Order.DOUBLE);
        element.add(522);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("S");
        element.add("S");
        element.add(Order.SINGLE);
        element.add(226);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("S");
        element.add("S");
        element.add(Order.DOUBLE);
        element.add(425);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("S");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(284);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("S");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(255);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Se");
        element.add("Se");
        element.add(Order.SINGLE);
        element.add(172);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Se");
        element.add("Se");
        element.add(Order.DOUBLE);
        element.add(272);
        bondEngergies.put(key++, element);

        //===================Group 17=================

        element = new ArrayList<Object>();
        element.add("F");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(155);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Cl");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(240);
        bondEngergies.put(key++, element);
        element = new ArrayList<Object>();
        element.add("Br");
        element.add("Br");
        element.add(Order.SINGLE);
        element.add(190);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("I");
        element.add("I");
        element.add(Order.SINGLE);
        element.add(148);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("At");
        element.add("At");
        element.add(Order.SINGLE);
        element.add(116);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("I");
        element.add("O");
        element.add(Order.SINGLE);
        element.add(201);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("I");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(273);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("I");
        element.add("Cl");
        element.add(Order.SINGLE);
        element.add(208);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("I");
        element.add("Br");
        element.add(Order.SINGLE);
        element.add(175);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Se");
        element.add("Se");
        element.add(Order.SINGLE);
        element.add(272);
        bondEngergies.put(key++, element);

        //===================Group 18=================

        element = new ArrayList<Object>();
        element.add("Kr");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(50);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Xe");
        element.add("O");
        element.add(Order.SINGLE);
        element.add(84);
        bondEngergies.put(key++, element);

        element = new ArrayList<Object>();
        element.add("Xe");
        element.add("F");
        element.add(Order.SINGLE);
        element.add(130);
        bondEngergies.put(key++, element);
    }

    /**
     * 
     * @param sourceAtom First atom
     * @param targetAtom Second atom
     * @param bondOrder (single, double etc)
     * @return bond energy
     */
    public Integer getEnergies(IAtom sourceAtom, IAtom targetAtom, Order bondOrder) {
        Integer D_kJ_per_mol = -1;

        for (List<Object> atom : bondEngergies.values()) {

            String atom1 = (String) atom.get(0);
            String atom2 = (String) atom.get(1);
            if (
                    (atom1.equalsIgnoreCase(sourceAtom.getSymbol())
                    && atom2.equalsIgnoreCase(targetAtom.getSymbol()))
                    || (atom2.equalsIgnoreCase(sourceAtom.getSymbol())
                    && atom1.equalsIgnoreCase(targetAtom.getSymbol()))){

                Order order = (Order) atom.get(2);
                if (order.compareTo(bondOrder) == 0) {

                    D_kJ_per_mol = (Integer) atom.get(3);
                }
            }

        }

        return D_kJ_per_mol;

    }
}
