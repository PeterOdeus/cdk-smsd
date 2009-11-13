/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.core.tools;

import java.util.List;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond.Order;

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
public class EBIBondEnergies {

    private static TreeMap<Integer, List<Object>> bondEngergies = null;
    private static EBIBondEnergies _instance = null;
    private Integer key = null;

    /**
     * 
     * @return
     * @throws CDKException 
     */
    public synchronized static EBIBondEnergies getInstance()
            throws CDKException {
        if (null == _instance) {
            _instance = new EBIBondEnergies();
        }

        return _instance;
    }

    protected EBIBondEnergies() {
        key = 1;


        // =========Hydrogen Block==============

        bondEngergies = new TreeMap<Integer, List<Object>>();
        List<Object> v = new Vector<Object>();
        v.add("H");
        v.add("H");
        v.add(Order.SINGLE);
        v.add(432);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("H");
        v.add("B");
        v.add(Order.SINGLE);
        v.add(389);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("H");
        v.add("C");
        v.add(Order.SINGLE);
        v.add(411);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("H");
        v.add("Si");
        v.add(Order.SINGLE);
        v.add(318);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("H");
        v.add("Ge");
        v.add(Order.SINGLE);
        v.add(288);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("H");
        v.add("Sn");
        v.add(Order.SINGLE);
        v.add(251);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("H");
        v.add("N");
        v.add(Order.SINGLE);
        v.add(386);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("H");
        v.add("P");
        v.add(Order.SINGLE);
        v.add(322);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("H");
        v.add("As");
        v.add(Order.SINGLE);
        v.add(247);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("H");
        v.add("O");
        v.add(Order.SINGLE);
        v.add(459);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("H");
        v.add("S");
        v.add(Order.SINGLE);
        v.add(363);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("H");
        v.add("Se");
        v.add(Order.SINGLE);
        v.add(276);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("H");
        v.add("Te");
        v.add(Order.SINGLE);
        v.add(238);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("H");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(565);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("H");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(428);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("H");
        v.add("Br");
        v.add(Order.SINGLE);
        v.add(362);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("H");
        v.add("I");
        v.add(Order.SINGLE);
        v.add(295);
        bondEngergies.put(key++, v);


        //==================Group 13=================

        v = new Vector<Object>();
        v.add("B");
        v.add("B");
        v.add(Order.SINGLE);
        v.add(293);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("B");
        v.add("O");
        v.add(Order.SINGLE);
        v.add(536);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("B");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(613);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("B");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(456);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("B");
        v.add("Br");
        v.add(Order.SINGLE);
        v.add(377);
        bondEngergies.put(key++, v);


        //===================Group 14 Part 1=================

        v = new Vector<Object>();
        v.add("C");
        v.add("C");
        v.add(Order.SINGLE);
        v.add(346);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("C");
        v.add(Order.DOUBLE);
        v.add(602);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("C");
        v.add(Order.TRIPLE);
        v.add(835);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("Si");
        v.add(Order.SINGLE);
        v.add(318);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("Ge");
        v.add(Order.SINGLE);
        v.add(238);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("Sn");
        v.add(Order.SINGLE);
        v.add(130);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("N");
        v.add(Order.SINGLE);
        v.add(305);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("N");
        v.add(Order.DOUBLE);
        v.add(615);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("N");
        v.add(Order.TRIPLE);
        v.add(887);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("P");
        v.add(Order.SINGLE);
        v.add(264);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("O");
        v.add(Order.SINGLE);
        v.add(358);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("C");
        v.add("O");
        v.add(Order.DOUBLE);
        v.add(799);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("C");
        v.add("O");
        v.add(Order.TRIPLE);
        v.add(1072);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("B");
        v.add(Order.SINGLE);
        v.add(356);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("C");
        v.add("S");
        v.add(Order.SINGLE);
        v.add(272);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("C");
        v.add("S");
        v.add(Order.DOUBLE);
        v.add(573);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(485);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(327);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("Br");
        v.add(Order.SINGLE);
        v.add(285);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("C");
        v.add("I");
        v.add(Order.SINGLE);
        v.add(213);
        bondEngergies.put(key++, v);

        //===================Group 14 Part 2=================

        v = new Vector<Object>();
        v.add("Si");
        v.add("Si");
        v.add(Order.SINGLE);
        v.add(222);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Si");
        v.add("N");
        v.add(Order.SINGLE);
        v.add(355);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Si");
        v.add("O");
        v.add(Order.SINGLE);
        v.add(452);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Si");
        v.add("S");
        v.add(Order.SINGLE);
        v.add(293);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Si");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(222);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Si");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(565);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Si");
        v.add("Br");
        v.add(Order.SINGLE);
        v.add(310);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Si");
        v.add("I");
        v.add(Order.SINGLE);
        v.add(234);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Ge");
        v.add("Ge");
        v.add(Order.SINGLE);
        v.add(188);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Ge");
        v.add("N");
        v.add(Order.SINGLE);
        v.add(257);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Ge");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(470);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Ge");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(349);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Ge");
        v.add("Br");
        v.add(Order.SINGLE);
        v.add(276);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Ge");
        v.add("I");
        v.add(Order.SINGLE);
        v.add(212);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Sn");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(414);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Sn");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(323);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("Sn");
        v.add("Br");
        v.add(Order.SINGLE);
        v.add(273);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Sn");
        v.add("I");
        v.add(Order.SINGLE);
        v.add(205);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Pb");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(331);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Pb");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(243);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("Pb");
        v.add("Br");
        v.add(Order.SINGLE);
        v.add(201);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Pb");
        v.add("I");
        v.add(Order.SINGLE);
        v.add(142);
        bondEngergies.put(key++, v);


        v = new Vector<Object>();
        v.add("N");
        v.add("N");
        v.add(Order.SINGLE);
        v.add(167);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("N");
        v.add("N");
        v.add(Order.DOUBLE);
        v.add(418);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("N");
        v.add("N");
        v.add(Order.TRIPLE);
        v.add(942);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("N");
        v.add("O");
        v.add(Order.SINGLE);
        v.add(201);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("N");
        v.add("O");
        v.add(Order.DOUBLE);
        v.add(607);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("N");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(283);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("N");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(313);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("P");
        v.add("P");
        v.add(Order.SINGLE);
        v.add(201);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("P");
        v.add("O");
        v.add(Order.SINGLE);
        v.add(335);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("P");
        v.add("O");
        v.add(Order.DOUBLE);
        v.add(544);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("P");
        v.add("S");
        v.add(Order.DOUBLE);
        v.add(335);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("P");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(490);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("P");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(326);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("P");
        v.add("Br");
        v.add(Order.SINGLE);
        v.add(264);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("P");
        v.add("I");
        v.add(Order.SINGLE);
        v.add(184);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("As");
        v.add("As");
        v.add(Order.SINGLE);
        v.add(146);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("As");
        v.add("O");
        v.add(Order.SINGLE);
        v.add(301);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("As");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(484);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("As");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(322);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("As");
        v.add("Br");
        v.add(Order.SINGLE);
        v.add(458);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("As");
        v.add("I");
        v.add(Order.SINGLE);
        v.add(200);
        bondEngergies.put(key++, v);


        v = new Vector<Object>();
        v.add("Sb");
        v.add("Sb");
        v.add(Order.SINGLE);
        v.add(121);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Sb");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(440);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Sb");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(248);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Sb");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(315);
        bondEngergies.put(key++, v);



        //===================Group 16=================

        v = new Vector<Object>();
        v.add("O");
        v.add("O");
        v.add(Order.SINGLE);
        v.add(142);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("O");
        v.add("O");
        v.add(Order.DOUBLE);
        v.add(494);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("O");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(190);
        bondEngergies.put(key++, v);


        v = new Vector<Object>();
        v.add("S");
        v.add("O");
        v.add(Order.SINGLE);
        v.add(365);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("S");
        v.add("O");
        v.add(Order.DOUBLE);
        v.add(522);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("S");
        v.add("S");
        v.add(Order.SINGLE);
        v.add(226);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("S");
        v.add("S");
        v.add(Order.DOUBLE);
        v.add(425);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("S");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(284);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("S");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(255);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Se");
        v.add("Se");
        v.add(Order.SINGLE);
        v.add(172);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Se");
        v.add("Se");
        v.add(Order.DOUBLE);
        v.add(272);
        bondEngergies.put(key++, v);

        //===================Group 17=================

        v = new Vector<Object>();
        v.add("F");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(155);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Cl");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(240);
        bondEngergies.put(key++, v);
        v = new Vector<Object>();
        v.add("Br");
        v.add("Br");
        v.add(Order.SINGLE);
        v.add(190);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("I");
        v.add("I");
        v.add(Order.SINGLE);
        v.add(148);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("At");
        v.add("At");
        v.add(Order.SINGLE);
        v.add(116);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("I");
        v.add("O");
        v.add(Order.SINGLE);
        v.add(201);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("I");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(273);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("I");
        v.add("Cl");
        v.add(Order.SINGLE);
        v.add(208);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("I");
        v.add("Br");
        v.add(Order.SINGLE);
        v.add(175);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Se");
        v.add("Se");
        v.add(Order.SINGLE);
        v.add(272);
        bondEngergies.put(key++, v);

        //===================Group 18=================

        v = new Vector<Object>();
        v.add("Kr");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(50);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Xe");
        v.add("O");
        v.add(Order.SINGLE);
        v.add(84);
        bondEngergies.put(key++, v);

        v = new Vector<Object>();
        v.add("Xe");
        v.add("F");
        v.add(Order.SINGLE);
        v.add(130);
        bondEngergies.put(key++, v);
    }

    /**
     * 
     * @param a1 First atom
     * @param a2 Second atom
     * @param bondOrder (single, double etc)
     * @return bond energy
     */
    public Integer getEnergies(IAtom a1, IAtom a2, Order bondOrder) {
        Integer D_kJ_per_mol = -1;

        for (List<Object> V : bondEngergies.values()) {

            String atom1 = (String) V.get(0);
            String atom2 = (String) V.get(1);
            if (atom1.equalsIgnoreCase(a1.getSymbol()) && atom2.equalsIgnoreCase(a2.getSymbol())) {

                Order o = (Order) V.get(2);
                if (o.compareTo(bondOrder) == 0) {

                    D_kJ_per_mol = (Integer) V.get(3);
                }
            } else if (atom2.equalsIgnoreCase(a1.getSymbol()) && atom1.equalsIgnoreCase(a2.getSymbol())) {

                Order o = (Order) V.get(2);
                if (o.compareTo(bondOrder) == 0) {

                    D_kJ_per_mol = (Integer) V.get(3);
                }

            }

        }

        return D_kJ_per_mol;

    }
}
