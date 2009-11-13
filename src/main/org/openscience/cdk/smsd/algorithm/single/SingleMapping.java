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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.smsd.algorithm.single;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.smsd.helper.FinalMappings;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;

/**
 * @cdk.module smsd
 */
public class SingleMapping {

    private IAtomContainer ac1 = null;
    private IAtomContainer ac2 = null;

       /**
     *
     * @param Molecule1
     * @param Molecule2
     */
    private void printMolecules(IAtomContainer Molecule1, IAtomContainer Molecule2) {

        System.out.println("Molecule 1: " + Molecule1.getAtomCount());

        for (int i = 0; i < Molecule1.getAtomCount(); i++) {

            System.out.print(Molecule1.getAtom(i).getSymbol() + " ");
        }

        System.out.println();
        System.out.println("Molecule 2: " + Molecule2.getAtomCount());
        for (int i = 0; i < Molecule2.getAtomCount(); i++) {

            System.out.print(Molecule2.getAtom(i).getSymbol() + " ");
        }
        System.out.println();

    }

    public void getOverLaps(IAtomContainer ac1, IAtomContainer ac2, boolean removeHydrogen) {

        this.ac1 = ac1;
        this.ac2 = ac2;

//        printMolecules(ac1, ac2);

        Vector<TreeMap<Integer, Integer>> _mappings = new Vector<TreeMap<Integer, Integer>>();
        Map<Integer, Integer> connectedBondOrder = new TreeMap<Integer, Integer>();
        //Vector<Vector<Integer>> _mappingFiltered = new Vector<Vector<Integer>>();

        //System.out.println("Hello in Single getOverlaps");

        int minOrder = 9999;

        if (removeHydrogen) {
            mapWithoutH(_mappings, connectedBondOrder, minOrder);
        } else {
            mapWithH(_mappings, connectedBondOrder, minOrder);
        }


        postFilter(_mappings, connectedBondOrder, minOrder);

        _mappings.clear();
    }

    private void mapWithoutH(Vector<TreeMap<Integer, Integer>> _mappings, Map<Integer, Integer> BondOrder, int minOrder) {
        int counter = 0;

        if ((ac1.getAtomCount() == 1) && (!ac1.getAtom(0).getSymbol().equals("H"))) {
            //System.out.println("Hello in Single getOverlaps-> Cond1");
            for (int i = 0; i < ac2.getAtomCount(); i++) {


                TreeMap<Integer, Integer> mapAtoms = new TreeMap<Integer, Integer>();

                if (ac1.getAtom(0).getSymbol().equalsIgnoreCase(ac2.getAtom(i).getSymbol())) {
                    //System.out.println("Atom 1: " + ac1.getAtom(0).getSymbol());
                    //System.out.println("Atom 2: " + ac2.getAtom(i).getSymbol());
                    mapAtoms.put(0, i);


                    IAtom atom = ac2.getAtom(i);
                    List<IBond> Bonds = ac2.getConnectedBondsList(atom);

                    int totalOrder = 0;
                    for (IBond B : Bonds) {

                        Order BO = B.getOrder();
                        totalOrder += BO.ordinal();
                    }

                    if (totalOrder < minOrder) {
                        minOrder = totalOrder;
                    }

                    BondOrder.put(counter, totalOrder);
                    _mappings.add(counter++, mapAtoms);

                    //_mappings.addElement(mapAtoms);


                }

                //System.out.println("Hello in Single getOverlaps Mapping Size: " + mapAtoms.size());
            }
        } else if ((ac2.getAtomCount() == 1) && (!ac2.getAtom(0).getSymbol().equals("H"))) {
            //System.out.println("Hello in Single getOverlaps-> Cond2");


            for (int i = 0; i < ac1.getAtomCount(); i++) {
                TreeMap<Integer, Integer> mapAtoms = new TreeMap<Integer, Integer>();

                if (ac2.getAtom(0).getSymbol().equalsIgnoreCase(ac1.getAtom(i).getSymbol())) {
                    mapAtoms.put(i, 0);

                    IAtom atom = ac1.getAtom(i);
                    List<IBond> Bonds = ac1.getConnectedBondsList(atom);

                    int totalOrder = 0;
                    for (IBond B : Bonds) {

                        Order BO = B.getOrder();
                        totalOrder += BO.ordinal();
                    }
                    if (totalOrder < minOrder) {
                        minOrder = totalOrder;
                    }

                    BondOrder.put(counter, totalOrder);
                    _mappings.add(counter++, mapAtoms);
                    //_mappings.addElement(mapAtoms);

                }

//                System.out.println("Hello in Single getOverlaps Mapping Size: " + mapAtoms.size());
            }

        } else {
            System.err.println("skippping Hydrogen mapping");
//            try {
//                
//                throw new CDKException("Hydrogens can't be Mapped");
//            } catch (CDKException ex) {
//                Logger.getLogger(SingleMapping.class.getName()).log(Level.SEVERE, null, ex);
//            }
        }
    }

    private void mapWithH(Vector<TreeMap<Integer, Integer>> _mappings, Map<Integer, Integer> BondOrder, int minOrder) {
        int counter = 0;

        if (ac1.getAtomCount() == 1) {
            //System.out.println("Hello in Single getOverlaps-> Cond1");
            for (int i = 0; i < ac2.getAtomCount(); i++) {


                TreeMap<Integer, Integer> mapAtoms = new TreeMap<Integer, Integer>();

                if (ac1.getAtom(0).getSymbol().equalsIgnoreCase(ac2.getAtom(i).getSymbol())) {
                    //System.out.println("Atom 1: " + ac1.getAtom(0).getSymbol());
                    //System.out.println("Atom 2: " + ac2.getAtom(i).getSymbol());
                    mapAtoms.put(0, i);


                    IAtom atom = ac2.getAtom(i);
                    List<IBond> Bonds = ac2.getConnectedBondsList(atom);

                    int totalOrder = 0;
                    for (IBond B : Bonds) {

                        Order BO = B.getOrder();
                        totalOrder += BO.ordinal();
                    }

                    if (totalOrder < minOrder) {
                        minOrder = totalOrder;
                    }

                    BondOrder.put(counter, totalOrder);
                    _mappings.add(counter++, mapAtoms);

                    //_mappings.addElement(mapAtoms);


                }

//                System.out.println("Hello in Single getOverlaps Mapping Size: " + mapAtoms.size());
            }
        } else if (ac2.getAtomCount() == 1) {
            //System.out.println("Hello in Single getOverlaps-> Cond2");


            for (int i = 0; i < ac1.getAtomCount(); i++) {
                TreeMap<Integer, Integer> mapAtoms = new TreeMap<Integer, Integer>();

                if (ac2.getAtom(0).getSymbol().equalsIgnoreCase(ac1.getAtom(i).getSymbol())) {
                    mapAtoms.put(i, 0);

                    IAtom atom = ac1.getAtom(i);
                    List<IBond> Bonds = ac1.getConnectedBondsList(atom);

                    int totalOrder = 0;
                    for (IBond B : Bonds) {

                        Order BO = B.getOrder();
                        totalOrder += BO.ordinal();
                    }
                    if (totalOrder < minOrder) {
                        minOrder = totalOrder;
                    }

                    BondOrder.put(counter, totalOrder);
                    _mappings.add(counter++, mapAtoms);
                    //_mappings.addElement(mapAtoms);

                }

//                System.out.println("Hello in Single getOverlaps Mapping Size: " + mapAtoms.size());
            }

        } else {
            System.err.println("This is not a single mapping case!");

        }
    }

    private void postFilter(Vector<TreeMap<Integer, Integer>> _mapping, Map<Integer, Integer> BondOrder, int minOrder) {

        for (Map.Entry<Integer, Integer> map : BondOrder.entrySet()) {

            if (map.getValue() > minOrder) {
                removedMap(_mapping, map.getKey());
            }

        }

//        System.out.println("Hello in Single getOverlaps Mapping Size: " + _mapping.firstElement().size());
        FinalMappings final_MAPPINGS = FinalMappings.getInstance();

        final_MAPPINGS.set(new Vector<TreeMap<Integer, Integer>>(_mapping));

//        System.out.println("Hello in Single getOverlaps Final Mapping Size: " + final_MAPPINGS.getSize());


    }

    private boolean removedMap(Vector<TreeMap<Integer, Integer>> maps, Integer Key) {
        boolean flag = false;
        for (TreeMap<Integer, Integer> map : maps) {

            if (map.containsKey(Key)) {

                flag = maps.remove(map);
            }
        }
        return flag;
    }
}
