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

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.smsd.helper.FinalMappings;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;

/**
 * @cdk.module smsd
 */
public class SingleMapping {

    private IAtomContainer source = null;
    private IAtomContainer target = null;

    /**
     *
     * @param source
     * @param target
     * @param removeHydrogen
     */
    protected void getOverLaps(IAtomContainer source, IAtomContainer target, boolean removeHydrogen) {

        this.source = source;
        this.target = target;
        List<TreeMap<Integer, Integer>> _mappings = new ArrayList<TreeMap<Integer, Integer>>();
        Map<Integer, Integer> connectedBondOrder = new TreeMap<Integer, Integer>();

        int minOrder = 9999;

        if (source.getAtomCount() == 1) {
            setSourceSingleAtomMap(removeHydrogen, _mappings, connectedBondOrder, minOrder);
        }
        if (target.getAtomCount() == 1) {
            setTargetSingleAtomMap(removeHydrogen, _mappings, connectedBondOrder, minOrder);
        }
        postFilter(_mappings, connectedBondOrder, minOrder);

        _mappings.clear();
    }

    private void setSourceSingleAtomMap(boolean removeHydrogen, List<TreeMap<Integer, Integer>> _mappings, Map<Integer, Integer> BondOrder, int minOrder) {
        int counter = 0;

        if ((removeHydrogen && !source.getAtom(0).getSymbol().equals("H")) || (!removeHydrogen)) {
            for (int i = 0; i < target.getAtomCount(); i++) {
                TreeMap<Integer, Integer> mapAtoms = new TreeMap<Integer, Integer>();

                if (source.getAtom(0).getSymbol().equalsIgnoreCase(target.getAtom(i).getSymbol())) {
                    mapAtoms.put(0, i);
                    IAtom atom = target.getAtom(i);
                    List<IBond> Bonds = target.getConnectedBondsList(atom);

                    int totalOrder = 0;
                    for (IBond bond : Bonds) {

                        Order bondOrder = bond.getOrder();
                        totalOrder += bondOrder.ordinal();
                    }

                    if (totalOrder < minOrder) {
                        minOrder = totalOrder;
                    }

                    BondOrder.put(counter, totalOrder);
                    _mappings.add(counter++, mapAtoms);
                }

                //System.out.println("Hello in Single getOverlaps Mapping Size: " + mapAtoms.size());
            }
        } else {
            System.err.println("Skippping Hydrogen mapping or This is not a single mapping case!");
        }
    }

    private void setTargetSingleAtomMap(boolean removeHydrogen, List<TreeMap<Integer, Integer>> _mappings, Map<Integer, Integer> BondOrder, int minOrder) {
        int counter = 0;

        if ((removeHydrogen && !target.getAtom(0).getSymbol().equals("H")) || (!removeHydrogen)) {
            for (int i = 0; i < source.getAtomCount(); i++) {
                TreeMap<Integer, Integer> mapAtoms = new TreeMap<Integer, Integer>();

                if (target.getAtom(0).getSymbol().equalsIgnoreCase(source.getAtom(i).getSymbol())) {
                    mapAtoms.put(i, 0);

                    IAtom atom = source.getAtom(i);
                    List<IBond> Bonds = source.getConnectedBondsList(atom);

                    int totalOrder = 0;
                    for (IBond bond : Bonds) {

                        Order bondOrder = bond.getOrder();
                        totalOrder += bondOrder.ordinal();
                    }
                    if (totalOrder < minOrder) {
                        minOrder = totalOrder;
                    }

                    BondOrder.put(counter, totalOrder);
                    _mappings.add(counter++, mapAtoms);
                }

//                System.out.println("Hello in Single getOverlaps Mapping Size: " + mapAtoms.size());
            }

        } else {
            System.err.println("Skippping Hydrogen mapping or This is not a single mapping case!");

        }
    }

    private void postFilter(List<TreeMap<Integer, Integer>> _mapping, Map<Integer, Integer> BondOrder, int minOrder) {

        for (Map.Entry<Integer, Integer> map : BondOrder.entrySet()) {

            if (map.getValue() > minOrder) {
                removedMap(_mapping, map.getKey());
            }
        }
        FinalMappings final_MAPPINGS = FinalMappings.getInstance();

        final_MAPPINGS.set(new ArrayList<TreeMap<Integer, Integer>>(_mapping));
    }

    private boolean removedMap(List<TreeMap<Integer, Integer>> maps, Integer Key) {
        boolean flag = false;
        for (TreeMap<Integer, Integer> map : maps) {

            if (map.containsKey(Key)) {

                flag = maps.remove(map);
            }
        }
        return flag;
    }
}
