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
    private List<TreeMap<Integer, Integer>> _mappings = null;
    private Map<Integer, Integer> connectedBondOrder = null;

    /**
     *
     * @param source
     * @param target
     * @param removeHydrogen
     */
    protected void getOverLaps(IAtomContainer source, IAtomContainer target, boolean removeHydrogen) {


        _mappings = new ArrayList<TreeMap<Integer, Integer>>();
        connectedBondOrder = new TreeMap<Integer, Integer>();
        this.source = source;
        this.target = target;

        int minOrder = 9999;

        if (source.getAtomCount() == 1) {
            minOrder = setSourceSingleAtomMap(removeHydrogen);
        }
        if (target.getAtomCount() == 1) {
            minOrder = setTargetSingleAtomMap(removeHydrogen);
        }
        postFilter(minOrder);

        _mappings.clear();
    }

    private int setSourceSingleAtomMap(boolean removeHydrogen) {
        int counter = 0;
        int minOrder = 9999;
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

                    connectedBondOrder.put(counter, totalOrder);
                    _mappings.add(counter++, mapAtoms);
                }

                //System.out.println("Hello in Single getOverlaps Mapping Size: " + mapAtoms.size());
            }
        } else {
            System.err.println("Skippping Hydrogen mapping or This is not a single mapping case!");
        }

        return minOrder;
    }

    private int setTargetSingleAtomMap(boolean removeHydrogen) {
        int counter = 0;
        int minOrder = 9999;
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

                    connectedBondOrder.put(counter, totalOrder);
                    _mappings.add(counter++, mapAtoms);
                }

//                System.out.println("Hello in Single getOverlaps Mapping Size: " + mapAtoms.size());
            }

        } else {
            System.err.println("Skippping Hydrogen mapping or This is not a single mapping case!");
        }
        return minOrder;
    }

    private void postFilter(int minOrder) {
        for (Map.Entry<Integer, Integer> map : connectedBondOrder.entrySet()) {
            if (map.getValue() > minOrder) {
                removedMap(map.getKey());
            }
        }
        FinalMappings final_MAPPINGS = FinalMappings.getInstance();
        final_MAPPINGS.set(new ArrayList<TreeMap<Integer, Integer>>(_mappings));
    }

    private void removedMap(Integer Key) {
        List<TreeMap<Integer, Integer>> nonRedundantMap = new ArrayList<TreeMap<Integer, Integer>>(_mappings);
        for (TreeMap<Integer, Integer> map : nonRedundantMap) {
            if (map.containsKey(Key)) {
                _mappings.remove(map);
            }
        }
    }
}
