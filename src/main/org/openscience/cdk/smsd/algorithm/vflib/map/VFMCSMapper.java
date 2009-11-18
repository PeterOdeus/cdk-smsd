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
package org.openscience.cdk.smsd.algorithm.vflib.map;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IMapper;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.INode;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IQuery;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IState;
import org.openscience.cdk.smsd.algorithm.vflib.query.TemplateCompiler;
import org.openscience.cdk.smsd.algorithm.vflib.validator.VFMatch;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @cdk.module smsd
 */
public class VFMCSMapper implements IMapper {

    private IQuery query;
    private List<Map<INode, IAtom>> maps;
    private int StoredMCSSize = -1;

    /**
     *
     * @param query
     */
    public VFMCSMapper(IQuery query) {
        this.query = query;
        this.maps = new Vector<Map<INode, IAtom>>();
    }

    /**
     *
     * @param molecule
     */
    public VFMCSMapper(IAtomContainer molecule) {
        this.query = TemplateCompiler.compile(molecule);
        this.maps = new Vector<Map<INode, IAtom>>();
    }

    /**
     *
     * @param target
     * @return
     */
    @Override
    public boolean hasMap(IAtomContainer target) {
        VFMCSState state = new VFMCSState(query, target);
        maps.clear();

        return mapFirst(state);
    }

    @Override
    public List<Map<INode, IAtom>> getMaps(IAtomContainer target) {
        VFMCSState state = new VFMCSState(query, target);
        maps.clear();
        mapAll(state);
        return new Vector<Map<INode, IAtom>>(maps);
    }

    /**
     *
     * @param target
     * @return
     */
    @Override
    public Map<INode, IAtom> getFirstMap(IAtomContainer target) {
        VFMCSState state = new VFMCSState(query, target);
        maps.clear();
        mapFirst(state);
        return maps.isEmpty() ? new HashMap<INode, IAtom>() : maps.get(0);
    }

    /**
     *
     * @param target
     * @return
     */
    @Override
    public int countMaps(IAtomContainer target) {
        VFMCSState state = new VFMCSState(query, target);
        maps.clear();
        mapAll(state);
        return maps.size();
    }

    private void mapAll(IState state) {
        if (state.isDead()) {
            return;
        }

        if (state.isGoal()) {

            Map<INode, IAtom> map = state.getMap();
            if (!hasMap(map) && isMCS(map)) {
                maps.add(state.getMap());
            }
            return;
        } else {

            Map<INode, IAtom> map = state.getMap();
            if (!map.isEmpty() && !hasMap(map) && isMCS(map)) {
                maps.add(state.getMap());
            }
        }


        while (state.hasNextCandidate()) {
            VFMatch candidate = state.nextCandidate();

            if (state.isMatchFeasible(candidate)) {
                IState nextState = state.nextState(candidate);
                mapAll(nextState);
                nextState.backTrack();

            }
        }
    }

    private boolean mapFirst(IState state) {
        if (state.isDead()) {
            return false;
        }

        if (state.isGoal()) {
            maps.add(state.getMap());
            return true;
        } else {

            Map<INode, IAtom> map = state.getMap();
            if (!map.isEmpty() && isMCS(map) && !hasMap(map)) {//!hasSubGraph(map)) {
                maps.add(state.getMap());
            }
        }

        boolean found = false;

        while (!found && state.hasNextCandidate()) {
            VFMatch candidate = state.nextCandidate();
            if (state.isMatchFeasible(candidate)) {
                IState nextState = state.nextState(candidate);
                found = mapFirst(nextState);
                nextState.backTrack();

            }
        }

        return found;
    }

//    //Method added by Asad
//    private boolean hasSubGraph(Map<INode, IAtom> map) {
//
//        for (Map<INode, IAtom> storedMappings : maps) {
//            boolean MatchFlag = true;
//            for (Map.Entry<INode, IAtom> mapping : storedMappings.entrySet()) {
//                if (!map.containsKey(mapping.getKey()) && !map.containsValue(mapping.getValue())) {
//                    MatchFlag = false;
//                    break;
//                }
//            }
//
//            if (MatchFlag) {
//
//                return true;
//            }
//
//        }
//
//        return false;
//    }
    //Method added by Asad
    private boolean isMCS(Map<INode, IAtom> map) {

        boolean flag = true;
        int mapSize = map.size();

        if (!maps.isEmpty() && StoredMCSSize > mapSize) {

            flag = false;
        }
        //Comment this if to get all the subgraphs
        if (mapSize > StoredMCSSize) {

            StoredMCSSize = mapSize;
            maps.clear();
        }

        return flag;
    }

    private boolean hasMap(Map<INode, IAtom> map) {
        for (Map<INode, IAtom> test : maps) {
            if (test.equals(map)) {
                return true;
            }
        }

        return false;
    }
}
