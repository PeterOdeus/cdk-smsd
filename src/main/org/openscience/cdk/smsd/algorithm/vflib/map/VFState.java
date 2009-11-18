/*
 * MX Cheminformatics Tools for Java
 *
 * Copyright (c) 2007-2009 Metamolecular, LLC
 *
 * http://metamolecular.com
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.openscience.cdk.smsd.algorithm.vflib.map;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IEdge;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.INode;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IQuery;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IState;
import org.openscience.cdk.smsd.algorithm.vflib.validator.VFMatch;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @author Richard L. Apodaca <rapodaca at metamolecular.com>
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk> (modified the orignal code)
 * 
 * @cdk.module smsd
 */
public class VFState implements IState {

    private List<VFMatch> candidates;
    private IQuery query;
    private IAtomContainer target;
    private List<INode> queryPath;
    private List<IAtom> targetPath;
    private Map<INode, IAtom> map;

    public VFState(IQuery query, IAtomContainer target) {
        this.map = new HashMap<INode, IAtom>();
        this.queryPath = new ArrayList<INode>();
        this.targetPath = new ArrayList<IAtom>();

        this.query = query;
        this.target = target;
        candidates = new ArrayList<VFMatch>();

        loadRootCandidates();
    }

    private VFState(VFState state, VFMatch match) {
        candidates = new ArrayList<VFMatch>();
        this.queryPath = new ArrayList<INode>(state.queryPath);
        this.targetPath = new ArrayList<IAtom>(state.targetPath);
        this.map = state.map;
        this.query = state.query;
        this.target = state.target;

        map.put(match.getQueryNode(), match.getTargetAtom());
        queryPath.add(match.getQueryNode());
        targetPath.add(match.getTargetAtom());

        loadCandidates(match);
    }

    @Override
    public void backTrack() {
        if (queryPath.isEmpty() || isGoal()) {
            map.clear();

            return;
        }

        if (isHeadMapped()) {
            return;
        }

        map.clear();

        for (int i = 0; i < queryPath.size() - 1; i++) {
            map.put(queryPath.get(i), targetPath.get(i));
        }

    }

    @Override
    public Map<INode, IAtom> getMap() {
        return new HashMap<INode, IAtom>(map);
    }

    @Override
    public boolean hasNextCandidate() {
        return !candidates.isEmpty();
    }

    @Override
    public boolean isDead() {
        return query.countNodes() > target.getAtomCount();
    }

    @Override
    public boolean isGoal() {
        return map.size() == query.countNodes();
    }

    @Override
    public boolean isMatchFeasible(VFMatch match) {
        if (map.containsKey(match.getQueryNode()) || map.containsValue(match.getTargetAtom())) {
            return false;
        }

        if (!matchAtoms(match)) {
            return false;
        }

        if (!matchBonds(match)) {
            return false;
        }

        return true;
    }

    @Override
    public VFMatch nextCandidate() {
        return candidates.remove(candidates.size() - 1);
    }

    @Override
    public IState nextState(VFMatch match) {
        return new VFState(this, match);
    }

    private void loadRootCandidates() {
        for (int i = 0; i < query.countNodes(); i++) {
            for (int j = 0; j < target.getAtomCount(); j++) {
                VFMatch match = new VFMatch(query.getNode(i), target.getAtom(j));
                candidates.add(match);
            }
        }


    }

    private void loadCandidates(VFMatch lastMatch) {
        IAtom a = lastMatch.getTargetAtom();
        List<IAtom> targetNeighbors = target.getConnectedAtomsList(a);

        for (INode q : lastMatch.getQueryNode().neighbors()) {
            for (IAtom t : targetNeighbors) {
                VFMatch match = new VFMatch(q, t);

                if (candidateFeasible(match)) {
                    candidates.add(match);
                }
            }
        }


    }

    private boolean candidateFeasible(VFMatch candidate) {
        for (INode queryAtom : map.keySet()) {
            if (queryAtom.equals(candidate.getQueryNode()) ||
                    map.get(queryAtom).equals(candidate.getTargetAtom())) {
                return false;
            }
        }

        return true;
    }

    private boolean matchAtoms(VFMatch match) {
        IAtom a = match.getTargetAtom();
        if (match.getQueryNode().countNeighbors() > target.getConnectedAtomsCount(a)) {
            return false;
        }


        return match.getQueryNode().getAtomMatcher().matches(a);
    }

    private boolean matchBonds(VFMatch match) {
        if (queryPath.isEmpty()) {
            return true;
        }

        if (!matchBondsToHead(match)) {
            return false;
        }

        for (int i = 0; i < queryPath.size() - 1; i++) {
            IEdge queryBond = query.getEdge(queryPath.get(i), match.getQueryNode());
            IBond targetBond = target.getBond(targetPath.get(i), match.getTargetAtom());

            if (queryBond == null) {
                continue;
            }

            if (targetBond == null) {
                return false;
            }

            if (!matchBond(queryBond, targetBond)) {
                return false;
            }
        }

        return true;
    }

    private boolean matchBondsToHead(VFMatch match) {
        INode queryHead = getQueryPathHead();
        IAtom targetHead = getTargetPathHead();

        IEdge queryBond = query.getEdge(queryHead, match.getQueryNode());
        IBond targetBond = target.getBond(targetHead, match.getTargetAtom());

        if (queryBond == null || targetBond == null) {
            return false;
        }

        return matchBond(queryBond, targetBond);
    }

    private INode getQueryPathHead() {
        return queryPath.get(queryPath.size() - 1);
    }

    private IAtom getTargetPathHead() {
        return targetPath.get(targetPath.size() - 1);
    }

    private boolean matchBond(IEdge edge, IBond targetBond) {


        return edge.getBondMatcher().matches(targetBond);

    }

    private boolean isHeadMapped() {
        INode head = queryPath.get(queryPath.size() - 1);
        for (INode neighbor : head.neighbors()) {
            if (!map.containsKey(neighbor)) {
                return false;
            }
        }

        return true;
    }
}
