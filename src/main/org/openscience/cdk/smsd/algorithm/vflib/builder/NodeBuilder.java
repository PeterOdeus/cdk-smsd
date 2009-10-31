/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.algorithm.vflib.builder;


import java.util.ArrayList;
import java.util.List;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IAtomMatcher;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IEdge;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.INode;

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
public class NodeBuilder implements INode {

    private List<INode> neighbors;
    private List<IEdge> edges;
    private IAtomMatcher matcher;

    public NodeBuilder(IAtomMatcher matcher) {
        edges = new ArrayList<IEdge>();
        neighbors = new ArrayList<INode>();
        this.matcher = matcher;
    }

    @Override
    public int countNeighbors() {
        return neighbors.size();
    }

    @Override
    public Iterable<INode> neighbors() {
        return neighbors;
    }

    @Override
    public IAtomMatcher getAtomMatcher() {
        return matcher;
    }

    /**
     * 
     * @return
     */
    @Override
    public List<IEdge> getEdges() {
        return edges;
    }

    /**
     *
     * @param edge
     */
    @Override
    public void addEdge(EdgeBuilder edge) {
        edges.add(edge);
    }

    /**
     *
     * @param node
     */
    @Override
    public void addneighbor(NodeBuilder node) {
        neighbors.add(node);
    }
}

