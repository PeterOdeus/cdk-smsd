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
package org.openscience.cdk.smsd.algorithm.vflib.builder;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IEdge;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.INode;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IQuery;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @cdk.module smsd
 */
public class VFQueryBuilder implements IQuery {

    private List<INode> nodes;
    private List<IEdge> edges;
    private Map<INode, IAtom> NodesBonds;

    public VFQueryBuilder() {
        nodes = new ArrayList<INode>();
        edges = new ArrayList<IEdge>();
        NodesBonds = new HashMap<INode, IAtom>();
    }

    /**
     *
     * @return
     */
    @Override
    public Iterable<IEdge> edges() {
        return edges;
    }

    /**
     *
     * @return
     */
    @Override
    public Iterable<INode> nodes() {
        return nodes;
    }

    /**
     *
     * @param index
     * @return
     */
    @Override
    public INode getNode(int index) {
        return nodes.get(index);
    }

    /**
     *
     * @param atom 
     * @return
     */
    public INode getNode(IAtom atom) {



        for (Map.Entry<INode, IAtom> v : NodesBonds.entrySet()) {
            if (v.getValue().equals(atom)) {
                return v.getKey();
            }
        }

        return null;
    }

    /**
     *
     * @param index
     * @return
     */
    @Override
    public IEdge getEdge(int index) {
        return edges.get(index);
    }

    /**
     *
     * @param source
     * @param target
     * @return
     */
    @Override
    public IEdge getEdge(INode source, INode target) {
        if (source == target) {
            return null;
        }

        NodeBuilder sourceImpl = (NodeBuilder) source;

        for (IEdge edge : sourceImpl.getEdges()) {
            if (edge.getSource() == target || edge.getTarget() == target) {
                return edge;
            }
        }

        return null;
    }

   
    /**
     *
     * @param matcher
     * @param atom
     * @return
     */
    public INode addNode(IQueryAtom matcher, IAtom atom) {
        NodeBuilder node = new NodeBuilder(matcher);
        nodes.add(node);
        NodesBonds.put(node, atom);
        return node;
    }

     /**
     *
     * @param node
     * @return
     */
    @Override
    public IAtom getAtom(INode node) {

        return NodesBonds.get(node);
    }

    @Override
    public int countNodes() {
        return nodes.size();
    }

    /**
     *
     * @return
     */
    @Override
    public int countEdges() {
        return edges.size();
    }

    /**
     * 
     * @param source
     * @param target
     * @param matcher
     * @return
     */
    public IEdge connect(INode source, INode target, IQueryBond matcher) {
        NodeBuilder sourceImpl = (NodeBuilder) source;
        NodeBuilder targetImpl = (NodeBuilder) target;
        EdgeBuilder edge = new EdgeBuilder(sourceImpl, targetImpl, matcher);

        sourceImpl.addneighbor(targetImpl);
        targetImpl.addneighbor(sourceImpl);

        sourceImpl.addEdge(edge);
        targetImpl.addEdge(edge);

        edges.add(edge);
        return edge;
    }
}
 