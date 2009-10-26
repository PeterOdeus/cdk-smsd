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
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
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

