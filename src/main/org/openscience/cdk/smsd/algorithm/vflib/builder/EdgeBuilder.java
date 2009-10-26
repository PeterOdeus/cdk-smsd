/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.algorithm.vflib.builder;

import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IBondMatcher;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IEdge;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.INode;


/**
 *
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class EdgeBuilder implements IEdge {

    private NodeBuilder source;
    private NodeBuilder target;
    private IBondMatcher matcher;

    /**
     * 
     * @param source
     * @param target
     * @param matcher
     */
    public EdgeBuilder(NodeBuilder source, NodeBuilder target, IBondMatcher matcher) {
        this.source = source;
        this.target = target;
        this.matcher = matcher;
    }

    
    /**
     * 
     * @return
     */
    @Override
    public INode getSource() {
        return source;
    }

    
    /**
     *
     * @return
     */
    @Override
    public INode getTarget() {
        return target;
    }

    
    /**
     *
     * @return
     */
    @Override
    public IBondMatcher getBondMatcher() {
        return matcher;
    }
}
