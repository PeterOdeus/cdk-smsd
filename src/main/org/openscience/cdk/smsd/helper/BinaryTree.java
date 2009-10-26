/*
 * BinaryTree.java
 *
 * Created on January 27, 2007, 9:36 PM
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @Revised by Asad, EBI, Java 1.5 Compatable
 */
package org.openscience.cdk.smsd.helper;

//the second part of the program extents the mapping by the McGregor algorithm in case
//that not all atoms of molecule A and molecule B are mapped by the clique approach
public class BinaryTree {

    /** Creates a new instance of BinaryTree
     * @param value
     */
    public BinaryTree(int value) {
        this.Value = value;
    }
    /**
     * equal is initialized as null
     */
    public BinaryTree equal = null;
    /**
     * not equal is initialized as null
     */
    public BinaryTree not_equal = null;
    private int Value = -1;

    /**
     * 
     * @return get the value of the current node
     */
    public int getValue() {
        return this.Value;
    }
}
