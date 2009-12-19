
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
package org.openscience.cdk.smsd.helper;


/**
 * @cdk.module smsd
 */
public class BinaryTree {

    /**
     * Creates a new instance of BinaryTree
     * the second part of the program extents the mapping by the McGregor algorithm in case
     * that not all atoms of molecule A and molecule B are mapped by the clique approach
     * @param value
     */
    public BinaryTree(int value) {
        this.value = value;
    }
    /**
     * equal is initialized as null
     */
    public BinaryTree equal = null;
    /**
     * not equal is initialized as null
     */
    public BinaryTree notEqual = null;
    private int value = -1;

    /**
     * 
     * @return get the value of the current node
     */
    public int getValue() {
        return this.value;
    }
}
