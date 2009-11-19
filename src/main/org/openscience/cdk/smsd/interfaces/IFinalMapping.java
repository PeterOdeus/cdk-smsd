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
package org.openscience.cdk.smsd.interfaces;

import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

/**
 * @cdk.module smsd
 */
public interface IFinalMapping {

   
    /**
     * 
     * @param mapping List of all MCS mapping between a given
     * reactant and product 
     */
    public void add(TreeMap<Integer, Integer> mapping);

    /**
     * 
     * @param mappings List of all MCS mapping between a given
     * reactant and product 
     */
    public void set(List<TreeMap<Integer, Integer>> mappings);

    /**
     * 
     * @return Iterator of mappings
     */
    public Iterator getIterator();

    /**
     * clear the maping
     */
    public void clear();

    /**
     * 
     * @return get of MCS mapping List
     */
    public List<TreeMap<Integer, Integer>> getFinalMapping();

    /**
     * 
     * @return size of the mapping
     */
    public int getSize();
}
