/*
 * IFinalMapping.java
 *
 * Created on Sep 12, 2007, 8:58:45 PM
 *
 */
package org.openscience.cdk.smsd.interfaces;

import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

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
public interface IFinalMapping {

   
    /**
     * 
     * @param v List of all MCS mapping between a given
     * reactant and product 
     * @throws java.lang.Exception
     */
    public void add(TreeMap<Integer, Integer> v) throws Exception;

    /**
     * 
     * @param v List of all MCS mapping between a given
     * reactant and product 
     */
    public void set(List<TreeMap<Integer, Integer>> v);

    /**
     * 
     * @return Iterator of mappings
     */
    public Iterator getIterator();

    /**
     * clear the maping
     */
    public void Clear();

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
