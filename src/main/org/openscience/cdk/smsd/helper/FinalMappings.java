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

import java.util.ArrayList;
import org.openscience.cdk.smsd.interfaces.IFinalMapping;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

/**
 * @cdk.module smsd
 */
public class FinalMappings implements IFinalMapping {

    private List<TreeMap<Integer, Integer>> _mappings;
    private static FinalMappings INSTANCE = null;

    protected FinalMappings() {
        _mappings = new ArrayList<TreeMap<Integer, Integer>>();
    }

    synchronized public static FinalMappings getInstance() {
        if (INSTANCE == null) {

            INSTANCE = new FinalMappings();

        }
        return INSTANCE;
    }

    @Override
    synchronized public void add(TreeMap<Integer, Integer> v) {
        _mappings.add(v);
    }

    /**
     * 
     * @param v
     */
    @Override
    synchronized public final void set(List<TreeMap<Integer, Integer>> v) {

        _mappings.clear();
        for (TreeMap<Integer, Integer> M : v) {
            _mappings.add(M);
        }
    }

    @Override
    synchronized public Iterator<TreeMap<Integer, Integer>> getIterator() {
        Iterator<TreeMap<Integer, Integer>> it = _mappings.iterator();
        return it;
    }

    @Override
    synchronized public void clear() {
        _mappings.clear();
    }

    @Override
    synchronized public List<TreeMap<Integer, Integer>> getFinalMapping() {
        return _mappings;
    }

    @Override
    synchronized public int getSize() {
        return _mappings.size();
    }
}
