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
package org.openscience.cdk.smsd.filters;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.smsd.helper.FinalMappings;

/**
 * @cdk.module smsd
 */
public class PostFilter {

    /** 
     * 
     * Creates a new instance of PostFilter
     * 
     * @param _mappings 
     */
    public static void filter(List<List<Integer>> _mappings) {
        FinalMappings final_MAPPINGS = FinalMappings.getInstance();
        if (_mappings != null && !_mappings.isEmpty()) {
            final_MAPPINGS.set(removeRedundantMapping(_mappings));
            _mappings.clear();
        } else {
            final_MAPPINGS.set(new ArrayList<TreeMap<Integer, Integer>>());
        }
    }

    private static boolean hasMap(TreeMap<Integer, Integer> newMap, List<TreeMap<Integer, Integer>> nonRedundantMapping) {
        boolean flag = false;

        for (Map<Integer, Integer> map : nonRedundantMapping) {
            if (map.equals(newMap)) {
                flag = true;
                break;
            } else {
                flag = false;
            }
        }
        return flag;
    }

    /**
     *
     * @param mapping_org
     * @return
     */
    public static List<TreeMap<Integer, Integer>> removeRedundantMapping(List<List<Integer>> mapping_org) {

        List<TreeMap<Integer, Integer>> nonRedundantMapping = new ArrayList<TreeMap<Integer, Integer>>();
        for (List<Integer> M : mapping_org) {
            TreeMap<Integer, Integer> newMap = getMappingMapFromList(M);
            if (!hasMap(newMap, nonRedundantMapping)) {
                nonRedundantMapping.add(new TreeMap<Integer, Integer>(newMap));
            }

        }

//        System.out.println("nonRedundantMapping Solutions " + nonRedundantMapping);
        return nonRedundantMapping;
    }

    private static TreeMap<Integer, Integer> getMappingMapFromList(List<Integer> map) {
        TreeMap<Integer, Integer> newMap = new TreeMap<Integer, Integer>();
        for (int index = 0; index < map.size(); index += 2) {
            newMap.put(map.get(index), map.get(index + 1));
        }
        return newMap;
    }
}
