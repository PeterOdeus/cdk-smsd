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

import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
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

//        System.out.println("Total Number of Solutions " + _mappings.size());

        if (_mappings != null && !_mappings.isEmpty()) {

            final_MAPPINGS.set(removeRedundantMapping(_mappings));
            _mappings.clear();
        } else {

            final_MAPPINGS.set(new Vector<TreeMap<Integer, Integer>>());
        }
//        System.out.println("Total Number of Solutions " + final_MAPPINGS.getFinalMapping().size());
        return;
    }

    private static boolean hasMap(TreeMap<Integer, Integer> newMap, List<TreeMap<Integer, Integer>> nonRedundantMapping) {
        boolean flag = false;

        for (Map<Integer, Integer> map : nonRedundantMapping) {
            if (map.entrySet().equals(newMap)) {
                flag = true;
            } else {
                flag = false;
                break;
            }
        }



        return flag;
    }

    private static List<TreeMap<Integer, Integer>> removeRedundantMapping(List<List<Integer>> mapping_org) {

//        System.out.println("Total Number of Solutions " + mapping_org.size());
        List<TreeMap<Integer, Integer>> nonRedundantMapping = new Vector<TreeMap<Integer, Integer>>();
        for (List<Integer> M : mapping_org) {
//            System.out.println("Solutions " + M);
            TreeMap<Integer, Integer> newMap = new TreeMap<Integer, Integer>();
            for (int index = 0; index < M.size(); index += 2) {
                newMap.put(M.get(index), M.get(index + 1));
            }

            if (!hasMap(newMap, nonRedundantMapping)) {

                nonRedundantMapping.add(new TreeMap<Integer, Integer>(newMap));
            }

        }

//        System.out.println("nonRedundantMapping Solutions " + nonRedundantMapping);
        return nonRedundantMapping;
    }
}
