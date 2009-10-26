/*
 * PostFilter.java
 *
 * Created on 06 March 2007, 14:52
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
package org.openscience.cdk.smsd.filters;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.smsd.helper.FinalMappings;


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
