/*
 * PostFilter.java
 *
 * Created on 06 March 2007, 14:52
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package org.openscience.cdk.smsd.filters;

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
    public static void filter(Vector<Vector<Integer>> _mappings) {
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

    private static boolean hasMap(TreeMap<Integer, Integer> newMap, Vector<TreeMap<Integer, Integer>> nonRedundantMapping) {
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

    private static Vector<TreeMap<Integer, Integer>> removeRedundantMapping(Vector<Vector<Integer>> mapping_org) {

//        System.out.println("Total Number of Solutions " + mapping_org.size());
        Vector<TreeMap<Integer, Integer>> nonRedundantMapping = new Vector<TreeMap<Integer, Integer>>();
        for (Vector<Integer> M : mapping_org) {
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
