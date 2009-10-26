/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.filters;

import java.util.List;
import java.util.Vector;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
public class ExactMapping {

    //extract atom mapping from the clique vector and store it in vector clique_MAPPING
    /**
     * 
     * @param comp_graph_nodes
     * @param clique_vector_org
     * @return 
     */
    private static List<Integer> extract_clique_mapping(List<Integer> comp_graph_nodes, List<Integer> clique_vector_org) {

        List<Integer> clique_mapping = new Vector<Integer>();
        List<Integer> clique_vector = new Vector<Integer>(clique_vector_org);
        int clique_siz = clique_vector.size();
        int vec_size = comp_graph_nodes.size();
        //System.out.println("VEC  SIZE " + vec_size);
        for (int a = 0; a < clique_siz; a++) {
            for (int b = 0; b < vec_size; b = b + 3) {
                if (clique_vector.get(a) == comp_graph_nodes.get(b + 2)) {
                    clique_mapping.add(comp_graph_nodes.get(b));
                    clique_mapping.add(comp_graph_nodes.get(b + 1));
                }
            }
        }

        return clique_mapping;
    }

    //extract atom mapping from the clique vector and print it on the screen
    /**
     * 
     * @param _mappings
     * @param comp_graph_nodes
     * @param clique_vector_org
     * @return
     */
    public static List<List<Integer>> extract_mapping(List<List<Integer>> _mappings, List<Integer> comp_graph_nodes, List<Integer> clique_vector_org) {

        try {

            //System.out.println("clique_vector_org Size: " + clique_vector_org);
            //System.out.println("FinalMapping Size in extract_mapping: " + FinalMapping.getInstance().getSize());

            //int size=clique_vector_org.size();
            //GlobalVariableContainer.getInstance().setBestCliqueSize(size);


            List<Integer> clique_vector = extract_clique_mapping(comp_graph_nodes, clique_vector_org);
            _mappings.add(clique_vector);
            
        } catch (Exception e) {
            System.err.println("Error in FinalMapping Vector: " + e.getCause());
            e.printStackTrace();
            System.exit(1);
        }
        return _mappings;
    }
    
    
}
