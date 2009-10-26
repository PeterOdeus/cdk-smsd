/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.filters;

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
    private static Vector<Integer> extract_clique_mapping(Vector<Integer> comp_graph_nodes, Vector<Integer> clique_vector_org) {

        Vector<Integer> clique_mapping = new Vector<Integer>();
        Vector<Integer> clique_vector = new Vector<Integer>(clique_vector_org);
        int clique_siz = clique_vector.size();
        int vec_size = comp_graph_nodes.size();
        //System.out.println("VEC  SIZE " + vec_size);
        for (int a = 0; a < clique_siz; a++) {
            for (int b = 0; b < vec_size; b = b + 3) {
                if (clique_vector.elementAt(a) == comp_graph_nodes.elementAt(b + 2)) {
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
    public static Vector<Vector<Integer>> extract_mapping(Vector<Vector<Integer>> _mappings, Vector<Integer> comp_graph_nodes, Vector<Integer> clique_vector_org) {

        try {

            //System.out.println("clique_vector_org Size: " + clique_vector_org);
            //System.out.println("FinalMapping Size in extract_mapping: " + FinalMapping.getInstance().getSize());

            //int size=clique_vector_org.size();
            //GlobalVariableContainer.getInstance().setBestCliqueSize(size);


            Vector<Integer> clique_vector = extract_clique_mapping(comp_graph_nodes, clique_vector_org);
            _mappings.addElement(clique_vector);
            
            

            /*Vector<Integer> temp_vector = new Vector<Integer>();

            int clique_siz = clique_vector.size();
            int vec_size = comp_graph_nodes.size();
            for (int a = 0; a < clique_siz; a++) {
            int _cliqueElementVal = clique_vector.get(a);
            for (int b = 0; b < vec_size; b = b + 3) {
            int val = comp_graph_nodes.get(b + 2);
            if (_cliqueElementVal == val) {
            temp_vector.add(comp_graph_nodes.get(b));
            temp_vector.add(comp_graph_nodes.get(b + 1));
            //System.out.println("Added to temp_vector: ");
            }
            }
            }
             _mappings.addElement(temp_vector);*/
            
            
        //System.out.println("In exact Mapping Best Mapping Size set To : " +  temp_vector.size()/2);
        //GlobalVariableContainer.getInstance().setBestMappingSize(temp_vector.size()/2);


        /*System.out.println("temp_vector Inhalt: ");
        for (int c = 0; c < clique_siz; c++) {
        System.out.println(temp_vector.get(c * 2 + 0) + " " + temp_vector.get(c * 2 + 1) + " ");
        }
        System.out.println("");
        System.out.println("");
         */


        //Print MCS_Search solutions on Screen before removing of redundant mappings
			/*int best_MAPPING_size = GlobalVariableContainer.getInstance().getBestMappingSize();
        int count_final_solutions = 1;
        System.out.println("Ausgabe von final_MAPPINGS !!!!!!!!!!!!!: ");
        try {
        Iterator<Vector<Integer>> MAP_iterator = FinalMapping.getInstance().getIterator();
        while (MAP_iterator.hasNext()) {
        Vector<Integer> final_solution = MAP_iterator.next();
        System.out.println("FinalMapping Size: " + best_MAPPING_size + ", Final mapping Nr. " + count_final_solutions++);
        for (int a = 0; a < best_MAPPING_size; a = a + 2) {
        System.out.println(final_solution.get(a) + " " + final_solution.get(a + 1) + "  |  ");
        }
        System.out.println("");
        }
        System.out.println("");
        } catch (Exception ex) {
        ex.printStackTrace();
        }*/


        } catch (Exception e) {
            System.err.println("Error in FinalMapping Vector: " + e.getCause());
            e.printStackTrace();
            System.exit(1);
        }
        return _mappings;
    }
    
    
}
