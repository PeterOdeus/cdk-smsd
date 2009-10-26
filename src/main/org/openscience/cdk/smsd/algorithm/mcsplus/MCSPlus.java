/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.algorithm.mcsplus;


import java.io.IOException;
import java.util.Stack;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.smsd.algorithm.mcgregor.McGregor;
import org.openscience.cdk.smsd.core.tools.EBIException;
import org.openscience.cdk.smsd.core.tools.EBITimeManager;
import org.openscience.cdk.smsd.filters.ExactMapping;
import org.openscience.cdk.smsd.global.TimeOut;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK  
 * @contact asad@ebi.ac.uk
 */
public class MCSPlus {

    public static double timeout = TimeOut.getInstance().getTimeOut();
    private static EBITimeManager TM;
    private static boolean timeoutFlag = false;

    /**
     * 
     * @param ac1
     * @param ac2
     * @param removeHydrogen
     * @return
     * @throws EBIException
     */
    public Vector<Vector<Integer>> getOverlaps(IAtomContainer ac1, IAtomContainer ac2, boolean removeHydrogen) throws EBIException {
        Stack<Vector<Integer>> Max_Cliques_Set = null;
        Vector<Vector<Integer>> _mappings = new Vector<Vector<Integer>>();


        try {

            GenerateCompatibilityGraph gcg = new GenerateCompatibilityGraph(ac1, ac2, removeHydrogen);


//            System.out.println("Best Clique Size: " + best_clique_size);


            Vector<Integer> comp_graph_nodes = gcg.get_comp_graph_nodes();

            Vector<Integer> C_edges = gcg.get_C_egdes();
            Vector<Integer> D_edges = gcg.get_D_egdes();

            if (D_edges.size() > 99999 && C_edges.size() > 2000) {
                System.err.println("D-edges Size " + D_edges.size() + " > : " + 99999);
                return null;
            }
//            System.out.println("C-edges: " + C_edges);
//            System.out.println("D-edges: " + D_edges);
//            System.out.println("Comp Graph: " + comp_graph_nodes);
//            System.out.println("C_edges: " + C_edges.size());
//            System.out.println("D_edges: " + D_edges.size());
            TM = new EBITimeManager();

            // build the CDKRGraph corresponding to this problem

            if (timeout > -1 && TM.getElapsedTimeInMinutes() > timeout) {
                //System.out.println("|Hello|");
                timeoutFlag = true;
                throw new EBIException("Timeout exceeded in getOverlaps");
            }
            BKKCKCF init = new BKKCKCF(comp_graph_nodes, C_edges, D_edges);
//            BronKerboschKochCliqueFinder init = new BronKerboschKochCliqueFinder(comp_graph_nodes, C_edges, D_edges);


            Max_Cliques_Set = init.getMaxCliqueSet(); // build the CDKRGraph corresponding to this problem

            if (timeout > -1 && TM.getElapsedTimeInMinutes() > timeout) {
                //System.out.println("|Hello|");
                timeoutFlag = true;
                throw new EBIException("Timeout exceeded in getOverlaps");
            }

            int solution_size_pass1 = Max_Cliques_Set.size();
            int max_map_size_allowed = init.getBestCliqueSize();
//            System.err.println("best_clique_size: vorher " + init.getBestCliqueSize());
//            System.err.println("Max_Cliques_Set: " + solution_size_pass1);
//            System.out.println("Maximum Mapping size allowed " + max_map_size_allowed);
//            System.err.println("Cliques: " + init.getMaxCliqueSet());

            //Clear all the compatibility graph content
            gcg.Clear();
//            System.gc();

            int clique_number = 1;
            while (!Max_Cliques_Set.empty()) {

                // build the CDKRGraph corresponding to this problem

                if (timeout > -1 && TM.getElapsedTimeInMinutes() > timeout) {
                    //System.out.println("|Hello|");
                    timeoutFlag = true;
                    throw new EBIException("Timeout exceeded in getOverlaps");
                }


//                System.out.println("Clique number: " + clique_number);
                Vector<Integer> clique_vector = Max_Cliques_Set.peek();

                int clique_size = clique_vector.size();

                //Uncommented by Asad for Debugging
                /*System.out.println(" Degug Session in Clique extension ");
                System.out.println("Test Loop Asad");
                for (int a = 0; a < clique_size; a++) {
                System.out.print(clique_vector.get(a) + " ");
                }
                System.out.println(" ");*/

                //Is the number of mappings smaller than the number of atomsMCS of molecule A and B?
                //In this case the clique is given to the McGregorMCSPlus algorithm
                if (clique_size < ac1.getAtomCount() && clique_size < ac2.getAtomCount()) {
//
//                    System.out.println("clique_size: " + clique_size + " atom_number1: " + ac1.getAtomCount() + " atom_number2: " + ac2.getAtomCount());
//                    System.out.println("McGregorMCSPlus");
//                    System.out.println("comp_graph_nodes " + comp_graph_nodes);
//                    System.out.println("clique_vector " + clique_vector);

                    McGregor mgit = new McGregor(ac1, ac2, _mappings);

                    //McGregorOptimized mgit = new McGregorOptimized(comp_graph_nodes, _mappings);
                    mgit.McGregor_IterationStart(mgit.getMCSSize(), clique_vector, comp_graph_nodes); //Start McGregorMCSPlus search

                    _mappings = mgit.getMappings();

                    mgit = null;

//                    System.out.println("_mappings " + _mappings);
//                    System.out.println("_mappings : " + _mappings.size());
//                    System.out.println("McGregorMCSPlus Over");
//                    System.out.println();


                } else {

                    //Vector<Integer> clique_mapping = new Vector<Integer>();
                    _mappings = ExactMapping.extract_mapping(_mappings, comp_graph_nodes, clique_vector); //erfolgt KEINE McGregorMCSPlus-Suche->speichere L�ungen in final_MAPPINGS

                }
                Max_Cliques_Set.pop();
                clique_number++;


            }
//
//            int best_MAPPING_size = GlobalVariableContainer.getInstance().getMappingSize();
//            System.out.println("best_clique_size:  " + best_clique_size);
//            System.out.println("best_MAPPING_size: " + best_MAPPING_size);
//            System.out.println("**************************************************");

//            System.out.println("Entering INTO PostFilter Stage");
//            System.out.println("Number of Solutions: " + _mappings.size());

        } catch (IOException ex) {
            Logger.getLogger(MCSPlus.class.getName()).log(Level.SEVERE, null, ex);
        }

        return _mappings;
    }

    /**
     *
     * @return
     */
    public static boolean getTimeOutFlag() {

        boolean flag = false;

        if (timeoutFlag) {
            flag = true;
        }
        System.out.println("Flag: " + flag);
        return flag;
    }
}
