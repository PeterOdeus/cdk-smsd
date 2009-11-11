/*
 * BronKerboschCazalsKarandeKochCliqueFinder.java
 *
 */

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
package org.openscience.cdk.smsd.algorithm.mcsplus;

import java.util.List;
import java.util.Stack;
import java.util.Vector;

/**
 * This class implements Bron-Kerbosch clique detection algorithm as it is
 * described in [F. Cazals, R. Karande: An Algorithm for reporting maximal c-cliques;
 * T.Comp. Sc. (2005); vol 349; pp.
 * 484-490]
 * 
 * @cdk.module smsd
 */
public class BKKCKCF {

    private Vector<Vector<Integer>> Max_Cliques_Set;
    /***********************************************************************/
    private List<Integer> C_edges;
    private List<Integer> D_edges;
    private int best_clique_size;
    private List<Integer> comp_graph_nodes;
    private double D_edge_Iteration_size = 0;
    private double C_edge_Iteration_size = 0;

    /**
     * Creates a new instance of BronKerboschKochCliqueFinder
     * @param comp_graph_nodes_org
     * @param C_edges_org C-Edges set of allowed edges
     * @param D_edges_org D-Edges set of prohibited edges
     */
    public BKKCKCF(List<Integer> comp_graph_nodes_org, List<Integer> C_edges_org, List<Integer> D_edges_org) {

        this.comp_graph_nodes = comp_graph_nodes_org;
        this.C_edges = C_edges_org;
        this.D_edges = D_edges_org;
        best_clique_size = 0;

//        System.out.println("C Edges Size: " + C_edges.size());
//        System.out.println("D Edges Size:" + D_edges.size());



        //Orignal assignment as per paper
        D_edge_Iteration_size = D_edges.size()/2;
        //Heuristic introduced by Asad

        //Orignal assignment as per paper
        C_edge_Iteration_size = C_edges.size()/2;
        //Heuristic introduced by Asad


        boolean d_edgeFlag = false;

        if (D_edges.size() > 10000000 && D_edges.size() > C_edges.size() && C_edges.size() > 100000) {
            D_edge_Iteration_size = (float) D_edges.size() * 0.000001;
            d_edgeFlag = true;

        } else if (D_edges.size() > 10000000 && D_edges.size() > C_edges.size() && 5000 < C_edges.size()) {
            D_edge_Iteration_size = (float) D_edges.size() * 0.001;
            d_edgeFlag = true;

        }
//        else if (D_edges.size() > 5000000 && D_edges.size() > C_edges.size()) {
//            D_edge_Iteration_size = (float) D_edges.size() * 0.0001;
//            d_edgeFlag = true;
//
//        } else if (D_edges.size() > 100000 && D_edges.size() > C_edges.size()) {
//            D_edge_Iteration_size = (float) D_edges.size() * 0.1;
//            d_edgeFlag = true;
//        }

//        } else if (D_edges.size() >= 10000 && 500 >= C_edges.size()) {
//            D_edge_Iteration_size = (float) D_edges.size() * 0.1;
//            d_edgeFlag = true;
//
//        }
//
//
//
        if (D_edge_Iteration_size < 1 && d_edgeFlag && C_edges.size() <= 5000) {

            D_edge_Iteration_size = 2;
        }

        if (D_edge_Iteration_size < 1 && d_edgeFlag) {

            D_edge_Iteration_size = 1;
        }

//        System.out.println("Optimised D-edges: " + (D_edge_Iteration_size));

        //Initialization Max_Cliques_Set

        Max_Cliques_Set = new Vector<Vector<Integer>>();

        Init_Algorithm();

        /*System.out.println();
        System.out.println("Clique Search Over");
        System.out.println("_______________________________");
         */

    }

    /*
     * Call the wrapper for ENUMERATE_CLIQUES
     *
     */
    private void Init_Algorithm() {


        /********************************************************************/
        /*
         *T: is a set of vertices which have already been used for the
         * initialization of ENUMERATE_CLIQUES
         */
        Vector<Integer> T = new Vector<Integer>(); //Initialize the T Vector;

        /*
         *V: stored all the vertices for the Graph G
         * V[G]
         *nodes of vector comp_graph_nodes are stored in V
         */

        Vector<Integer> V = new Vector<Integer>(); //Initialization of Vector V

        int V_set_size = comp_graph_nodes.size() / 3;

        //System.out.println("Vector V is initialized");
        for (int a = 0; a < V_set_size; a++) {
            V.add(comp_graph_nodes.get(a * 3 + 2));
            //System.out.print("V[" + a + "]: " + comp_graph_nodes.get(a * 3 + 2) + " ");
        }
        //System.out.println();

        V.add(0);
        // System.out.println("Vector V :" + V);


        /*
         * R: set of vertices belonging to the current clique
         */
        Vector<Integer> R = new Vector<Integer>();
        /*
         *P: is a set of vertices which <b>can</b> be added to R, because they are
         * neighbours of vertex u via <i>c-edges</i>
         */
        Stack<Integer> P = new Stack<Integer>();
        /*
         *Q: is a set of vertices which <b>cannot</b> be added to R, because they are
         * neighbours of vertex u via <i>d-edges</i>
         */

        Vector<Integer> Q = new Vector<Integer>();
        /*
         *X: set of vertices which are not allowed to be added
         * to R
         */
        Vector<Integer> X = new Vector<Integer>();


        /*
         *Y: set of vertices which are not allowed to be added
         * to C
         */

        Vector<Integer> Y = new Vector<Integer>();



        /*
         * N[u]: set of neighbours of vertex u in Graph G
         *
         */

        Vector<Integer> N = new Vector<Integer>();

        int b = 0;

        /*
         * Let T be the set of Nodes already been used in the initialization
         *
         */

        T.clear();

        while (V.get(b) != 0) {


            // V[b] is node u, v belogs to V[G]


            //System.out.println();
            //System.out.println("#########################################");
            //System.out.println("Central node " + V.get(b));


            int central_node = V.get(b);


            P.clear();
            Q.clear();
            X.clear();
            R.clear();

            //find the neighbors of the central node from V
            N = find_neighbors(central_node);

            //System.out.println("N-Neigh: " + N.size());



            //gehe Nachbarn durch und ordne sie zu X, P oder Q

            for (int c = 0; c < N.size(); c = c + 2) {
                // N[c] is node v
                // System.out.println("N[" + c  +  "]= " + N.get(c) + " ");
                //Grouping of the neighbors in X,P and Q

                /*
                 * u and v are adjacent via a R-edge
                 */
                int N_at_c = N.get(c);


                //find respective neighbor position in P, which is needed for the deletion from V
                //delete neighbor from set V


                if (N.elementAt(c + 1) == 1) {


                    if (T.contains(N_at_c)) {
                        X.add(N_at_c);
                    } else {
                        P.push(N_at_c);
                    }

                } else if (N.elementAt(c + 1) == 2) {
                    // u and v are adjacent via a Q-edge
                    //System.out.println("u and v are adjacent via a Q-edge: " + N.elementAt(c));

                    if (T.contains(N_at_c)) {
                        Y.add(N_at_c);
                    } else {
                        Q.add(N_at_c);
                    }
                }

                if (V.indexOf(N_at_c) <= b && V.indexOf(N_at_c) > -1) {
                    --b;

                }
                V.removeElement(N_at_c);
                //System.out.println("Elements Removed from V:" + N_at_c);



            }

            P.add(0);
            R.add(central_node);

            /*if (checkRepeatFlag) {
            checkRepeatFlag = false;
            break;
            }*///Added by asad for huristic check

            //System.out.println(" Calling Enumerate_Cliques, Vector Size in CliquesGenerator: ");

            //System.out.println("C: " + R.size() + ", P: " + P.size() + ", D: " + Q.size() + ", S: " + X.size());

            Enumerate_Cliques(R, P, Q, X, Y);
            //Enumerate_Cliques(R, P, Q, X);
            T.add(central_node);

            b++;

            //System.out.println("HELLO, B: " + b + " V Size: " + V.size());
            //checkRepeatFlag = false;
        }
        //System.out.println("Max_Cliques_Set: " + Max_Cliques_Set);

    }

    private int Enumerate_Cliques(Vector<Integer> R, Stack<Integer> P, Vector<Integer> Q, Vector<Integer> X, Vector<Integer> Y) {
        //private int Enumerate_Cliques(Vector<Integer> R, Stack<Integer> P, Vector<Integer> Q, Vector<Integer> X) {


        Vector<Integer> N = new Vector<Integer>(); ////Initialization Vector N
        Stack<Integer> ut_set = new Stack<Integer>();//Defined as P' in the paper


        for (Integer I : P) {
            ut_set.add(I);
        }

        if (P.size() == 1) {
            //if (P.size() == 0) {
            if (X.size() == 0) {

                //store best solutions in stack Max_Cliques_Set
                int clique_size = R.size();



                if (clique_size >= best_clique_size) {
                    if (clique_size > best_clique_size) {

                        Max_Cliques_Set.clear();
                        best_clique_size = clique_size;

                    }
                    if (clique_size == best_clique_size) {
                        //System.out.println("R-Clique " + R);
                        Max_Cliques_Set.addElement(R);

//                       System.out.println("Best Cliques Size: " + best_clique_size + " " + clique_size);
                    }

                }

                return 0;
            }
        }
        //Added by Asad
        /*if(best_clique_size==GlobalVariableContainer.getInstance().getReactantAtomSize()||
        best_clique_size==GlobalVariableContainer.getInstance().getProductAtomSize()){
        System.out.println("Found");
        return 0;
        }*/

        int a = 0;

        while (ut_set.elementAt(a) != 0) {

            // P[a] is node ut
            //  System.out.println("Central recursion node " + ut_set.elementAt(a));
            int ui = ut_set.get(a);
            //remove ut_set[a] from P
            //find position of ut_set node in P

            //int ui = ut_set.peek();

            P.removeElement(ui);



            Vector<Integer> R_copy = new Vector<Integer>(R);
            Stack<Integer> P_copy = new Stack<Integer>();
            Stack<Integer> Q_copy = new Stack<Integer>();
            Vector<Integer> X_copy = new Vector<Integer>(X);
            Vector<Integer> Y_copy = new Vector<Integer>(Y);
            //Vector<Integer> Y_copy = new Vector<Integer>();

            N.clear();


            for (Integer obj : P) {
                P_copy.add(obj);
            }
            for (Integer obj : Q) {
                Q_copy.add(obj);
            }



            P_copy.pop();
            //find the neighbors of the central node from P
            //System.out.println("ut_set.elementAt(a): " + ut_set.elementAt(a));

            N = find_neighbors(ui);

            int N_size = N.size();
            /*
            System.out.println("FIND NEIGHBORS!!!!!!!!!!!!!!");
            for(int gh=0; gh<N_size; gh++){
            System.out.print(N.get(gh) + " ");
            }
            //            System.out.println();
            //            System.out.println();
             */

            //System.out.println("Neighbors: ");

            for (int b = 0; b < N_size; b = b + 2) {
                // N[b] is node v
                // System.out.print( N.get(b) + " ");
                //Grouping of the neighbors:


                int Nelement_at_b = N.get(b);

                //   System.out.println("N["+ b + "]: " + N.elementAt(b) + " " + "Q[" + c + "]: " + Q.elementAt(c));
                if (N.elementAt(b + 1) == 1) {
                    //u and v are adjacent via a C-edge

                    /*if (X.contains(Nelement_at_b)) {
                    X_copy.addElement(Nelement_at_b);
                    }else*/
                    if (Q.contains(Nelement_at_b)) {

                        P_copy.push(Nelement_at_b);
                        //delete N[b] bzw. Q[c] from set Q_copy, remove C-edges
                        Q_copy.removeElement(Nelement_at_b);

                    }
                    if (Y.contains(Nelement_at_b)) {
                        if (X.contains(Nelement_at_b)) {
                            X_copy.addElement(Nelement_at_b);
                        }
                        Y_copy.removeElement(Nelement_at_b);
                    }
                }



                //find respective neighbor position in ut_set, which is needed for the deletion from ut_set

                if (ut_set.indexOf(Nelement_at_b) <= a && ut_set.indexOf(Nelement_at_b) > -1) {
                    --a;
                }

                ut_set.removeElement(Nelement_at_b);

            }
            Stack<Integer> P_copy_N_intersec = new Stack<Integer>();
            Vector<Integer> Q_copy_N_intersec = new Vector<Integer>();
            Vector<Integer> X_copy_N_intersec = new Vector<Integer>();
            Vector<Integer> Y_copy_N_intersec = new Vector<Integer>();
            // System.out.println();
            // System.out.println("P_Copy, N: " + P_copy + " " + N);

            int nElement = -1;


            for (int sec = 0; sec < N_size; sec = sec + 2) {

                nElement = N.get(sec);

                if (P_copy.contains(nElement)) {
                    P_copy_N_intersec.push(nElement);
                }
                if (Q_copy.contains(nElement)) {
                    Q_copy_N_intersec.add(nElement);
                }
                if (X_copy.contains(nElement)) {
                    X_copy_N_intersec.add(nElement);
                }
                if (Y_copy.contains(nElement)) {
                    Y_copy_N_intersec.add(nElement);
                }
            }


            P_copy_N_intersec.push(0);
            R_copy.add(ui);
            Enumerate_Cliques(R_copy, P_copy_N_intersec, Q_copy_N_intersec, X_copy_N_intersec, Y_copy_N_intersec);
            X.add(ui);
            a++;
        }
        return 0;
    }

    private Vector<Integer> find_neighbors(int central_node) {

        Vector<Integer> neighbor_vec = new Vector<Integer>();

        //  System.out.println("C_edge Size: " + C_edges.size());
//        int C_edge_number = C_edges.size() / 2;

        //    System.out.println("C_edge.size/2: " + C_edge_number);
        //    System.out.println("");
        //    System.out.println("C_edges: ");
        for (int a = 0; a < C_edge_Iteration_size; a++) {
            if (C_edges.get(a * 2 + 0) == central_node) {
                //          System.out.println( C_edges.get(a*2+0) + " " + C_edges.get(a*2+1));
                neighbor_vec.add(C_edges.get(a * 2 + 1));
                neighbor_vec.add(1); // 1 means: is connected via C-edge
            }

            if (C_edges.get(a * 2 + 1) == central_node) {
                //           System.out.println(C_edges.get(a*2+0) + " " + C_edges.get(a*2+1));
                neighbor_vec.add(C_edges.get(a * 2 + 0));
                neighbor_vec.add(1); // 1 means: is connected via C-edge
            }

        }



        //System.out.println("");
//        System.out.println("D_edges Size: " + D_edges.size());
//        System.out.println("Reduced D_edges Size: " + D_edge_number);
        for (int a = 0; a < D_edge_Iteration_size; a++) {
            if (D_edges.get(a * 2 + 0) == central_node) {
                //       System.out.println( D_edges.get(a*2+0) + " " + D_edges.get(a*2+1));
                neighbor_vec.add(D_edges.get(a * 2 + 1));
                neighbor_vec.add(2); // 2 means: is connected via D-edge
            }

            if (D_edges.get(a * 2 + 1) == central_node) {
                //        System.out.println(D_edges.get(a*2+0) + " " + D_edges.get(a*2+1));
                neighbor_vec.add(D_edges.get(a * 2 + 0));
                neighbor_vec.add(2); // 2 means: is connected via D-edge
            }

        }

        //System.out.println("Neighbor Edges: " + neighbor_vec);

        return neighbor_vec;
    }

    public int getBestCliqueSize() {
        //=
        //System.out.println(" best_clique_size " +  best_clique_size);
        return best_clique_size;
    }

    public Stack<List<Integer>> getMaxCliqueSet() {
        //System.out.println("Max_Cliques_Set: " + Max_Cliques_Set.size());
        Stack<List<Integer>> solution = new Stack<List<Integer>>();
        solution.addAll(Max_Cliques_Set);
        return solution;
    }
}
