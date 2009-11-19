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

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

/**
 * This class implements Bron-Kerbosch clique detection algorithm as it is
 * described in [F. Cazals, R. Karande: An Algorithm for reporting maximal c-cliques;
 * T.Comp. Sc. (2005); vol 349; pp.
 * 484-490]
 * 
 * @cdk.module smsd
 */
public class BKKCKCF {

    private List<List<Integer>> maxCliquesSet;
    /***********************************************************************/
    private List<Integer> cEdges;
    private List<Integer> dEdges;
    private int bestCliqueSize;
    private List<Integer> comp_graph_nodes;
    private double dEdgeIterationSize = 0;
    private double cEdgeIterationSize = 0;

    /**
     * Creates a new instance of BronKerboschCazalsKarandeKochCliqueFinder
     * @param comp_graph_nodes_org
     * @param C_edges_org C-Edges set of allowed edges
     * @param D_edges_org D-Edges set of prohibited edges
     */
    protected BKKCKCF(List<Integer> comp_graph_nodes_org, List<Integer> C_edges_org, List<Integer> D_edges_org) {

        this.comp_graph_nodes = comp_graph_nodes_org;
        this.cEdges = C_edges_org;
        this.dEdges = D_edges_org;
        bestCliqueSize = 0;
        //Orignal assignment as per paper
        dEdgeIterationSize = dEdges.size() / 2;

        //Orignal assignment as per paper
        cEdgeIterationSize = cEdges.size() / 2;


        boolean d_edgeFlag = false;

        if (dEdges.size() > 10000000 && dEdges.size() > cEdges.size() && cEdges.size() > 100000) {
            dEdgeIterationSize = (float) dEdges.size() * 0.000001;
            d_edgeFlag = true;

        } else if (dEdges.size() > 10000000 && dEdges.size() > cEdges.size() && 5000 < cEdges.size()) {
            dEdgeIterationSize = (float) dEdges.size() * 0.001;
            d_edgeFlag = true;

        }
//        else if (dEdges.size() > 5000000 && dEdges.size() > cEdges.size()) {
//            dEdgeIterationSize = (float) dEdges.size() * 0.0001;
//            d_edgeFlag = true;
//
//        } else if (dEdges.size() > 100000 && dEdges.size() > cEdges.size()) {
//            dEdgeIterationSize = (float) dEdges.size() * 0.1;
//            d_edgeFlag = true;
//        }

//        } else if (dEdges.size() >= 10000 && 500 >= cEdges.size()) {
//            dEdgeIterationSize = (float) dEdges.size() * 0.1;
//            d_edgeFlag = true;
//
//        }
//
//
//
        if (dEdgeIterationSize < 1 && d_edgeFlag && cEdges.size() <= 5000) {

            dEdgeIterationSize = 2;
        }

        if (dEdgeIterationSize < 1 && d_edgeFlag) {

            dEdgeIterationSize = 1;
        }

        //Initialization maxCliquesSet

        maxCliquesSet = new ArrayList<List<Integer>>();

        init();

    }

    /*
     * Call the wrapper for ENUMERATE_CLIQUES
     *
     */
    private void init() {


        /********************************************************************/
        /*
         *T: is a set of vertices which have already been used for the
         * initialization of ENUMERATE_CLIQUES
         */
        List<Integer> T = new ArrayList<Integer>(); //Initialize the T ArrayList;

        /*
         *V: stored all the vertices for the Graph G
         * V[G]
         *nodes of vector comp_graph_nodes are stored in V
         */

        List<Integer> V = new ArrayList<Integer>(); //Initialization of ArrayList V

        int V_set_size = comp_graph_nodes.size() / 3;

        //System.out.println("ArrayList V is initialized");
        for (int a = 0; a < V_set_size; a++) {
            V.add(comp_graph_nodes.get(a * 3 + 2));
            //System.out.print("V[" + a + "]: " + comp_graph_nodes.get(a * 3 + 2) + " ");
        }
        //System.out.println();

        V.add(0);
        // System.out.println("ArrayList V :" + V);

        /*
         * R: set of vertices belonging to the current clique
         */
        List<Integer> R = new ArrayList<Integer>();
        /*
         *P: is a set of vertices which <b>can</b> be added to R, because they are
         * neighbours of vertex u via <i>c-edges</i>
         */
        Stack<Integer> P = new Stack<Integer>();
        /*
         *Q: is a set of vertices which <b>cannot</b> be added to R, because they are
         * neighbours of vertex u via <i>d-edges</i>
         */

        List<Integer> Q = new ArrayList<Integer>();
        /*
         *X: set of vertices which are not allowed to be added
         * to R
         */
        List<Integer> X = new ArrayList<Integer>();


        /*
         *Y: set of vertices which are not allowed to be added
         * to C
         */

        List<Integer> Y = new ArrayList<Integer>();

        /*
         * N[u]: set of neighbours of vertex u in Graph G
         *
         */

        List<Integer> N = new ArrayList<Integer>();

        int b = 0;

        /*
         * Let T be the set of Nodes already been used in the initialization
         *
         */

        T.clear();

        while (V.get(b) != 0) {


            int central_node = V.get(b);


            P.clear();
            Q.clear();
            X.clear();
            R.clear();

            //find the neighbors of the central node from V
            N = findNeighbors(central_node);

            for (int c = 0; c < N.size(); c = c + 2) {

                /*
                 * u and v are adjacent via a R-edge
                 */
                Integer N_at_c = N.get(c);


                //find respective neighbor position in P, which is needed for the deletion from V
                //delete neighbor from set V


                if (N.get(c + 1) == 1) {


                    if (T.contains(N_at_c)) {
                        X.add(N_at_c);
                    } else {
                        P.push(N_at_c);
                    }

                } else if (N.get(c + 1) == 2) {
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
                V.remove(N_at_c);
                //System.out.println("Elements Removed from V:" + N_at_c);
            }

            P.add(0);
            R.add(central_node);

            enumerateCliques(R, P, Q, X, Y);
            //enumerateCliques(R, P, Q, X);
            T.add(central_node);

            b++;
        }
        //System.out.println("maxCliquesSet: " + maxCliquesSet);

    }

    private int enumerateCliques(List<Integer> R, Stack<Integer> P, List<Integer> Q, List<Integer> X, List<Integer> Y) {
        List<Integer> N = new ArrayList<Integer>(); ////Initialization ArrayList N
        Stack<Integer> ut_set = new Stack<Integer>();//Defined as P' in the paper


        for (Integer I : P) {
            ut_set.add(I);
        }

        if ((P.size() == 1) && (X.size() == 0)) {

            //store best solutions in stack maxCliquesSet
            int clique_size = R.size();



            if (clique_size >= bestCliqueSize) {
                if (clique_size > bestCliqueSize) {

                    maxCliquesSet.clear();
                    bestCliqueSize = clique_size;

                }
                if (clique_size == bestCliqueSize) {
                    //System.out.println("R-Clique " + R);
                    maxCliquesSet.add(R);
                }

            }

            return 0;

        }


        int a = 0;

        while (ut_set.elementAt(a) != 0) {

            int ui = ut_set.get(a);

            P.removeElement(ui);

            List<Integer> R_copy = new ArrayList<Integer>(R);
            Stack<Integer> P_copy = new Stack<Integer>();
            Stack<Integer> Q_copy = new Stack<Integer>();
            List<Integer> X_copy = new ArrayList<Integer>(X);
            List<Integer> Y_copy = new ArrayList<Integer>(Y);

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

            N = findNeighbors(ui);

            int N_size = N.size();
            
            //System.out.println("Neighbors: ");

            for (int b = 0; b < N_size; b = b + 2) {
                // N[b] is node v
                //Grouping of the neighbors:


                int Nelement_at_b = N.get(b);

                if (N.get(b + 1) == 1) {
                    //u and v are adjacent via a C-edge

                    if (Q.contains(Nelement_at_b)) {

                        P_copy.push(Nelement_at_b);
                        //delete N[b] bzw. Q[c] from set Q_copy, remove C-edges
                        Q_copy.removeElement(Nelement_at_b);

                    }
                    if (Y.contains(Nelement_at_b)) {
                        if (X.contains(Nelement_at_b)) {
                            X_copy.add(Nelement_at_b);
                        }
                        Y_copy.remove(Nelement_at_b);
                    }
                }

                //find respective neighbor position in ut_set, which is needed for the deletion from ut_set

                if (ut_set.indexOf(Nelement_at_b) <= a && ut_set.indexOf(Nelement_at_b) > -1) {
                    --a;
                }

                ut_set.removeElement(Nelement_at_b);

            }
            Stack<Integer> P_copy_N_intersec = new Stack<Integer>();
            List<Integer> Q_copy_N_intersec = new ArrayList<Integer>();
            List<Integer> X_copy_N_intersec = new ArrayList<Integer>();
            List<Integer> Y_copy_N_intersec = new ArrayList<Integer>();

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
            enumerateCliques(R_copy, P_copy_N_intersec, Q_copy_N_intersec, X_copy_N_intersec, Y_copy_N_intersec);
            X.add(ui);
            a++;
        }
        return 0;
    }

    private List<Integer> findNeighbors(int central_node) {

        List<Integer> neighbor_vec = new ArrayList<Integer>();

        for (int a = 0; a < cEdgeIterationSize; a++) {
            if (cEdges.get(a * 2 + 0) == central_node) {
                //          System.out.println( cEdges.get(a*2+0) + " " + cEdges.get(a*2+1));
                neighbor_vec.add(cEdges.get(a * 2 + 1));
                neighbor_vec.add(1); // 1 means: is connected via C-edge
            }

            if (cEdges.get(a * 2 + 1) == central_node) {
                //           System.out.println(cEdges.get(a*2+0) + " " + cEdges.get(a*2+1));
                neighbor_vec.add(cEdges.get(a * 2 + 0));
                neighbor_vec.add(1); // 1 means: is connected via C-edge
            }

        }
        for (int a = 0; a < dEdgeIterationSize; a++) {
            if (dEdges.get(a * 2 + 0) == central_node) {
                //       System.out.println( dEdges.get(a*2+0) + " " + dEdges.get(a*2+1));
                neighbor_vec.add(dEdges.get(a * 2 + 1));
                neighbor_vec.add(2); // 2 means: is connected via D-edge
            }

            if (dEdges.get(a * 2 + 1) == central_node) {
                //        System.out.println(dEdges.get(a*2+0) + " " + dEdges.get(a*2+1));
                neighbor_vec.add(dEdges.get(a * 2 + 0));
                neighbor_vec.add(2); // 2 means: is connected via D-edge
            }

        }

        return neighbor_vec;
    }

    protected int getBestCliqueSize() {
        return bestCliqueSize;
    }

    protected Stack<List<Integer>> getMaxCliqueSet() {
        Stack<List<Integer>> solution = new Stack<List<Integer>>();
        solution.addAll(maxCliquesSet);
        return solution;
    }
}
