/*
 * BronKerboschKochCliqueFinder.java
 *
 * Created on October 17, 2007, 4:42 PM
 *
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package org.openscience.cdk.smsd.algorithm.mcsplus;

import java.util.Stack;
import java.util.Vector;

/**
 * This class implements Bron-Kerbosch clique detection algorithm as it is
 * described in [Ina Koch: Enumerating all connected maximal common subgraphs
 * in two graphs; T.Comp. Sc. (2001); vol 250; pp.
 * 1-30]
 *
 * @author Syed Asad Rahman, EMBL-EBI
 */
public class BronKerboschKochCliqueFinder {

    private Stack<Vector<Integer>> Max_Cliques_Set;
    /********************************************************************/
    /*
     *T: is a set of vertices which have already been used for the
     * initialization of ENUMERATE_CLIQUES
     */
    private Vector<Integer> T;
    /*
     *V: stored all the vertices for the Graph G
     * V[G]
     *nodes of vector comp_graph_nodes are stored in V
     */
    private static Stack<Integer> V;
    /***********************************************************************/
    private Vector<Integer> C_edges;
    private Vector<Integer> D_edges;
    private int best_clique_size;
    private Vector<Integer> comp_graph_nodes;
    private double D_edge_Iteration_size = 0;

    /**
     * Creates a new instance of BronKerboschKochCliqueFinder
     * @param comp_graph_nodes_org
     * @param C_edges_org C Edges
     * @param D_edges_org D Edges
     */
    public BronKerboschKochCliqueFinder(Vector<Integer> comp_graph_nodes_org, Vector<Integer> C_edges_org, Vector<Integer> D_edges_org) {

        this.comp_graph_nodes = comp_graph_nodes_org;
        this.C_edges = C_edges_org;
        this.D_edges = D_edges_org;
        best_clique_size = 0;

        //System.out.println("C Edges Size: " + C_edges.size());
        //System.out.println("D Edges Size:" + D_edges.size());




        //Orignal assignment as per paper
        D_edge_Iteration_size = D_edges.size() / 2;
        //Heuristic introduced by Asad


        if (D_edges.size() > 10000000 && D_edges.size() > C_edges.size()) {
            D_edge_Iteration_size = (float) D_edges.size() * 0.00000001;

        } else if (D_edges.size() > 5000000 && D_edges.size() > C_edges.size()) {
            D_edge_Iteration_size = (float) D_edges.size() * 0.0000001;

        } else if (D_edges.size() > 1000000 && D_edges.size() > C_edges.size()) {
            D_edge_Iteration_size = (float) D_edges.size() * 0.00001;

        } else if (D_edges.size() > 100000 && D_edges.size() > C_edges.size()) {
            D_edge_Iteration_size = (float) D_edges.size() * 0.001;
        }
//        else if (D_edges.size() > 10000 && D_edges.size() > C_edges.size()) {
//            D_edge_Iteration_size = (float) D_edges.size() * 0.01;
//
//        } else if (D_edges.size() > 1000 && D_edges.size() > C_edges.size()) {
//            D_edge_Iteration_size = (float) D_edges.size() * 0.1;
//
//        }

        if (D_edge_Iteration_size > 0 && D_edge_Iteration_size < 1 && C_edges.size() < 3000) {

            D_edge_Iteration_size = 1;
        }

//        System.out.println("Optimised D-edges: " + (D_edge_Iteration_size));



        //Initialization Max_Cliques_Set

        Max_Cliques_Set = new Stack<Vector<Integer>>();

        T = new Vector<Integer>(); //Initialize the T Vector

        V = new Stack<Integer>(); //Initialization of Vector V


        //Vector<Integer> V = new Vector<Integer>(); //Initialization of Vector V
        int V_set_size = comp_graph_nodes.size() / 3;

        //System.out.println("Vector V is initialized");
        for (int a = 0; a < V_set_size; a++) {
            V.add(comp_graph_nodes.get(a * 3 + 2));
            //System.out.print(comp_graph_nodes.get(a * 3 + 2) + " ");
        }
        // System.out.println();

        V.add(0); //Anh�gen einer Endekennung

        // System.out.println("Vector V :" + V);

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


        /*
         * C: set of vertices belonging to the current clique
         */
        Vector<Integer> C = new Vector<Integer>();
        /*
         *P: is a set of vertices which <b>can</b> be added to C, because they are
         * neighbours of vertex u via <i>c-edges</i>
         */
        Stack<Integer> P = new Stack<Integer>();
        /*
         *D: is a set of vertices which <b>cannot</b> be added to C, because they are
         * neighbours of vertex u via <i>d-edges</i>
         */

        Vector<Integer> D = new Vector<Integer>();
        /*
         *S: set of vertices which are not allowed to be added
         * to C
         */
        Vector<Integer> S = new Vector<Integer>();

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


            /*System.out.println();
            System.out.println("#########################################");
            System.out.println("Central node " + V.get(b));
             */

            int central_node = V.get(b);

            P.clear();
            D.clear();
            S.clear();
            C.clear();

            //find the neighbors of the central node from V
            N = find_neighbors(V.get(b));

            //System.out.println("N-Neigh: " + N.size());

            //int N_size=N.size();

            //gehe Nachbarn durch und ordne sie zu S, P oder D



            for (int c = 0; c < N.size(); c = c + 2) {
                // N[c] is node v
                // System.out.println("N[" + c  +  "]= " + N.get(c) + " ");
                //Grouping of the neighbors in S,P and D

                /*
                 * u and v are adjacent via a C-edge
                 */
                if (N.elementAt(c + 1) == 1) {

                    /*
                    boolean Nc_belongs_to_T = false;
                    int T_size = T.size();
                    // System.out.println("T_size " +T_size);
                    for (int e = 0; e < T_size; e++) {
                    //    System.out.println("T[" + e  +  "]= " +T.elementAt(e));
                    if (N.elementAt(c) == T.elementAt(e)) {
                    T.contains(N.elementAt(c));
                    S.add(N.get(c));
                    Nc_belongs_to_T = true;
                    }
                    }*/



                    if (T.contains(N.elementAt(c))) {
                        S.add(N.get(c));
                    } else {
                        P.push(N.get(c));
                    }

                } else if (N.elementAt(c + 1) == 2) {
                    // u and v are adjacent via a D-edge
                    //System.out.println("u and v are adjacent via a D-edge: " + N.elementAt(c));
                    D.add(N.get(c));
                }
                //find respective neighbor position in P, which is needed for the deletion from V
                //int V_size = V.size();
                int neighbor_position = -1;

                /* System.out.println();
                System.out.println( "--------------------------");
                System.out.println("V Size: "+ V.size());
                 **/

                //Bug solved due to CPP to Java
                int elementAtC = N.get(c);

                for (int d = 0; d < V.size(); d++) {
                    //System.out.println(" N[c]: " + N.elementAt(c)+ " , V[" +  d + "]: " + V.elementAt(d));
                    if (elementAtC == V.elementAt(d)) {
                        neighbor_position = d;
                        // System.out.println("Nachbarposition " + neighbor_position );
                    }
                }

                //delete neighbor from set V
                if (neighbor_position != -1) {
                    //System.out.println("neighbor_position : " + neighbor_position);
                    for (int e = neighbor_position; e < V.size() - 1; e++) {
                        V.set(e, V.get(e + 1));
                    }


                    /* System.out.println("Vector V :" + V);
                    System.out.println("");
                    System.out.println("V element before POP:" + V.lastElement());
                     */
                    V.pop(); //ask Srikant about v.pop() in C++ Check ASAD

                    /*
                    System.out.println("V size after POP:" + V.size());
                     */

                    if (neighbor_position < b) {
                        b = b - 1;
                        //System.out.println("b: " + b);
                    }
                }
            }
            /*
            System.out.println();
            System.out.println("Sets vor Funktionsaufruf: ");
            int Ssize = S.size();
            System.out.print("Set S: ");
            for(int f=0; f<Ssize; f++){
            System.out.print( S.get(f) + " ");
            }
            System.out.println("");
            int Psize = P.size();
            System.out.print("Set P: ");
            for(int f=0; f<Psize; f++){
            System.out.print(  P.get(f) + " ");
            }*/
            /*
            System.out.println( );
            int Dsize = D.size();
            System.out.print(  "Set D: ");
            for(int f=0; f<Dsize; f++){
            System.out.print(  D.get(f) + " ");
            }
            System.out.println( );
             */
            /*
            int Tsize = T.size();
            System.out.print(  "Set T: ");
            for(int f=0; f<Tsize; f++){
            System.out.print( T.get(f) + " ");
            }
            System.out.println( );
            System.out.println( );
             */
            P.add(0);
            C.add(central_node);


            //System.out.println(" Calling Enumerate_Cliques, Vector Size in CliquesGenerator: ");

            //System.out.println("C: " + C.size() + ", P: " + P.size() + ", D: " + D.size() + ", S: " + S.size());

            Enumerate_Cliques(C, P, D, S);
            T.add(V.get(b));
            b++;

            // System.out.println("HELLO, B: " + b + " V Size: " + V.size() );
        }
        //System.out.println("Max_Cliques_Set: " + Max_Cliques_Set);

    }

    private int Enumerate_Cliques(Vector<Integer> C_org, Stack<Integer> P_org, Vector<Integer> D_org, Vector<Integer> S_org) {

        Vector<Integer> N = new Vector<Integer>(); ////Initialization Vector N
        Stack<Integer> ut_set = new Stack<Integer>();//Defined as P' in the paper

        Vector<Integer> C = new Vector<Integer>(C_org);
        Stack<Integer> P = new Stack<Integer>();
        for (Integer I : P_org) {
            P.add(I);
        }
        Vector<Integer> D = new Vector<Integer>(D_org);
        Vector<Integer> S = new Vector<Integer>(S_org);

        Vector<Integer> C_copy = new Vector<Integer>();
        Stack<Integer> P_copy = new Stack<Integer>();
        Stack<Integer> D_copy = new Stack<Integer>();
        Vector<Integer> S_copy = new Vector<Integer>();


        for (Integer I : P) {
            ut_set.add(I);
        }



        //System.out.println(" Vector C Size: " + C.size());
        /*
        System.out.println("ut-set Funktionsbeginn!!!!!!!!!!! " );
        int utt_count = ut_set.size();
        for(int tuc=0; tuc<utt_count; tuc++){
        System.out.print(ut_set.get(tuc) + " ");
        }
        System.out.println("");
        System.out.println("");
        System.out.println("P-set Funktionsbeginn! " );
        int Pss_count = P.size();
        for(int tuc=0; tuc<Pss_count; tuc++){
        System.out.print(P.get(tuc) + " ");
        }
        System.out.println();
        System.out.println();
        System.out.println( "D-set Funktionsbeginn! " );
        int Dss_count = D.size();
        for(int tuc=0; tuc<Dss_count; tuc++){
        System.out.print( D.get(tuc) + " ");
        }
        System.out.println();
        System.out.println();
        System.out.println( "S-set Funktionsbeginn! ");
        int Sss_count = S.size();
        for(int tuc=0; tuc<Sss_count; tuc++){
        System.out.print(S.get(tuc)+" ");
        }
        System.out.println();
        System.out.println();
         **/

        //breche Suche in Pfad ab, wenn Anzahl der Elemente in C,P und D kleiner ist als best_clique_size
        /* int cur_C_size = C.size();
        int cur_P_size = P.size();
        int cur_D_size = D.size();
        int cur_solution_size = cur_C_size + cur_P_size + cur_D_size;
        if(cur_solution_size < best_clique_size){
        System.out.println("Breche ab, da Clique mit " + cur_solution_size
        + " kleiner als "
        + best_clique_size);
        return 0;
        }*/


        if (P.size() == 1) {
            if (S.size() == 0) {


                /* System.out.println("Clique gefunden!!!!!!!!!!!!!!!!");
                System.out.println( "C Vector Size: " + C.size());
                for(int cl=0; cl<C.size(); cl++){
                System.out.print( C.get(cl) + " ");
                }
                System.out.println();
                 */
                //store best solutions in stack Max_Cliques_Set
                int clique_size = C.size();

                if (clique_size >= best_clique_size) {
                    if (clique_size > best_clique_size) {
                        while (!Max_Cliques_Set.empty()) {
                            Max_Cliques_Set.pop();
                        }
                        best_clique_size = clique_size;
                        //System.out.println("Best Cliques Size: " + best_clique_size + " " + clique_size );
                    }
                    if (clique_size == best_clique_size) {
                        Max_Cliques_Set.push(C);
                    }
                }

                return 0;
            }
        }

        /* int UTT_size = ut_set.size();
        System.out.println("UT SET vor SCHLEIFE!");
        for (int ghh = 0; ghh < UTT_size; ghh++) {
        System.out.print(ut_set.get(ghh) + " ");
        }
        System.out.println("");
        System.out.println("");
         */
        int a = 0;

        while (ut_set.elementAt(a) != 0) {
            // P[a] is node ut
            //  System.out.println("Central recursion node " + ut_set.elementAt(a));
            int ui = ut_set.get(a); //P und ut_set �dern sich -> brauch festen Knoten fr C und S
            //remove ut_set[a] from P
            //find position of ut_set node in P
            int P_size = P.size();
            Integer ut_node_pos = 100000;
            for (int counter = 0; counter < P_size - 1; counter++) {
                //-1 wegen Endekennung

                if (P.elementAt(counter) == ui) {
                    ut_node_pos = counter;
                }
            }
            if (ut_node_pos == 100000) {
                System.out.println("ut_node_pos = 100000");
            }
            //delete ut_set node in P
            for (int counter = ut_node_pos; counter < P_size - 1; counter++) {
                P.set(counter, P.get(counter + 1));
            }

            P.pop(); //POP in CPP Asad

            C_copy.clear();
            P_copy.clear();
            D_copy.clear();
            S_copy.clear();
            N.clear();


            for (Integer obj : C) {
                C_copy.add(obj);
            }

            for (Integer obj : P) {
                P_copy.add(obj);
            }
            for (Integer obj : D) {
                D_copy.add(obj);
            }
            for (Integer obj : S) {
                S_copy.add(obj);
            }

            /*
            C_copy=C;
            P_copy=P;
            D_copy=D;
            S_copy=S;
             */

            /*
            System.out.println();
            System.out.println("P Set: " + P);
             */

            /*
            System.out.println("************");
            System.out.println("Size");
            System.out.println("C_copy: " + C_copy.size());
            System.out.println("D_copy: " + D_copy.size());
            System.out.println("P_copy: " + P_copy.size());
            System.out.println("S_copy: " + S_copy.size());
            System.out.println("************");
             */

            P_copy.pop(); //Entferne Endekennung bei P_copy POP method() in CPP

            /*
            int Pcize = P_copy.size();
            System.out.println("P_copy!");
            for(int gt=0; gt<Pcize; gt++){
            System.out.print( P_copy.get(gt) + " ");
            }
            System.out.println();
            System.out.println();
             */
            //find the neighbors of the central node from P
            //System.out.println("ut_set.elementAt(a): " + ut_set.elementAt(a));

            N = find_neighbors(ut_set.get(a));

            //gehe Nachbarn durch und ordne sie zu S, P oder D
            int N_size = N.size();
            /*
            System.out.println("FIND NEIGHBORS!!!!!!!!!!!!!!");
            for(int gh=0; gh<N_size; gh++){
            System.out.print(N.get(gh) + " ");
            }
            System.out.println();
            System.out.println();
             */

            //System.out.println("Neighbors: ");

            for (int b = 0; b < N_size; b = b + 2) {
                // N[b] is node v
                // System.out.print( N.get(b) + " ");
                //Grouping of the neighbors:

                /*int P_size = P.size();  //Abschnitt sinnlos, denn war etwas in P bleibt es in P'
                for (int c = 0; c < P_size - 1; c++) {
                if (N[b] == P[c]) { //is neighbor N[b] an element of P?
                P_copy.add(N[b]); //so it will be remain in set P_copy
                }
                }*/

                int D_set_size = D.size(); //gehe Nachbarn durch, die ber D-edge mit central_node verbunden sind

                int Nelement_at_b = N.elementAt(b);

                for (int c = 0; c < D_set_size; c++) {

                    //   System.out.println("N["+ b + "]: " + N.elementAt(b) + " " + "D[" + c + "]: " + D.elementAt(c));
                    if (Nelement_at_b == D.elementAt(c)) {
                        //is neighbor N[b] an element of D?
                        //if((303==N[b])&&(303==D[c]))
                        //cout<<N[b]<<" IST NACHBARRR UNNDDDD IINNNN DDDDD SETT !!!"<<endl;
                        if (N.elementAt(b + 1) == 1) {
                            //u and v are adjacent via a C-edge


                            if (T.contains(Nelement_at_b)) {
                                S_copy.add(N.get(b));
                            } else {
                                P_copy.push(N.get(b));
                            }


                            //delete N[b] bzw. D[c] from set D_copy
                            int D_copy_size = D_copy.size();
                            int Nb_position = 10000;
                            for (int e = 0; e < D_copy_size; e++) {
                                //suche Knoten N[b] in D_copy
                                if (Nelement_at_b == D_copy.elementAt(e)) {
                                    Nb_position = e;
                                }
                            }
                            for (int e = Nb_position; e < D_copy_size - 1; e++) {
                                D_copy.set(e, D_copy.get(e + 1));
                                //delete N[b] bzw. D[c] from set D_copy
                            }

                            D_copy.pop();
                        }
                        /*//Abschnitt sinnlos, denn wenn etwas in S war ist, es nach S' kopiert worden
                        if(N[b+1] == 2){     //u and v are adjacent via a D-edge
                        if().....
                        }*/
                    }
                }
                //find respective neighbor position in ut_set, which is needed for the deletion from ut_set
                int ut_set_size = ut_set.size();
                int neighbor_position = -1;
                for (int e = 0; e < ut_set_size; e++) {
                    if (Nelement_at_b == ut_set.elementAt(e)) {
                        neighbor_position = e;
                        //System.out.println ("Nachbarposition " + neighbor_position );
                    }
                }
                if (neighbor_position != -1) {
                    //neighbor der nicht in P vorkommt, wrde sonst Elemente l�chen
                    //delete neighbor from set P
                    for (int e = neighbor_position; e < ut_set_size - 1; e++) {
                        ut_set.set(e, ut_set.get(e + 1));
                    }
                    ut_set.pop(); //TODO:Check removeElementsAt to see whether size returns number of elements or index value
                    if (neighbor_position < a) {
                        a = a - 1;
                    } //wenn neighbor_position<a, wrde sonst ein ut bersprungen
                }
            }
            /*
            int dcopy = D_copy.size();
            cout << "D_copy vor Schnittmengenbildung:" << endl;
            for(int dc=0; dc<dcopy; dc++){
            cout << D_copy[dc] << " ";
            }
            cout << endl << endl;
             */
            //Erstellen der Schnittmengen:
            Stack<Integer> P_copy_N_intersec = new Stack<Integer>();
            Vector<Integer> D_copy_N_intersec = new Vector<Integer>();
            Vector<Integer> S_copy_N_intersec = new Vector<Integer>();



            // System.out.println();
            // System.out.println("P_Copy, N: " + P_copy + " " + N);

            int nElement = -1;

            for (int sec = 0; sec < N_size; sec = sec + 2) {

                nElement = N.get(sec);

                if (P_copy.contains(nElement)) {
                    P_copy_N_intersec.push(nElement);
                }
                if (D_copy.contains(nElement)) {
                    D_copy_N_intersec.add(nElement);
                }
                if (S_copy.contains(nElement)) {
                    S_copy_N_intersec.add(nElement);
                }

                /*
                for (int pc = 0; pc < P_copy.size(); pc++) {
                if (P_copy.elementAt(pc) == nElement) {
                P_copy_N_intersec.push(P_copy.get(pc));
                //System.out.println("HEllo 220");
                //System.out.println("Adding: " + P_copy_N_intersec);
                }
                }
                for (int pd = 0; pd < D_copy.size(); pd++) {
                if (D_copy.elementAt(pd) == nElement) {
                D_copy_N_intersec.add(D_copy.get(pd));
                //System.out.println("Adding: " + D_copy_N_intersec);
                }
                }
                for (int ps = 0; ps < S_copy.size(); ps++) {
                if (S_copy.elementAt(ps) == nElement) {
                S_copy_N_intersec.add(S_copy.get(ps));
                }
                }*/

            }
            P_copy_N_intersec.add(0); // Anh�gen der Endekennung 0 bei P
            C_copy.add(ui);

            /*
            System.out.println("");
            System.out.println("Sets vor Funktionsaufruf: ");
            System.out.print("Set P-N: ");
            int PE_size = P_copy_N_intersec.size();
            for (int f = 0; f < PE_size; f++) {
            System.out.print(P_copy_N_intersec.get(f) + " ");
            }
            System.out.println();
            System.out.print("Set D-N: ");
            int DE_size = D_copy_N_intersec.size();
            System.out.println(DE_size);
            for(int f=0; f<DE_size; f++){
            System.out.print( D_copy_N_intersec.get(f) + " ");
            //System.out.print( D_copy_N_intersec.size() + " ");
            }
            System.out.println();
            System.out.print("Set S-N: ");
            int SE_size = S_copy_N_intersec.size();
            for(int f=0; f<SE_size; f++){
            System.out.print( S_copy_N_intersec.get(f) + " ");
            }
            System.out.println();
            System.out.println();
            System.out.println();
             */


            Enumerate_Cliques(C_copy, P_copy_N_intersec, D_copy_N_intersec, S_copy_N_intersec);
            S.add(ui);
            a++;
        }

        /*cout << "ut-set VORHER!!!!!!!!!!! " << endl;
        int ut_count = ut_set.size();
        for(int tuc=0; tuc<ut_count; tuc++){
        cout << ut_set[tuc] << " ";
        }
        cout << endl << endl;
        cout << "P VORHER! " << endl;
        int Ps_count = P.size();
        for(int tuc=0; tuc<Ps_count; tuc++){
        cout << P[tuc] << " ";
        }
        cout << endl << endl;
        cout << "D-set VORHER! " << endl;
        int Ds_count = D.size();
        for(int tuc=0; tuc<Ds_count; tuc++){
        cout << D[tuc] << " ";
        }
        cout << endl << endl;
        cout << "S-set VORHER! " << endl;
        int Ss_count = S.size();
        for(int tuc=0; tuc<Ss_count; tuc++){
        cout << S[tuc] << " ";
        }
        cout << endl << endl;*/



        return 0;
    }

    private Vector<Integer> find_neighbors(int central_node) {
        /*
        System.out.println();
        System.out.println("Knoten:  " + central_node + "  hat die Nachbarn: ");
         */
        Vector<Integer> neighbor_vec = new Vector<Integer>();

        //  System.out.println("C_edge Size: " + C_edges.size());
        int C_edge_number = C_edges.size() / 2;

        //    System.out.println("C_edge.size/2: " + C_edge_number);
        //    System.out.println("");
        //    System.out.println("C_edges: ");
        for (int a = 0; a < C_edge_number; a++) {
            if (C_edges.elementAt(a * 2 + 0) == central_node) {
                //          System.out.println( C_edges.get(a*2+0) + " " + C_edges.get(a*2+1));
                neighbor_vec.add(C_edges.get(a * 2 + 1));
                neighbor_vec.add(1); // 1 means: is connected via C-edge
            }
            if (C_edges.elementAt(a * 2 + 1) == central_node) {
                //           System.out.println(C_edges.get(a*2+0) + " " + C_edges.get(a*2+1));
                neighbor_vec.add(C_edges.get(a * 2 + 0));
                neighbor_vec.add(1); // 1 means: is connected via C-edge
            }
        }


        //System.out.println("");
//        System.out.println("D_edges Size: " + D_edges.size());
//        System.out.println("Reduced D_edges Size: " + D_edge_number);
        for (int a = 0; a < D_edge_Iteration_size; a++) {
            if (D_edges.elementAt(a * 2 + 0) == central_node) {
                //       System.out.println( D_edges.get(a*2+0) + " " + D_edges.get(a*2+1));
                neighbor_vec.add(D_edges.get(a * 2 + 1));
                neighbor_vec.add(2); // 2 means: is connected via D-edge
            }

            if (D_edges.elementAt(a * 2 + 1) == central_node) {
                //        System.out.println(D_edges.get(a*2+0) + " " + D_edges.get(a*2+1));
                neighbor_vec.add(D_edges.get(a * 2 + 0));
                neighbor_vec.add(2); // 2 means: is connected via D-edge
            }

        }

//		int D_edge_number = D_edges.size() / 2;
//
//		//System.out.println("");
//		//  System.out.println("D_edges Size: "+ D_edges.size());
//		for (int a = 0; a < D_edge_number; a++) {
//			if (D_edges.elementAt(a * 2 + 0) == central_node) {
//				//       System.out.println( D_edges.get(a*2+0) + " " + D_edges.get(a*2+1));
//				neighbor_vec.add(D_edges.get(a * 2 + 1));
//				neighbor_vec.add(2); // 2 means: is connected via D-edge
//			}
//			if (D_edges.elementAt(a * 2 + 1) == central_node) {
//				//        System.out.println(D_edges.get(a*2+0) + " " + D_edges.get(a*2+1));
//				neighbor_vec.add(D_edges.get(a * 2 + 0));
//				neighbor_vec.add(2); // 2 means: is connected via D-edge
//			}
//		}

        //System.out.println("Neighbor Edges: " + neighbor_vec);

        return neighbor_vec;
    }

    /**
     * 
     * @return
     */
    public int getBestCliqueSize() {
        return best_clique_size;
    }

    /**
     * 
     * @return
     */
    public Stack<Vector<Integer>> getMaxCliqueSet() {
        //System.out.println("Max_Cliques_Set: " + Max_Cliques_Set.size());
        return Max_Cliques_Set;
    }
}
