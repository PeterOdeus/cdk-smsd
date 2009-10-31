/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.algorithm.mcsplus;


import java.io.IOException;
import java.util.List;
import java.util.Vector;
import org.openscience.cdk.smsd.global.BondType;
import org.openscience.cdk.smsd.helper.LabelContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;

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
public class GenerateCompatibilityGraph {

    private Vector<Integer> comp_graph_nodes = new Vector<Integer>();
    private Vector<Integer> comp_graph_nodes_C_zero = new Vector<Integer>();
    private Vector<Integer> C_edges = new Vector<Integer>();
    private Vector<Integer> D_edges = new Vector<Integer>();
    private int C_edges_size = 0;
    private int D_edges_size = 0;
    private boolean removeHydrogen = false;
    private IAtomContainer ac1 = null;
    private IAtomContainer ac2 = null;
//    private static int min_atoms_size = 0;
    private boolean bondTypeFlag = BondType.getInstance().getBondSensitiveFlag();

    /**
     * 
     * @param ac1
     * @param ac2
     * @param HFlag 
     * @throws java.io.IOException
     */
    public GenerateCompatibilityGraph(IAtomContainer ac1, IAtomContainer ac2, boolean HFlag) throws IOException {



//        System.out.println("\n\nHello I am in generate_compatibility_graph\n");
        this.ac1 = ac1;
        this.ac2 = ac2;
        this.removeHydrogen = HFlag;

        generate_compatibility_graph_nodes();

        if (bondTypeFlag) {
            generate_compatibility_graph_BS();
        } else {
            generate_compatibility_graph_BIS();
        }

        if (get_C_edges_size() == 0) {
//            System.out.println("C-edges are Zero");
            Clear_CompGraphNodes();

            Clear_C_Egdes();
            Clear_D_Egdes();

            Reset_C_EdgesSize();
            Reset_D_EdgesSize();

            generate_compatibility_graph_nodes_if_C_edge_number_is_zero();
            if (bondTypeFlag) {
                generate_compatibility_graph_if_C_edge_number_is_zero_BS();
            } else {
                generate_compatibility_graph_if_C_edge_number_is_zero_BIS();
            }

            Clear_CompGraphNodes_C_Zero();
        }

    }

    private Vector<Vector<Integer>> label_atoms(IAtomContainer ac) {
        Vector<Vector<Integer>> label_list = new Vector<Vector<Integer>>();


        for (int i = 0; i < ac.getAtomCount(); i++) {
            Vector<Integer> label = new Vector<Integer>();

            LabelContainer labelContainer = LabelContainer.getInstance();

            label.setSize(7);

            for (int a = 0; a < 7; a++) {
                label.setElementAt(0, a);
            }

            IAtom refAtom = ac.getAtom(i);
            String atom1_type = refAtom.getSymbol();



            label.setElementAt(labelContainer.getLabelID(atom1_type), 0);

            int count_neighbors = 1;


            List<IAtom> connAtoms = ac.getConnectedAtomsList(refAtom);

            for (IAtom negAtom : connAtoms) {
                String atom2_type = negAtom.getSymbol();
                label.setElementAt(labelContainer.getLabelID(atom2_type), count_neighbors++);
            }

            bubble_sort(label);
            label_list.add(label);

        }

        return label_list;

    }

    private void bubble_sort(Vector<Integer> label) {

        boolean flag = true; // set flag to 1 to begin initial pass

        int temp; // holding variable

        for (int i = 0; i < 7 && flag; i++) {
            flag = false;
            for (int j = 0; j < 6; j++) {
                if (label.get(i) > label.get(j + 1)) {
                    // descending order simply changes to >
                    temp = label.get(i); // swap elements

                    label.setElementAt(label.get(j + 1), i);
                    label.setElementAt(temp, j + 1);
                    flag = true; // indicates that a swap occurred.

                }
            }
        }
        return; //arrays are passed to functions by address; they are not returned

    }

    private Vector<IAtom> reduce_atomset(IAtomContainer ac) {

        Vector<IAtom> basic_atoms = new Vector<IAtom>();
        for (IAtom atom : ac.atoms()) {
            if (removeHydrogen) {
                if (!atom.getSymbol().equalsIgnoreCase("H")) {
                    basic_atoms.addElement(atom);
                }
            } else {
                basic_atoms.addElement(atom);
            }
        }
        return basic_atoms;
    }

//    private Vector<IAtom> reduce_atomset(IAtomContainer ac) {
//
//        Vector<IAtom> phosphate_O_atoms = new Vector<IAtom>();
//        Vector<IAtom> H_atoms = new Vector<IAtom>();
//
//
////        System.out.println("atom_num: " + ac.getAtomCount());
//
//        for (IAtom atom : ac.atoms()) {
//            if (atom.getSymbol().compareToIgnoreCase("O") == 0) {
//                int O_neighbor_num = 0;
//                boolean P_neighbor = false;
//
//                for (IBond bond : ac.bonds()) {
//                    if (atom.equals(bond.getAtom(0))) {
//                        O_neighbor_num++;
//                        if (((bond.getAtom(1)).getSymbol().compareToIgnoreCase("P") == 0) && (bond.getOrder().compareTo(Order.DOUBLE) == 0)) {
//                            P_neighbor = true;
//                        }
//                    }
//                    if (atom.equals(bond.getAtom(1))) {
//                        O_neighbor_num++;
//                        if (((bond.getAtom(0)).getSymbol().compareToIgnoreCase("P") == 0) && (bond.getOrder().compareTo(Order.DOUBLE) == 0)) {
//                            P_neighbor = true;
//                        }
//                    }
//                }
//                if ((O_neighbor_num == 1) && (P_neighbor)) {
//                    phosphate_O_atoms.addElement(atom);
//                }
//            }
//            if (removeHydrogen) {
//                if (atom.getSymbol().equalsIgnoreCase("H")) {
//                    H_atoms.addElement(atom);
//                }
//            }
//        }
//
//        Vector<IAtom> basic_atoms = new Vector<IAtom>();
//        int phosphate_O_atoms_size = phosphate_O_atoms.size();
//        int H_atoms_size = H_atoms.size();
//
//
////        System.out.println("phosphate_O_atoms_size: " + phosphate_O_atoms_size);
////        System.out.println("H_atoms_size: " + H_atoms_size);
//
//
//        for (IAtom atom : ac.atoms()) {
//            boolean no_P_O_atom = true;
//            for (int b = 0; b < phosphate_O_atoms_size; b++) {
//                if (atom.equals(phosphate_O_atoms.get(b))) {
//                    no_P_O_atom = false;
//                }
//            }
//
//            boolean no_H_atom = true;
//            for (int b = 0; b < H_atoms_size; b++) {
//                if (atom.equals(H_atoms.elementAt(b))) {
//                    no_H_atom = false;
//                }
//            }
//
//            if ((no_P_O_atom) && (no_H_atom)) {
//                basic_atoms.addElement(atom);
//            }
//        }
//
//
//
////        System.out.println("basic_atoms: " + basic_atoms.size());
////        Iterator it = basic_atoms.iterator();
////        while (it.hasNext()) {
////            System.out.print((Integer) it.next() + " ");
////        }
////        System.out.println();
//
//
//        return basic_atoms;
//    }
    protected int generate_compatibility_graph_nodes() throws IOException {

        comp_graph_nodes.clear();

        //System.out.println("Hello I am in generate_compatibility_graph_nodes");

        Vector<IAtom> basic_atom_vec_A = null;
        Vector<IAtom> basic_atom_vec_B = null;

        /* ASAD Please check CPP list to Java Iterator*/

        IAtomContainer reactant = ac1;
        IAtomContainer product = ac2;

//        System.out.println("ac1: " + ac1.getAtomCount());
//        System.out.println("ac2: " + ac2.getAtomCount());

        basic_atom_vec_A = reduce_atomset(reactant);
        basic_atom_vec_B = reduce_atomset(product);

        Vector<Vector<Integer>> label_list_molA = label_atoms(reactant);
        Vector<Vector<Integer>> label_list_molB = label_atoms(product);



        int molA_nodes = 0;
        int count_nodes = 1;

//        System.out.println("basic_atom_vec_A " + basic_atom_vec_A);

        for (Vector<Integer> labelA : label_list_molA) {
//
//            System.out.println("labelA Vector: " + labelA.size());
//            System.out.println("Size of the Label A: ");
//            System.out.println(labelA);

            int molB_nodes = 0;

            for (Vector<Integer> labelB : label_list_molB) {

//                System.out.println("With : " + labelB.size() + " Label B: " + labelB);

                if (labelA.equals(labelB)) {
//
//                    System.out.println("Node Counter: " + count_nodes);
//                    System.out.println("labelA == labelB ");
//                    System.out.println("Molecule 1 Label Vector: ");
//                    System.out.println(labelA);
//                    System.out.println("Molecule 2 Label Vector: ");
//                    System.out.println(labelB);
//                    System.out.println("molA_nodes  " + molA_nodes);
//                    System.out.println("molB_nodes  " + molB_nodes);
//
//                    System.out.println("basic_atom_vec_A  " + basic_atom_vec_A.size());
//                    System.out.println("basic_atom_vec_B " + basic_atom_vec_B.size());
////
//
//                    System.out.println("basic_atom_vec_A.get(molA_nodes)  " + basic_atom_vec_A.get(molA_nodes).getSymbol());
//                    System.out.println("basic_atom_vec_B.get(molB_nodes) " + basic_atom_vec_B.get(molB_nodes).getSymbol());

                    comp_graph_nodes.addElement(reactant.getAtomNumber(basic_atom_vec_A.get(molA_nodes)));
                    comp_graph_nodes.addElement(product.getAtomNumber(basic_atom_vec_B.get(molB_nodes)));
                    comp_graph_nodes.addElement(count_nodes++);
                    
                }

                molB_nodes++;

            }
            molA_nodes++;
        }
//        System.out.println("Map size: " + map.size());
//        System.out.println("count_nodes size: " + count_nodes);
//        System.out.println("comp_graph_nodes size: " + comp_graph_nodes.size() / 3);
//
//        map.clear();
////

        return 0;
    }
//
//    protected int generate_compatibility_graph_nodes() throws IOException {
//
//
//        Vector<String> map = new Vector<String>();
//        comp_graph_nodes.clear();
//
//
//        //Asad has rewritten this loop and it gives better result than the older one
//        int count_nodes = 1;
//        for (int i = 0; i < ac1.getAtomCount(); i++) {
//            IAtom atom1 = ac1.getAtom(i);
//            for (int j = 0; j < ac2.getAtomCount(); j++) {
//                IAtom atom2 = ac2.getAtom(j);
//
//                //You can also check object equal or charge, hydrogen count etc
//
//                if (atom1.getSymbol().equalsIgnoreCase(atom2.getSymbol())) {
//                    if (removeHydrogen) {
//                        if (!atom1.getSymbol().equalsIgnoreCase("H") || !atom2.getSymbol().equalsIgnoreCase("H")) {
//                            if (!map.contains(i + "_" + j)) {
//
//                                comp_graph_nodes.addElement(i);
//                                comp_graph_nodes.addElement(j);
//                                comp_graph_nodes.addElement(count_nodes++);
//                                map.add(i + "_" + j);
//                            }
//                        }
//                    } else {
//                        if (!map.contains(i + "_" + j)) {
//                            comp_graph_nodes.addElement(i);
//                            comp_graph_nodes.addElement(j);
//                            comp_graph_nodes.addElement(count_nodes++);
//                            map.add(i + "_" + j);
//
//                        }
//                    }
//                }
//            }
//        }
//
//        System.out.println("Map size: " + map.size());
//        System.out.println("count_nodes size: " + count_nodes);
//        System.out.println("comp_graph_nodes size: " + comp_graph_nodes.size() / 3);
//
//        map.clear();
//
//        return 0;
//    }

    protected int generate_compatibility_graph_BIS() throws IOException {

//        System.out.println("Hello I am in generate_compatibility_graph");


        int comp_graph_nodes_vector_size = comp_graph_nodes.size();



        //int ReactantBondType = 0;
        //int ProductBondType = 0;

//        System.out.println("comp_graph_nodes Size: " + comp_graph_nodes_vector_size);

        C_edges = new Vector<Integer>(); //Initialize the C_edges Vector

        D_edges = new Vector<Integer>(); //Initialize the D_edges Vector

        for (int a = 0; a < comp_graph_nodes_vector_size; a = a + 3) {


            int index_a = comp_graph_nodes.elementAt(a);
            int index_aPlus1 = comp_graph_nodes.elementAt(a + 1);

            for (int b = a + 3; b < comp_graph_nodes_vector_size; b = b + 3) {


                int index_b = comp_graph_nodes.elementAt(b);
                int index_bPlus1 = comp_graph_nodes.elementAt(b + 1);

                // if element ac !=b and atoms on the adjacent sides of the bonds are not equal
                if (a != b && index_a != index_b &&
                        index_aPlus1 != index_bPlus1) {


                    /*ASAD-Comment this part as this module does not care about bond type sensitivity*/

                    IBond ReactantBond = null;
                    IBond ProductBond = null;

                    ReactantBond = ac1.getBond(ac1.getAtom(index_a), ac1.getAtom(index_b));
                    ProductBond = ac2.getBond(ac2.getAtom(index_aPlus1), ac2.getAtom(index_bPlus1));

                    if (ReactantBond != null && ProductBond != null) {

                        C_edges.addElement((a / 3) + 1);
                        C_edges.addElement((b / 3) + 1);


                    } else if (ReactantBond == null && ProductBond == null) {

                        D_edges.addElement((a / 3) + 1);
                        D_edges.addElement((b / 3) + 1);
                    }

                }

                //print C and D edges of the compatibility graph
                C_edges_size = C_edges.size();
                D_edges_size = D_edges.size();

            }
        }
        return 0;
    }

    protected int generate_compatibility_graph_BS() throws IOException {

//        System.out.println("Hello I am in generate_compatibility_graph");

        int comp_graph_nodes_vector_size = comp_graph_nodes.size();

        C_edges = new Vector<Integer>(); //Initialize the C_edges Vector

        D_edges = new Vector<Integer>(); //Initialize the D_edges Vector

        for (int a = 0; a < comp_graph_nodes_vector_size; a = a + 3) {


            int index_a = comp_graph_nodes.elementAt(a);
            int index_aPlus1 = comp_graph_nodes.elementAt(a + 1);

            for (int b = a + 3; b < comp_graph_nodes_vector_size; b = b + 3) {

                // for (int b = a+3; b < comp_graph_nodes_vector_size; b = b + 3) {

                int index_b = comp_graph_nodes.elementAt(b);
                int index_bPlus1 = comp_graph_nodes.elementAt(b + 1);

                // if element a !=b and atoms on the adjacent sides of the bonds are not equal
                if (a != b &&
                        index_a != index_b &&
                        index_aPlus1 != index_bPlus1) {

                    IBond ReactantBond = null;
                    IBond ProductBond = null;

                    ReactantBond = ac1.getBond(ac1.getAtom(index_a), ac1.getAtom(index_b));
                    ProductBond = ac2.getBond(ac2.getAtom(index_aPlus1), ac2.getAtom(index_bPlus1));

                    if (ReactantBond != null && ProductBond != null) {

                        Order ReactantBondType = ReactantBond.getOrder();//.ordinal();

                        Order ProductBondType = ProductBond.getOrder();//.ordinal();

//                            System.out.println();
//                            System.out.println("ReactantBondNo: " + a + " Order:" + ReactantBondType.ordinal() + " ProductBondNo: " + b + " Order:" + ProductBondType.ordinal());
//                            System.out.println("+++++++++++++++++++++++++++");

                        if ((ReactantBond.getFlag(CDKConstants.ISAROMATIC) == ProductBond.getFlag(CDKConstants.ISAROMATIC)) && (ReactantBondType.equals(ProductBondType))) {
                            //System.out.println("C-Edges A=B " + a + " " +  ReactantBondType + " = " +  b + " " + ProductBondType);
                            C_edges.addElement((a / 3) + 1);
                            C_edges.addElement((b / 3) + 1);

//                                System.out.println("\nac1.getAtom(index_a) " + ac1.getAtom(index_a).getSymbol());
//                                System.out.println("ac1.getAtom(index_b) " + ac1.getAtom(index_b).getSymbol());
//                                System.out.println("ac2.getAtom(index_aPlus1) " + ac2.getAtom(index_aPlus1).getSymbol());
//                                System.out.println("ac2.getAtom(index_aPlus1) " + ac2.getAtom(index_aPlus1).getSymbol());
//                                System.out.println("-----------------------");
                        } else if (ReactantBond.getFlag(CDKConstants.ISAROMATIC) && ProductBond.getFlag(CDKConstants.ISAROMATIC)) {
                            //System.out.println("C-Edges isAromatic A=B " + a + " " +  ReactantBond.getFlag(CDKConstants.ISAROMATIC) + " = " +  b + " " + ProductBond.getFlag(CDKConstants.ISAROMATIC));

                            C_edges.addElement((a / 3) + 1);
                            C_edges.addElement((b / 3) + 1);

//                                System.out.println("\nac1.getAtom(index_a) " + ac1.getAtom(index_a).getSymbol());
//                                System.out.println("ac1.getAtom(index_b) " + ac1.getAtom(index_b).getSymbol());
//                                System.out.println("ac2.getAtom(index_aPlus1) " + ac2.getAtom(index_aPlus1).getSymbol());
//                                System.out.println("ac2.getAtom(index_aPlus1) " + ac2.getAtom(index_aPlus1).getSymbol());
//                                System.out.println("-----------------------");
                        } else {

                            D_edges.addElement((a / 3) + 1);
                            D_edges.addElement((b / 3) + 1);
                        }
//

                    }
                }

                //print C and D edges of the compatibility graph
                C_edges_size = C_edges.size();
                D_edges_size = D_edges.size();




            }
        }

        return 0;
    }

    //comp_graph_nodes_C_zero is used to build up of the edges of the compatibility graph
    protected Integer generate_compatibility_graph_nodes_if_C_edge_number_is_zero() throws IOException {

        int count_nodes = 1;
        Vector<String> map = new Vector<String>();
        comp_graph_nodes_C_zero = new Vector<Integer>(); //Initialize the comp_graph_nodes_C_zero Vector


        LabelContainer labelContainer = LabelContainer.getInstance();


        // resets the target graph.
        comp_graph_nodes.clear();

        for (int i = 0; i < ac1.getAtomCount(); i++) {
            for (int j = 0; j < ac2.getAtomCount(); j++) {
                IAtom atom1 = ac1.getAtom(i);
                IAtom atom2 = ac2.getAtom(j);

                //You can also check object equal or charge, hydrogen count etc

                if (atom1.getSymbol().equalsIgnoreCase(atom2.getSymbol())) {
                    if (removeHydrogen) {
                        if (!atom1.getSymbol().equalsIgnoreCase("H") || !atom2.getSymbol().equalsIgnoreCase("H")) {
                            if (!map.contains(i + "_" + j)) {
                                comp_graph_nodes_C_zero.addElement(i);
                                comp_graph_nodes_C_zero.addElement(j);
                                comp_graph_nodes_C_zero.addElement(labelContainer.getLabelID(atom1.getSymbol())); //i.e C is label 1
                                comp_graph_nodes_C_zero.addElement(count_nodes);
                                comp_graph_nodes.addElement(i);
                                comp_graph_nodes.addElement(j);
                                comp_graph_nodes.addElement(count_nodes++);
                                map.add(i + "_" + j);

                            }

                        }
                    } else {
                        if (!map.contains(i + "_" + j)) {
                            comp_graph_nodes_C_zero.addElement(i);
                            comp_graph_nodes_C_zero.addElement(j);
                            comp_graph_nodes_C_zero.addElement(labelContainer.getLabelID(atom1.getSymbol())); //i.e C is label 1

                            comp_graph_nodes_C_zero.addElement(count_nodes);
                            comp_graph_nodes.addElement(i);
                            comp_graph_nodes.addElement(j);
                            comp_graph_nodes.addElement(count_nodes++);
                            map.add(i + "_" + j);

                        }

                    }
                }
            }
        }

        map.clear();
        return count_nodes;
    }

    protected int generate_compatibility_graph_if_C_edge_number_is_zero_BS() throws IOException {


        int comp_graph_nodes_C_zero_vector_size = comp_graph_nodes_C_zero.size();



        //int ReactantBondType = 0;
        //int ProductBondType = 0;

        //System.out.println("Vector_Size: " + comp_graph_nodes_vector_size);

        C_edges = new Vector<Integer>(); //Initialize the C_edges Vector

        D_edges = new Vector<Integer>(); //Initialize the D_edges Vector

        for (int a = 0; a < comp_graph_nodes_C_zero_vector_size; a = a + 4) {


            int index_a = comp_graph_nodes_C_zero.elementAt(a);
            int index_aPlus1 = comp_graph_nodes_C_zero.elementAt(a + 1);



            for (int b = 0; b < comp_graph_nodes_C_zero_vector_size; b = b + 4) {


                int index_b = comp_graph_nodes_C_zero.elementAt(b);
                int index_bPlus1 = comp_graph_nodes_C_zero.elementAt(b + 1);

                // if element a !=b and atoms on the adjacent sides of the bonds are not equal
                if ((a != b) &&
                        (index_a != index_b) &&
                        (index_aPlus1 != index_bPlus1)) {


//
//                    if (molecule1_pair_connected && molecule2_pair_connected) {

                    int ReactantBondType = 0;
                    int ProductBondType = 0;

                    IBond ReactantBond = ac1.getBond(ac1.getAtom(index_a), ac1.getAtom(index_b));
                    IBond ProductBond = ac2.getBond(ac2.getAtom(index_aPlus1), ac2.getAtom(index_bPlus1));



                    //in case that both molecule pairs are connected a c-edge is generated

                    //The bond type check introduced by Asad

                    if (ReactantBond != null && ProductBond != null) {

                        ReactantBondType = ReactantBond.getOrder().ordinal();

                        ProductBondType = ProductBond.getOrder().ordinal();



                        if (ReactantBond.getFlag(CDKConstants.ISAROMATIC) == ProductBond.getFlag(CDKConstants.ISAROMATIC) && ReactantBondType == ProductBondType) {
                            C_edges.addElement((a / 4) + 1);
                            C_edges.addElement((b / 4) + 1);
                        } else if (ReactantBond.getFlag(CDKConstants.ISAROMATIC) && ProductBond.getFlag(CDKConstants.ISAROMATIC)) {


                            C_edges.addElement((a / 4) + 1);
                            C_edges.addElement((b / 4) + 1);
                        } else {

                            D_edges.addElement((a / 4) + 1);
                            D_edges.addElement((b / 4) + 1);
                        }
                    }

//                    } //This was commented by Asad as the bond type mathching is reduced
//                    else if (!molecule1_pair_connected && !molecule2_pair_connected) {
//                        D_edges.addElement((a / 4) + 1);
//                        D_edges.addElement((b / 4) + 1);
//                    }


                }
            }
        }

        //Size of C and D edges of the compatibility graph
        C_edges_size = C_edges.size();
        D_edges_size = D_edges.size();


        return 0;
    }

    protected int generate_compatibility_graph_if_C_edge_number_is_zero_BIS() throws IOException {

        int comp_graph_nodes_C_zero_vector_size = comp_graph_nodes_C_zero.size();



        //int ReactantBondType = 0;
        //int ProductBondType = 0;

        //System.out.println("Vector_Size: " + comp_graph_nodes_vector_size);

        C_edges = new Vector<Integer>(); //Initialize the C_edges Vector

        D_edges = new Vector<Integer>(); //Initialize the D_edges Vector

        for (int a = 0; a <
                comp_graph_nodes_C_zero_vector_size; a =
                        a + 4) {


            int index_a = comp_graph_nodes_C_zero.elementAt(a);
            int index_aPlus1 = comp_graph_nodes_C_zero.elementAt(a + 1);



            for (int b = a + 4; b < comp_graph_nodes_C_zero_vector_size; b =
                            b + 4) {


                int index_b = comp_graph_nodes_C_zero.elementAt(b);
                int index_bPlus1 = comp_graph_nodes_C_zero.elementAt(b + 1);

                // if element ac !=b and atoms on the adjacent sides of the bonds are not equal
                if ((a != b) &&
                        (index_a != index_b) &&
                        (index_aPlus1 != index_bPlus1)) {


                    IBond ReactantBond = null;
                    IBond ProductBond = null;

                    ReactantBond = ac1.getBond(ac1.getAtom(index_a), ac1.getAtom(index_b));
                    ProductBond = ac2.getBond(ac2.getAtom(index_aPlus1), ac2.getAtom(index_bPlus1));

                    if (ReactantBond != null && ProductBond != null) {

                        C_edges.addElement((a / 4) + 1);
                        C_edges.addElement((b / 4) + 1);


                    } else {

                        D_edges.addElement((a / 4) + 1);
                        D_edges.addElement((b / 4) + 1);
                    }


                }
            }
        }

        //Size of C and D edges of the compatibility graph
        C_edges_size = C_edges.size();
        D_edges_size = D_edges.size();


        return 0;
    }

    protected Vector<Integer> get_C_egdes() {

        return C_edges;
    }

    protected Vector<Integer> get_D_egdes() {

        return D_edges;
    }

    protected int get_C_edges_size() {
        return C_edges_size;
    }

    protected int get_D_edges_size() {
        return D_edges_size;
    }

    protected Vector<Integer> get_comp_graph_nodes() {
        return comp_graph_nodes;
    }

    protected Vector<Integer> get_comp_graph_nodes_C_zero() {
        return comp_graph_nodes_C_zero;
    }

    protected void Clear_C_Egdes() {

        C_edges.clear();
    }

    protected void Clear_D_Egdes() {

        D_edges.clear();
    }

    protected void Clear_CompGraphNodes() {
        comp_graph_nodes.clear();
    }

    protected void Clear_CompGraphNodes_C_Zero() {
        comp_graph_nodes_C_zero.clear();
    }

    protected void Reset_C_EdgesSize() {
        C_edges_size = 0;
    }

    protected void Reset_D_EdgesSize() {
        D_edges_size = 0;
    }

    protected void Clear() {
        C_edges = null;
        D_edges = null;
        comp_graph_nodes = null;
        comp_graph_nodes_C_zero = null;
    }
}
