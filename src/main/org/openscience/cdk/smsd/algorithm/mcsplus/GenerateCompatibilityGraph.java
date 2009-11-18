/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
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

import java.io.IOException;
import java.util.List;
import java.util.Vector;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.smsd.global.BondType;
import org.openscience.cdk.smsd.helper.LabelContainer;

/**
 * @cdk.module smsd
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
        this.ac1 = ac1;
        this.ac2 = ac2;
        this.removeHydrogen = HFlag;

        compatibilityGraphNodes();

        if (bondTypeFlag) {
            compatibilityGraphBS();
        } else {
            CompatibilityGraphBIS();
        }

        if (getCEdgesSize() == 0) {
//            System.out.println("C-edges are Zero");
            clearCompGraphNodes();

            clearCEgdes();
            clearDEgdes();

            resetCEdgesSize();
            resetDEdgesSize();

            generate_compatibility_graph_nodes_if_C_edge_number_is_zero();
            if (bondTypeFlag) {
                generate_compatibility_graph_if_C_edge_number_is_zero_BS();
            } else {
                generate_compatibility_graph_if_C_edge_number_is_zero_BIS();
            }

            clearCompGraphNodesCZero();
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

    /**
     * Generate Compatibility Graph Nodes
     *
     * @return
     * @throws IOException 
     */
    protected int compatibilityGraphNodes() throws IOException {

        comp_graph_nodes.clear();
        Vector<IAtom> basic_atom_vec_A = null;
        Vector<IAtom> basic_atom_vec_B = null;
        IAtomContainer reactant = ac1;
        IAtomContainer product = ac2;

        basic_atom_vec_A = reduce_atomset(reactant);
        basic_atom_vec_B = reduce_atomset(product);

        Vector<Vector<Integer>> label_list_molA = label_atoms(reactant);
        Vector<Vector<Integer>> label_list_molB = label_atoms(product);



        int molA_nodes = 0;
        int count_nodes = 1;

        for (Vector<Integer> labelA : label_list_molA) {

            int molB_nodes = 0;

            for (Vector<Integer> labelB : label_list_molB) {
                if (labelA.equals(labelB)) {
                    comp_graph_nodes.addElement(reactant.getAtomNumber(basic_atom_vec_A.get(molA_nodes)));
                    comp_graph_nodes.addElement(product.getAtomNumber(basic_atom_vec_B.get(molB_nodes)));
                    comp_graph_nodes.addElement(count_nodes++);

                }

                molB_nodes++;

            }
            molA_nodes++;
        }
        return 0;
    }

    /**
     * Generate Compatibility Graph Nodes Bond Insensitive
     *
     * @return
     * @throws IOException
     */
    protected int CompatibilityGraphBIS() throws IOException {
        int comp_graph_nodes_vector_size = comp_graph_nodes.size();

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

    /**
     * Generate Compatibility Graph Graph Bond Sensitive
     *
     * @return
     * @throws IOException
     */
    protected int compatibilityGraphBS() throws IOException {
        int comp_graph_nodes_vector_size = comp_graph_nodes.size();

        C_edges = new Vector<Integer>(); //Initialize the C_edges Vector
        D_edges = new Vector<Integer>(); //Initialize the D_edges Vector

        for (int a = 0; a < comp_graph_nodes_vector_size; a = a + 3) {


            int index_a = comp_graph_nodes.elementAt(a);
            int index_aPlus1 = comp_graph_nodes.elementAt(a + 1);

            for (int b = a + 3; b < comp_graph_nodes_vector_size; b = b + 3) {

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

                        if ((ReactantBond.getFlag(CDKConstants.ISAROMATIC) == ProductBond.getFlag(CDKConstants.ISAROMATIC)) && (ReactantBondType.equals(ProductBondType))) {
                            C_edges.addElement((a / 3) + 1);
                            C_edges.addElement((b / 3) + 1);
                        } else if (ReactantBond.getFlag(CDKConstants.ISAROMATIC) && ProductBond.getFlag(CDKConstants.ISAROMATIC)) {
                            C_edges.addElement((a / 3) + 1);
                            C_edges.addElement((b / 3) + 1);
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

    protected Vector<Integer> getCEgdes() {
        return C_edges;
    }

    protected Vector<Integer> getDEgdes() {
        return D_edges;
    }

    protected int getCEdgesSize() {
        return C_edges_size;
    }

    protected int getDEdgesSize() {
        return D_edges_size;
    }

    protected Vector<Integer> getCompGraphNodes() {
        return comp_graph_nodes;
    }

    protected Vector<Integer> getCompGraphNodesCZero() {
        return comp_graph_nodes_C_zero;
    }

    protected void clearCEgdes() {

        C_edges.clear();
    }

    protected void clearDEgdes() {

        D_edges.clear();
    }

    protected void clearCompGraphNodes() {
        comp_graph_nodes.clear();
    }

    protected void clearCompGraphNodesCZero() {
        comp_graph_nodes_C_zero.clear();
    }

    protected void resetCEdgesSize() {
        C_edges_size = 0;
    }

    protected void resetDEdgesSize() {
        D_edges_size = 0;
    }

    protected void clear() {
        C_edges = null;
        D_edges = null;
        comp_graph_nodes = null;
        comp_graph_nodes_C_zero = null;
    }
}
