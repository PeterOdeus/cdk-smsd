/* Copyright (C) 2005-2006 Markus Leber
 *               2006-2009 Syed Asad Rahman {asad@ebi.ac.uk}
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
package org.openscience.cdk.smsd.algorithm.mcgregor;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.Vector;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smsd.global.BondType;
import org.openscience.cdk.smsd.helper.BinaryTree;

/**
 * @cdk.module smsd
 */
public class McGregor {

    private IAtomContainer ac1 = null;
    private IAtomContainer ac2 = null;
    private BinaryTree last = null;
    private BinaryTree first = null;
    private Stack<Vector<Integer>> BESTARCS = null;
    private Vector<Integer> MARCS = null;
    private List<Integer> i_globalA = null;
    private List<Integer> i_globalB = null;
    private List<String> c_globalA = null;
    private List<String> c_globalB = null;
    private Vector<String> c_tab1_copy = null;
    private Vector<String> c_tab2_copy = null;
    private Vector<Integer> i_bond_neighborsA = null;
    private Vector<String> c_bond_neighborsA = null;
    private int nNum_globalA = 0;
    private int nNum_globalB = 0;
    private int bestarcsleft = 0;
    private int globalMCSSize = 0;
    private List<List<Integer>> mappings = null;
    private int neighbor_bondnum_A = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
    private int set_bondnum_A = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
    private int neighbor_bondnum_B = 0; //number of remaining molecule B bonds after the clique search, which are neighbors of the MCS_1
    private int set_bondnum_B = 0; //number of remaining molecule B bonds after the clique search, which aren't neighbors
    /*This should be more or equal to all the atom types*/
    private String[] SignROW = {"$1", "$2", "$3", "$4", "$5", "$6", "$7", "$8", "$9", "$10", "$11", "$12",
        "$13", "$15", "$16", "$17", "$18", "$19", "$20", "$21", "$22", "$23", "$24",
        "$25", "$26", "$27", "$28", "$29", "$30", "$31", "$32", "$33", "$34", "$35", "$36",
        "$37", "$38", "$39", "$40", "$41", "$42", "$43", "$44", "$45", "$46",
        "$47", "$48", "$49", "$50", "$51", "$52", "$53", "$54", "$55"
    };
    protected boolean new_matrix = false;
    protected boolean bondTypeFlag = BondType.getInstance().getBondSensitiveFlag();

    /**
     * Creates a new instance of McGregor
     * @param ac1
     * @param ac2 
     * @param mappings
     */
    public McGregor(IAtomContainer ac1, IAtomContainer ac2, List<List<Integer>> _mappings) {


        this.ac1 = ac1;
        this.ac2 = ac2;
        this.mappings = _mappings;
        this.nNum_globalA = 0;
        this.nNum_globalB = 0;
        bestarcsleft = 0;

        if (_mappings.isEmpty()) {
            this.globalMCSSize = 0;
        } else {
            this.globalMCSSize = _mappings.get(0).size();
        }
        MARCS = new Vector<Integer>();

        BESTARCS = new Stack<Vector<Integer>>();

        //Initialization of global vectors
        i_globalA = new Vector<Integer>();
        i_globalB = new Vector<Integer>();
        c_globalA = new Vector<String>();
        c_globalB = new Vector<String>();

        c_tab1_copy = new Vector<String>();
        c_tab2_copy = new Vector<String>();

        new_matrix = false;



    }

    /**
     *
     * @param best_Mapping_size 
     * @param present_Mapping
     * @throws IOException
     */
    public void startMcGregorIteration(int best_Mapping_size, Map<Integer, Integer> present_Mapping) throws IOException {

        this.globalMCSSize = (best_Mapping_size / 2);
        c_tab1_copy.clear();
        generateCTab1Copy();


        c_tab2_copy.clear();
        generateCTab2Copy();


        //find mapped atoms of both molecules and store these in mapped_atoms
        List<Integer> mapped_atoms = new Vector<Integer>();
//        System.out.println("\nMapped Atoms");
        for (Map.Entry<Integer, Integer> map : present_Mapping.entrySet()) {
//            System.out.println("i:" + map.getKey() + " j:" + map.getValue());
            mapped_atoms.add(map.getKey());
            mapped_atoms.add(map.getValue());
        }
        int mapping_size = present_Mapping.size();


        i_bond_neighborsA = new Vector<Integer>();
        Vector<Integer> i_bond_setA = new Vector<Integer>();
        c_bond_neighborsA = new Vector<String>();
        Vector<String> c_bond_setA = new Vector<String>();

        Vector<Integer> i_bond_neighborsB = new Vector<Integer>();
        Vector<Integer> i_bond_setB = new Vector<Integer>();
        Vector<String> c_bond_neighborsB = new Vector<String>();
        Vector<String> c_bond_setB = new Vector<String>();

        //find unmapped atoms of molecule A
        Vector<Integer> unmapped_atoms_molA = new Vector<Integer>();

        int unmapped_numA = 0;
        boolean atomA_is_unmapped = true;

        for (int a = 0; a < ac1.getAtomCount(); a++) {
            //Atomic list are only numbers from 1 to atom_number1

            for (Integer key : present_Mapping.keySet()) {
                if (key == a) {
                    atomA_is_unmapped = false;
                }
            }


            if (atomA_is_unmapped) {
//                System.out.println("UnMapped Atoms: " + a +
//                        " (" + ac1.getAtom(a).getSymbol() + ")");
                unmapped_atoms_molA.add(unmapped_numA, a);
                unmapped_numA++;
            }
            atomA_is_unmapped = true;
        }


//        System.out.println("neighbor_bondnum_A Before:" + neighbor_bondnum_A);

//        System.out.println("clique Size " + mappingSize);

//        System.out.println("unmapped_atoms_molA: " + unmapped_atoms_molA.size());

        int SR_count = 0;

        QueryProcessor QP = new QueryProcessor(ac1, ac2, c_tab1_copy, c_tab2_copy, SignROW, neighbor_bondnum_A, set_bondnum_A);


        QP.process(unmapped_atoms_molA, mapping_size,
                i_bond_neighborsA, i_bond_setA,
                c_bond_neighborsA, c_bond_setA,
                mapped_atoms, SR_count);

        this.c_tab1_copy = QP.getCtab1();
        this.c_tab2_copy = QP.getCtab2();
        this.set_bondnum_A = QP.geBondnum_A();
        this.neighbor_bondnum_A = QP.getNeighbor_bondnum_A();
        this.SignROW = QP.getSigns();
        this.i_bond_neighborsA = QP.getIBondNeighboursA();
        this.c_bond_neighborsA = QP.getCBondNeighborsA();



        //find unmapped atoms of molecule B
        Vector<Integer> unmapped_atoms_molB = new Vector<Integer>();
        int unmapped_numB = 0;
        boolean atomB_is_unmapped = true;

//        System.out.println("neighbor_bondnum_A After:" + neighbor_bondnum_A);
//
//        System.out.println("\n---------------\n");
        for (int a = 0; a < ac2.getAtomCount(); a++) {
            for (Integer value : present_Mapping.values()) {

                if (a == value) {
                    atomB_is_unmapped = false;
                }
            }
            if (atomB_is_unmapped) {
//                System.out.println("UnMapped Atoms: " + a +
//                        " (" + ac2.getAtom(a).getSymbol() + ")");
                unmapped_atoms_molB.add(unmapped_numB, a);
                unmapped_numB++;
            }
            atomB_is_unmapped = true;
        }

//        System.out.println("unmapped_atoms_molB: " + unmapped_atoms_molB.size());

        //Extract bonds which are related with unmapped atoms of molecule B.
        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
        //new special signs -> the result are two vectors: c_bond_neighborsA and int_bonds_molB, which contain those
        //bonds of molecule B, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule A


        TargetProcessor TP = new TargetProcessor(ac2, c_tab1_copy, c_tab2_copy, SignROW, neighbor_bondnum_B, set_bondnum_B, neighbor_bondnum_A, i_bond_neighborsA, c_bond_neighborsA);


        TP.process(unmapped_atoms_molB,
                mapping_size,
                i_bond_neighborsB,
                i_bond_setB,
                c_bond_neighborsB,
                c_bond_setB,
                mapped_atoms,
                SR_count);

        this.c_tab1_copy = TP.getCtab1();
        this.c_tab2_copy = TP.getCtab2();
        this.set_bondnum_B = TP.geBondnum_B();
        this.neighbor_bondnum_B = TP.getNeighbor_bondnum_B();
        this.SignROW = TP.getSigns();

        boolean dummy = false;

//        System.out.println("neighbor_bondnum_A: " + neighbor_bondnum_A);
//        System.out.println("neighbor_bondnum_B: " + neighbor_bondnum_B);
//        System.out.println("set_bondnum_A: " + set_bondnum_A);
//        System.out.println("set_bondnum_B: " + set_bondnum_B);
//        System.out.println("unmapped_atoms_molA: " + unmapped_atoms_molA.size());
//        System.out.println("unmapped_atoms_molB: " + unmapped_atoms_molB.size());


        iterator(dummy, present_Mapping.size(), mapped_atoms, neighbor_bondnum_A, neighbor_bondnum_B, i_bond_neighborsA, i_bond_neighborsB, c_bond_neighborsA, c_bond_neighborsB, set_bondnum_A, set_bondnum_B, i_bond_setA, i_bond_setB, c_bond_setA, c_bond_setB);

        //System.exit(1); //uncomment to debug


    }

    /**
     * 
     * @param best_Mapping_size 
     * @param clique_vector
     * @param comp_graph_nodes
     * @throws IOException
     */
    public void startMcGregorIteration(int best_Mapping_size, List<Integer> clique_vector, List<Integer> comp_graph_nodes) throws IOException {
        this.globalMCSSize = (best_Mapping_size / 2);
        int SR_count = 0;


        c_tab1_copy.clear();
        generateCTab1Copy();


        c_tab2_copy.clear();
        generateCTab2Copy();


        //System.out.println("C_Tab1_Copy " + c_tab1_copy);
        //System.out.println("C_Tab2_Copy " + c_tab2_copy);

        //find mapped atoms of both molecules and store these in mapped_atoms
        Vector<Integer> mapped_atoms = new Vector<Integer>();

        int mapped_atom_number = 0;
//

        i_bond_neighborsA = new Vector<Integer>();
        Vector<Integer> i_bond_setA = new Vector<Integer>();
        c_bond_neighborsA = new Vector<String>();
        Vector<String> c_bond_setA = new Vector<String>();

        Vector<Integer> i_bond_neighborsB = new Vector<Integer>();
        Vector<Integer> i_bond_setB = new Vector<Integer>();
        Vector<String> c_bond_neighborsB = new Vector<String>();
        Vector<String> c_bond_setB = new Vector<String>();



        int clique_siz = clique_vector.size();
        int vec_size = comp_graph_nodes.size();

//        System.out.println("clique_siz: " + clique_siz);
        //System.out.println("vec_size: " + vec_size);

        int cliqueNumber = 0;

        for (int a = 0; a < clique_siz; a++) {
            //go through all clique nodes
            cliqueNumber = clique_vector.get(a);
            for (int b = 0; b < vec_size; b = b + 3) {
                //go through all nodes in the compatibility graph
                if (cliqueNumber == comp_graph_nodes.get(b + 2)) {
                    mapped_atoms.add(comp_graph_nodes.get(b));
                    mapped_atoms.add(comp_graph_nodes.get(b + 1));
                    mapped_atom_number++;
                }
            }
        }


        //find unmapped atoms of molecule A
        Vector<Integer> unmapped_atoms_molA = new Vector<Integer>();

        int unmapped_numA = 0;
        boolean atomA_is_unmapped = true;

//        System.out.println("Mapped Atoms: " + mapped_atoms);

        for (int a = 0; a < ac1.getAtomCount(); a++) {
            //Atomic list are only numbers from 1 to atom_number1

            for (int b = 0; b < clique_siz; b++) {
                //the number of nodes == number of assigned pairs
                //System.out.println(mapped_atoms.elementAt(b * 2)+" z:" + z);
                if (mapped_atoms.elementAt(b * 2) == a) {
                    atomA_is_unmapped = false;
                }
            }


            if (atomA_is_unmapped == true) {
                unmapped_atoms_molA.addElement(a);
                unmapped_numA++;
            }
            atomA_is_unmapped = true;
        }

//        System.out.println("Not mapped A: ");
//        for (int z = 0; z < unmapped_numA; z++) {
//            System.out.print(unmapped_atoms_molA.get(z) + " ");
//        }
//        System.out.println("");


        //Extract bonds which are related with unmapped atoms of molecule A.
        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
        //new special signs -> the result are two vectors: c_bond_neighborsA and int_bonds_molA, which contain those
        //bonds of molecule A, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule B


//        System.out.println("neighbor_bondnum_A Before:" + neighbor_bondnum_A);


//        System.out.println("clique Size " + clique_siz);


        QueryProcessor QP = new QueryProcessor(ac1, ac2, c_tab1_copy, c_tab2_copy, SignROW, neighbor_bondnum_A, set_bondnum_A);


        QP.process(unmapped_atoms_molA, clique_siz,
                i_bond_neighborsA, i_bond_setA,
                c_bond_neighborsA, c_bond_setA,
                mapped_atoms, SR_count);

        this.c_tab1_copy = QP.getCtab1();
        this.c_tab2_copy = QP.getCtab2();
        this.set_bondnum_A = QP.geBondnum_A();
        this.neighbor_bondnum_A = QP.getNeighbor_bondnum_A();
        this.i_bond_neighborsA = QP.getIBondNeighboursA();
        this.c_bond_neighborsA = QP.getCBondNeighborsA();
        this.SignROW = QP.getSigns();



        //find unmapped atoms of molecule B
        Vector<Integer> unmapped_atoms_molB = new Vector<Integer>();
        int unmapped_numB = 0;
        boolean atomB_is_unmapped = true;


//        System.out.println("neighbor_bondnum_A After:" + neighbor_bondnum_A);


        for (int a = 0; a < ac2.getAtomCount(); a++) {
            for (int b = 0; b < clique_siz; b++) {
                if (a == mapped_atoms.elementAt(b * 2 + 1)) {
                    atomB_is_unmapped = false;
                }
            }
            if (atomB_is_unmapped == true) {
                unmapped_atoms_molB.add(a);
                unmapped_numB++;
            }
            atomB_is_unmapped = true;
        }

//        System.out.println("Not mapped B: ");
//        for (int z = 0; z < unmapped_numB; z++) {
//            System.out.print(unmapped_atoms_molB.get(z) + " ");
//        }
//        System.out.println("");


        //Extract bonds which are related with unmapped atoms of molecule B.
        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
        //new special signs -> the result are two vectors: c_bond_neighborsA and int_bonds_molB, which contain those
        //bonds of molecule B, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule A



        TargetProcessor TP = new TargetProcessor(ac2, c_tab1_copy, c_tab2_copy, SignROW, neighbor_bondnum_B, set_bondnum_B, neighbor_bondnum_A, i_bond_neighborsA, c_bond_neighborsA);


        TP.process(unmapped_atoms_molB,
                clique_siz,
                i_bond_neighborsB,
                i_bond_setB,
                c_bond_neighborsB,
                c_bond_setB,
                mapped_atoms,
                SR_count);

        this.c_tab1_copy = TP.getCtab1();
        this.c_tab2_copy = TP.getCtab2();
        this.set_bondnum_B = TP.geBondnum_B();
        this.neighbor_bondnum_B = TP.getNeighbor_bondnum_B();
        this.SignROW = TP.getSigns();

        boolean dummy = false;

//        System.out.println("Calling iterator with mapped atoms number " + mapped_atom_number);
//        System.out.println("Neighbor");
//        System.out.println(neighbor_bondnum_A + " " + neighbor_bondnum_B);

//        System.out.println("Mapped Atoms before Iterator1: " + mapped_atoms);


        iterator(dummy, mapped_atom_number, mapped_atoms, neighbor_bondnum_A, neighbor_bondnum_B, i_bond_neighborsA, i_bond_neighborsB, c_bond_neighborsA, c_bond_neighborsB, set_bondnum_A, set_bondnum_B, i_bond_setA, i_bond_setB, c_bond_setA, c_bond_setB);


        // System.out.println("iterator Over");

        //System.exit(1); //uncomment to debug


    }

    private int iterator(boolean MAPPING_check,
            int mapped_atom_num,
            List<Integer> mapped_atoms_org,
            int neighbor_bondnum_A,
            int neighbor_bondnum_B,
            Vector<Integer> i_bond_neighbor_atoms_A,
            Vector<Integer> i_bond_neighbor_atoms_B,
            Vector<String> c_bond_neighborsA,
            Vector<String> c_bond_neighborsB,
            int set_num_A, int set_num_B,
            Vector<Integer> i_bond_setA,
            Vector<Integer> i_bond_setB,
            Vector<String> c_bond_setA,
            Vector<String> c_bond_setB) throws IOException {

//        System.out.println("iterator");
//        System.out.println("mapped_atom_num " + mapped_atom_num);
//        System.out.println("mapped_atoms " + mapped_atoms_org);
//        System.out.println("neighbor_bondnum_A " + neighbor_bondnum_A);
//        System.out.println("neighbor_bondnum_B " + neighbor_bondnum_B);
//        System.out.println("i_bond_neighbor_atoms_A " + i_bond_neighbor_atoms_A);
//        System.out.println("i_bond_neighbor_atoms_B " + i_bond_neighbor_atoms_B);
//        System.out.println("c_bond_neighborsA " + c_bond_neighborsA);
//        System.out.println("c_bond_neighborsB " + c_bond_neighborsB);
//        System.out.println("i_bond_setA " + i_bond_setA);
//        System.out.println("i_bond_setB " + i_bond_setB);
//        System.out.println("c_bond_setA " + c_bond_setA);
//        System.out.println("c_bond_setB " + c_bond_setB);





        List<Integer> mapped_atoms = new Vector<Integer>(mapped_atoms_org);

        //check possible mappings:
        boolean no_further_mapping_possible = true;

        for (int row = 0; row < neighbor_bondnum_A; row++) {

            for (int column = 0; column < neighbor_bondnum_B; column++) {
                String G1A = c_bond_neighborsA.get(row * 4 + 0);
                String G2A = c_bond_neighborsA.get(row * 4 + 1);
                String G1B = c_bond_neighborsB.get(column * 4 + 0);
                String G2B = c_bond_neighborsB.get(column * 4 + 1);

                //Get the atom index from the i_bond neighbor vactor

                int Index_I = i_bond_neighbor_atoms_A.get(row * 3 + 0);
                int Index_IPlus1 = i_bond_neighbor_atoms_A.get(row * 3 + 1);

                //Get the atoms
                IAtom R1_A = ac1.getAtom(Index_I);
                IAtom R2_A = ac1.getAtom(Index_IPlus1);
                IBond ReactantBond = ac1.getBond(R1_A, R2_A);

                //Get the atom index from the i_bond neighbor vactor
                int Index_J = i_bond_neighbor_atoms_B.get(column * 3 + 0);
                int Index_JPlus1 = i_bond_neighbor_atoms_B.get(column * 3 + 1);

                //Get the atoms
                IAtom P1_B = ac2.getAtom(Index_J);
                IAtom P2_B = ac2.getAtom(Index_JPlus1);
                IBond ProductBond = ac2.getBond(P1_B, P2_B);


                if (bondTypeFlag && bondMatch(ReactantBond, ProductBond)) {
                    if ((G1A.compareToIgnoreCase(G1B) == 0 && G2A.compareToIgnoreCase(G2B) == 0) || (G1A.compareToIgnoreCase(G2B) == 0 && G2A.compareToIgnoreCase(G1B) == 0)) {

                        no_further_mapping_possible = false;
                    }
//                      
                } else if (!bondTypeFlag) {
//                  System.out.println("Eins bei: " + G1A + " " + G2A + " " + G1B + " " + G2B);

                    if (((G1A.compareToIgnoreCase(G1B) == 0) && (G2A.compareToIgnoreCase(G2B) == 0)) || ((G1A.compareToIgnoreCase(G2B) == 0) && (G2A.compareToIgnoreCase(G1B) == 0))) {

//                    System.out.println("Eins bei: " + G1A + " " + G2A + " " + G1B + " " + G2B);
                        no_further_mapping_possible = false;
                    }


                }


            }
        }

        if (neighbor_bondnum_A == 0 || neighbor_bondnum_B == 0 || MAPPING_check || no_further_mapping_possible) {
//            System.out.println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
//            System.out.println("^^^^^^^^^^^^" + neighbor_bondnum_A +  "^^^^^^^^^^^^^^^^^^^");
//            System.out.println("^^^^^^^^^^^^" + neighbor_bondnum_B +  "^^^^^^^^^^^^^^^^^^^");
//            System.out.println("^^^^^^^^^^^^" + MAPPING_check +  "^^^^^^^^^^^^^^^^^^^");
//            System.out.println("^^^^^^^^^^^^" + no_further_mapping_possible + "^^^^^^^^^^^^^^^^^^^^");
//            System.out.println(" globalMCSSize " + globalMCSSize);
//            System.out.println("best_MAPPING_size: " + mapped_atom_num);
            try {
                if (mapped_atom_num >= globalMCSSize) {
//                    System.out.println("Hello-1");
                    if (mapped_atom_num > globalMCSSize) {
//                        System.out.println("Hello-2");
                        this.globalMCSSize = mapped_atom_num;
//                        System.out.println("best_MAPPING_size: " + globalMCSSize);
                        mappings.clear();
                        //mappings=new Vector<Vector<Integer>>();
                    }

                    mappings.add(mapped_atoms);
//                    System.out.println("mappings " + mappings);
                }




            } catch (Exception ex) {
                ex.printStackTrace();
            }
            //System.out.println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
            //NoMoreMappingFlag = false;
            return 0;
        }

        i_globalA.clear();
        i_globalB.clear();
        c_globalA.clear();
        c_globalB.clear();

        //redefining of global vectors and variables
        nNum_globalA = neighbor_bondnum_A; //N global variable defined
        nNum_globalB = neighbor_bondnum_B; //N global variable defined

        i_globalA.addAll(i_bond_neighbor_atoms_A);
        i_globalB.addAll(i_bond_neighbor_atoms_B);
        c_globalA.addAll(c_bond_neighborsA);
        c_globalB.addAll(c_bond_neighborsB);
        MARCS.clear();

        MARCS.setSize(neighbor_bondnum_A * neighbor_bondnum_B);
        for (int i = 0; i < neighbor_bondnum_A * neighbor_bondnum_B; i++) {

            MARCS.setElementAt(0, i);
        }

        for (int row = 0; row < neighbor_bondnum_A; row++) {
            for (int column = 0; column < neighbor_bondnum_B; column++) {

                String G1A = c_bond_neighborsA.get(row * 4 + 0);
                String G2A = c_bond_neighborsA.get(row * 4 + 1);
                String G1B = c_bond_neighborsB.get(column * 4 + 0);
                String G2B = c_bond_neighborsB.get(column * 4 + 1);


                if ((G1A.compareToIgnoreCase(G1B) == 0 && G2A.compareToIgnoreCase(G2B) == 0) || (G1A.compareToIgnoreCase(G2B) == 0 && G2A.compareToIgnoreCase(G1B) == 0)) {


                    int Index_I = i_bond_neighbor_atoms_A.get(row * 3 + 0);
                    int Index_IPlus1 = i_bond_neighbor_atoms_A.get(row * 3 + 1);

                    IAtom R1_A = ac1.getAtom(Index_I);
                    IAtom R2_A = ac1.getAtom(Index_IPlus1);
                    IBond ReactantBond = ac1.getBond(R1_A, R2_A);

                    int Index_J = i_bond_neighbor_atoms_B.get(column * 3 + 0);
                    int Index_JPlus1 = i_bond_neighbor_atoms_B.get(column * 3 + 1);

                    IAtom P1_B = ac2.getAtom(Index_J);
                    IAtom P2_B = ac2.getAtom(Index_JPlus1);
                    IBond ProductBond = ac2.getBond(P1_B, P2_B);
                    if (bondTypeFlag && bondMatch(ReactantBond, ProductBond)) {
                        MARCS.setElementAt(1, row * neighbor_bondnum_B + column);

                    } else if (!bondTypeFlag) {
                        MARCS.setElementAt(1, row * neighbor_bondnum_B + column);
                    }
                }


            }
        }
        first = last = new BinaryTree(-1);
        last.equal = null;
        last.not_equal = null;

        bestarcsleft = 0;

        startsearch();
        Stack<Vector<Integer>> BESTARCS_copy = new Stack<Vector<Integer>>();


        BESTARCS_copy.addAll(BESTARCS);
        while (!BESTARCS.empty()) {
            BESTARCS.pop();
        }

        while (!BESTARCS_copy.empty()) {

            Vector<Integer> MARCS_vector = new Vector<Integer>(BESTARCS_copy.peek());
            Vector<Integer> new_MAPPING = findMcGregorMapping(MARCS_vector, mapped_atom_num, mapped_atoms, neighbor_bondnum_A, i_bond_neighbor_atoms_A, neighbor_bondnum_B, i_bond_neighbor_atoms_B);

            int new_MAPPING_size = new_MAPPING.size() / 2;
            boolean no_further_MAPPINGS = false;
            if (mapped_atom_num == new_MAPPING_size) {
                no_further_MAPPINGS = true;
            }


            int new_neighbor_numA = 0; //instead of neighbor_bondnum_A

            int new_neighbor_numB = 0; //instead of neighbor_bondnum_B

            Vector<Integer> new_i_neighborsA = new Vector<Integer>(); //instead of i_bond_neighbor_atoms_A

            Vector<Integer> new_i_neighborsB = new Vector<Integer>(); //instead of i_bond_neighbor_atoms_B

            Vector<String> new_c_neighborsA = new Vector<String>(); //instead of c_bond_neighborsA

            Vector<String> new_c_neighborsB = new Vector<String>(); //instead of c_bond_neighborsB

            new_i_neighborsA.clear();
            new_i_neighborsB.clear();
            new_c_neighborsA.clear();
            new_c_neighborsB.clear();


            //new values for set_num_A + set_num_B
            //new arrays for i_bond_setA + i_bond_setB + c_bond_setB + c_bond_setB

            set_bondnum_A = 0; //instead of set_num_A

            set_bondnum_B = 0; //instead of set_num_B

            Vector<Integer> new_i_bond_setA = new Vector<Integer>(); //instead of i_bond_setA
            Vector<Integer> new_i_bond_setB = new Vector<Integer>(); //instead of i_bond_setB
            Vector<String> new_c_bond_setA = new Vector<String>(); //instead of c_bond_setA
            Vector<String> new_c_bond_setB = new Vector<String>(); //instead of c_bond_setB
            Vector<String> c_setB_copy = new Vector<String>();
            Vector<String> c_setA_copy = new Vector<String>();

            generateCSetACopy(set_num_A, c_bond_setA, c_setA_copy);
            generateCSetBCopy(set_num_B, c_bond_setB, c_setB_copy);

            //find unmapped atoms of molecule A
            Vector<Integer> unmapped_atoms_molA = new Vector<Integer>();
            int unmapped_numA = 0;
            boolean atomA_is_unmapped = true;

            for (int a = 0; a < ac1.getAtomCount(); a++) {
                for (int b = 0; b < new_MAPPING_size; b++) {
                    if (a == new_MAPPING.elementAt(b * 2 + 0)) {
                        atomA_is_unmapped = false;
                    }

                }
                if (atomA_is_unmapped) {
                    unmapped_atoms_molA.add(a);
                    unmapped_numA++;

                }


                atomA_is_unmapped = true;
            }


            //The special signs must be transfered to the corresponding atoms of molecule B

            int SR_count = 0;
            boolean bond_considered = false;
            boolean normal_bond = true;
            for (int a = 0; a < set_num_A; a++) {

                int _elementAt_a = i_bond_setA.get(a * 3 + 0).intValue();
                for (int b = 0; b < unmapped_numA; b++) {
                    Integer unMappedAtomIndex = unmapped_atoms_molA.elementAt(b);
                    if (unMappedAtomIndex == _elementAt_a) {
                        for (int c = 0; c < new_MAPPING_size; c++) {

                            if (new_MAPPING.elementAt(c * 2 + 0).equals(i_bond_setA.elementAt(a * 3 + 1))) {

                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 0));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 1));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 2));
                                new_c_neighborsA.add(c_setA_copy.get(a * 4 + 0));
                                if (c_setA_copy.elementAt(a * 4 + 3).compareToIgnoreCase("X") == 0) {

                                    new_c_neighborsA.add(SignROW[SR_count]);
                                    new_c_neighborsA.add("X");
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    changeCharBonds(i_bond_setA.get(a * 3 + 1), SignROW[SR_count], set_num_A, i_bond_setA, c_setA_copy);
                                    int cor_atom = searchCorrespondingAtom(new_MAPPING_size, i_bond_setA.get(a * 3 + 1), 1, new_MAPPING);
                                    changeCharBonds(cor_atom, SignROW[SR_count], set_num_B, i_bond_setB, c_setB_copy);
                                    SR_count++;

                                } else {

                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    new_c_neighborsA.add("X");
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 3));

                                }

                                normal_bond = false;
                                new_neighbor_numA++;

                            }
                        }

                        if (normal_bond) {

                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 0));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 1));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 2));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 0));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 1));
                            new_c_bond_setA.add("X");
                            new_c_bond_setA.add("X");
                            set_bondnum_A++;

                        }
                        normal_bond = true;
                        bond_considered = true;
                    }
                    if (unMappedAtomIndex.equals(i_bond_setA.elementAt(a * 3 + 1))) {
                        for (int c = 0; c < new_MAPPING_size; c++) {

                            if (new_MAPPING.elementAt(c * 2 + 0).equals(i_bond_setA.elementAt(a * 3 + 0))) {

                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 0));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 1));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 2));
                                if (c_setA_copy.elementAt(a * 4 + 2).compareToIgnoreCase("X") == 0) {

                                    new_c_neighborsA.add(SignROW[SR_count]);
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 0));
                                    new_c_neighborsA.add("X");
                                    changeCharBonds(i_bond_setA.get(a * 3 + 0), SignROW[SR_count], set_num_A, i_bond_setA, c_setA_copy);
                                    int cor_atom = searchCorrespondingAtom(new_MAPPING_size, i_bond_setA.get(a * 3 + 0), 1, new_MAPPING);
                                    changeCharBonds(cor_atom, SignROW[SR_count], set_num_B, i_bond_setB, c_setB_copy);
                                    SR_count++;

                                } else {

                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 0));
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 2));
                                    new_c_neighborsA.add("X");

                                }

                                normal_bond = false;
                                new_neighbor_numA++;

                            }


                        }
                        if (normal_bond) {

                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 0));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 1));
                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 2));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 0));
                            new_c_bond_setA.add(c_setA_copy.get(a * 4 + 1));
                            new_c_bond_setA.add("X");
                            new_c_bond_setA.add("X");
                            set_bondnum_A++;

                        }

                        normal_bond = true;
                        bond_considered = true;
                    }

                    if (bond_considered) {
                        break;
                    }

                }

                bond_considered = false;
            }

            //find unmapped atoms of molecule B

            Vector<Integer> unmapped_atoms_molB = new Vector<Integer>();
            int unmapped_numB = 0;
            boolean atomB_is_unmapped = true;

            for (int a = 0; a < ac2.getAtomCount(); a++) {
                for (int b = 0; b < new_MAPPING_size; b++) {
                    if (a == new_MAPPING.elementAt(b * 2 + 1)) {
                        atomB_is_unmapped = false;
                    }

                }
                if (atomB_is_unmapped) {
                    unmapped_atoms_molB.add(a);
                    unmapped_numB++;

                }


                atomB_is_unmapped = true;
            }

            //The special signs must be transfered to the corresponding atoms of molecule A

            bond_considered = false;
            normal_bond = true;
            for (int a = 0; a < set_num_B; a++) {
                for (int b = 0; b < unmapped_numB; b++) {
                    if (unmapped_atoms_molB.elementAt(b).equals(i_bond_setB.elementAt(a * 3 + 0))) {
                        for (int c = 0; c < new_MAPPING_size; c++) {
                            if (new_MAPPING.elementAt(c * 2 + 1).equals(i_bond_setB.elementAt(a * 3 + 1))) {
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 0));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 1));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 2));
                                if (c_setB_copy.elementAt(a * 4 + 3).compareToIgnoreCase("X") == 0) {
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add(SignROW[SR_count]);
                                    new_c_neighborsB.add("X");
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    changeCharBonds(i_bond_setB.get(a * 3 + 1), SignROW[SR_count], set_num_B, i_bond_setB, c_setB_copy);
                                    int cor_atom = searchCorrespondingAtom(new_MAPPING_size, i_bond_setB.get(a * 3 + 1), 2, new_MAPPING);
                                    changeCharBonds(cor_atom, SignROW[SR_count], new_neighbor_numA, new_i_neighborsA, new_c_neighborsA);
                                    SR_count++;

                                } else {
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    new_c_neighborsB.add("X");
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 3));
                                }

                                normal_bond = false;
                                new_neighbor_numB++;

                            }


                        }
                        if (normal_bond) {
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 0));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 1));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 2));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 0));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 1));
                            new_c_bond_setB.add("X");
                            new_c_bond_setB.add("X");
                            set_bondnum_B++;

                        }


                        normal_bond = true;
                        bond_considered = true;
                    }
                    if (unmapped_atoms_molB.elementAt(b).equals(i_bond_setB.elementAt(a * 3 + 1))) {
                        for (int c = 0; c < new_MAPPING_size; c++) {

                            if (new_MAPPING.elementAt(c * 2 + 1).equals(i_bond_setB.elementAt(a * 3 + 0))) {

                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 0));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 1));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 2));

                                if (c_setB_copy.elementAt(a * 4 + 2).compareToIgnoreCase("X") == 0) {
                                    new_c_neighborsB.add(SignROW[SR_count]);
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add("X");
                                    changeCharBonds(i_bond_setB.get(a * 3 + 0), SignROW[SR_count], set_num_B, i_bond_setB, c_setB_copy);
                                    int cor_atom = searchCorrespondingAtom(new_MAPPING_size, i_bond_setB.get(a * 3 + 0), 2, new_MAPPING);
                                    changeCharBonds(cor_atom, SignROW[SR_count], new_neighbor_numA, new_i_neighborsA, new_c_neighborsA);
                                    SR_count++;

                                } else {

                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 2));
                                    new_c_neighborsB.add("X");
                                }

                                normal_bond = false;
                                new_neighbor_numB++;

                            }


                        }

                        if (normal_bond) {

                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 0));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 1));
                            new_i_bond_setB.add(i_bond_setB.get(a * 3 + 2));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 0));
                            new_c_bond_setB.add(c_setB_copy.get(a * 4 + 1));
                            new_c_bond_setB.add("X");
                            new_c_bond_setB.add("X");
                            set_bondnum_B++;

                        }


                        normal_bond = true;
                        bond_considered = true;
                    }

                    if (bond_considered) {
                        break;
                    }

                }
                bond_considered = false;
            }


//             System.out.println("Mapped Atoms before Iterator2: " + mapped_atoms);
            iterator(no_further_MAPPINGS, new_MAPPING_size, new_MAPPING, new_neighbor_numA, new_neighbor_numB, new_i_neighborsA, new_i_neighborsB, new_c_neighborsA, new_c_neighborsB,
                    set_bondnum_A, set_bondnum_B, new_i_bond_setA, new_i_bond_setB, new_c_bond_setA, new_c_bond_setB);
            BESTARCS_copy.pop();
//            System.out.println("Schleife beendet in iterator!!!!");
        }

        //}
        //System.out.println("In the iterator Termination");
        //System.out.println("============+++++++++==============");

        //System.out.println("Mapped Atoms before iterator Over: " + mapped_atoms);
        return 0;
    }

    private int generateCTab1Copy() throws IOException {
        IAtomContainer reactant = ac1;
        for (int a = 0; a < reactant.getBondCount(); a++) {
            String AtomI = reactant.getBond(a).getAtom(0).getSymbol();
            String AtomJ = reactant.getBond(a).getAtom(1).getSymbol();
            c_tab1_copy.addElement(AtomI);
            c_tab1_copy.addElement(AtomJ);
            c_tab1_copy.addElement("X");
            c_tab1_copy.addElement("X");
        }
        return 0;
    }

    private int generateCTab2Copy() throws IOException {
        IAtomContainer product = ac2;
        for (int a = 0; a < product.getBondCount(); a++) {
            String AtomI = product.getBond(a).getAtom(0).getSymbol();
            String AtomJ = product.getBond(a).getAtom(1).getSymbol();
            c_tab2_copy.addElement(AtomI);
            c_tab2_copy.addElement(AtomJ);
            c_tab2_copy.addElement("X");
            c_tab2_copy.addElement("X");
        }

        return 0;

    }

    private int changeCharBonds(int corresponding_atom, String new_symbol, int neighbor_bondnum, Vector<Integer> i_bond_neighbors, Vector<String> c_bond_neighbors) {

        for (int a = 0; a < neighbor_bondnum; a++) {
            if ((i_bond_neighbors.elementAt(a * 3 + 0) == (corresponding_atom)) && (c_bond_neighbors.elementAt(a * 4 + 2).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 2, c_bond_neighbors.get(a * 4 + 0));
                c_bond_neighbors.set(a * 4 + 0, new_symbol);
            }

            if ((i_bond_neighbors.elementAt(a * 3 + 1) == (corresponding_atom)) && (c_bond_neighbors.elementAt(a * 4 + 3).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 3, c_bond_neighbors.get(a * 4 + 1));
                c_bond_neighbors.set(a * 4 + 1, new_symbol);
            }

        }

        return 0;
    }

    public boolean bondMatch(IBond ReactantBond, IBond ProductBond) {
        boolean Flag = false;
        int ReactantBondType = ReactantBond.getOrder().ordinal();
        int ProductBondType = ProductBond.getOrder().ordinal();
        if (bondTypeFlag) {
            if ((ReactantBond.getFlag(CDKConstants.ISAROMATIC) == ProductBond.getFlag(CDKConstants.ISAROMATIC)) && (ReactantBondType == ProductBondType)) {
                Flag = true;
            }
            if (ReactantBond.getFlag(CDKConstants.ISAROMATIC) && ProductBond.getFlag(CDKConstants.ISAROMATIC)) {
                Flag = true;
            }
        }
        return Flag;
    }

    private Vector<Integer> findMcGregorMapping(Vector<Integer> MARCS_vector_org, int mapped_atoms_num, List<Integer> current_MAPPING_org, int bondnum_A, Vector<Integer> i_bonds_A_org, int bondnum_B, Vector<Integer> i_bonds_B_org) {

        Vector<Integer> MARCS_vector = new Vector<Integer>(MARCS_vector_org);
        Vector<Integer> current_MAPPING = new Vector<Integer>(current_MAPPING_org);
        Vector<Integer> i_bonds_A = new Vector<Integer>(i_bonds_A_org);
        Vector<Integer> i_bonds_B = new Vector<Integer>(i_bonds_B_org);

        Vector<Integer> additional_mapping = new Vector<Integer>();



        for (int x = 0; x < bondnum_A; x++) {
            for (int y = 0; y < bondnum_B; y++) {


                if (MARCS_vector.elementAt(x * bondnum_B + y) == 1) {


                    int Atom1_moleculeA = i_bonds_A.get(x * 3 + 0);
                    int Atom2_moleculeA = i_bonds_A.get(x * 3 + 1);
                    int Atom1_moleculeB = i_bonds_B.get(y * 3 + 0);
                    int Atom2_moleculeB = i_bonds_B.get(y * 3 + 1);

                    IAtom R1_A = ac1.getAtom(Atom1_moleculeA);
                    IAtom R2_A = ac1.getAtom(Atom2_moleculeA);
                    IBond ReactantBond = ac1.getBond(R1_A, R2_A);

                    IAtom P1_B = ac2.getAtom(Atom1_moleculeB);
                    IAtom P2_B = ac2.getAtom(Atom2_moleculeB);
                    IBond ProductBond = ac2.getBond(P1_B, P2_B);


//                    //Bond Order Check Introduced by Asad

                    boolean bMatch = bondMatch(ReactantBond, ProductBond);

                    for (int z = 0; z < mapped_atoms_num; z++) {

                        int Mapped_Atom_1 = current_MAPPING.elementAt(z * 2 + 0);
                        int Mapped_Atom_2 = current_MAPPING.elementAt(z * 2 + 1);

                        if ((Mapped_Atom_1 == Atom1_moleculeA) && (Mapped_Atom_2 == Atom1_moleculeB)) {
                            if (bondTypeFlag && bMatch) {
                                additional_mapping.add(Atom2_moleculeA);
                                additional_mapping.add(Atom2_moleculeB);
                            } else if (!bondTypeFlag) {
                                additional_mapping.add(Atom2_moleculeA);
                                additional_mapping.add(Atom2_moleculeB);
                            }
                        }

                        if ((Mapped_Atom_1 == Atom1_moleculeA) && (Mapped_Atom_2 == Atom2_moleculeB)) {
                            if (bondTypeFlag && bMatch) {
                                additional_mapping.add(Atom2_moleculeA);
                                additional_mapping.add(Atom1_moleculeB);
                            } else if (!bondTypeFlag) {
                                additional_mapping.add(Atom2_moleculeA);
                                additional_mapping.add(Atom1_moleculeB);
                            }
                        }

                        if ((Mapped_Atom_1 == Atom2_moleculeA) && (Mapped_Atom_2 == Atom1_moleculeB)) {
                            if (bondTypeFlag && bMatch) {
                                additional_mapping.add(Atom1_moleculeA);
                                additional_mapping.add(Atom2_moleculeB);
                            } else if (!bondTypeFlag) {
                                additional_mapping.add(Atom1_moleculeA);
                                additional_mapping.add(Atom2_moleculeB);
                            }
                        }

                        if ((Mapped_Atom_1 == Atom2_moleculeA) && (Mapped_Atom_2 == Atom2_moleculeB)) {
                            if (bondTypeFlag && bMatch) {
                                additional_mapping.add(Atom1_moleculeA);
                                additional_mapping.add(Atom1_moleculeB);
                            } else if (!bondTypeFlag) {
                                additional_mapping.add(Atom1_moleculeA);
                                additional_mapping.add(Atom1_moleculeB);
                            }
                        }

                    }
                }
            }
        }


        int additional_mapping_size = additional_mapping.size();

        //add McGregorBondTypeInSensitive mapping to the Clique mapping
        for (int a = 0; a < additional_mapping_size; a = a + 2) {
            current_MAPPING.add(additional_mapping.get(a + 0));
            current_MAPPING.add(additional_mapping.get(a + 1));
        }
        //remove recurring mappings from current_MAPPING

        Vector<Integer> unique_MAPPING = removeRecurringMappings(current_MAPPING);

        return unique_MAPPING;
    }

    //Function compaires a structure array with itself. Sometimes a mapping occurs several times within the array.
//The function eliminates these recurring mappings. Function is called in function best_solution.
//The function is called by itself as long as the last list element is processed.
    private Vector<Integer> removeRecurringMappings(Vector<Integer> atom_mapping) {

        //System.out.println("Mapped Atoms removeRecurringMappings: " + atom_mapping);

        //Vector<Integer> atom_mapping = new Vector<Integer>(atom_mapping_org);

        boolean exist = true;
        Vector<Integer> temp_map = new Vector<Integer>();
        int temp_counter = 0;
        int atom_mapping_size = atom_mapping.size();
        for (int x = 0; x < atom_mapping_size; x = x + 2) {
            int atom = atom_mapping.get(x);
            for (int y = x + 2; y < atom_mapping_size; y = y + 2) {
                if (atom == atom_mapping.elementAt(y)) {
                    exist = false;
                }

            }
            if (exist == true) {
                temp_map.add(atom_mapping.get(x + 0));
                temp_map.add(atom_mapping.get(x + 1));
                temp_counter = temp_counter + 2;
            }

            exist = true;
        }


        return temp_map;
    }

    private void partsearch(int xstart, int ystart, Vector<Integer> TEMPMARCS_ORG) {


        int x = xstart;
        int y = ystart;

        Vector<Integer> TEMPMARCS = new Vector<Integer>(TEMPMARCS_ORG);

        if (TEMPMARCS.elementAt(xstart * nNum_globalB + ystart) == 1) {

            removeRedundantArcs(xstart, ystart, TEMPMARCS);
            int arcsleft = 0;

            for (int a = 0; a < nNum_globalA; a++) {
                for (int b = 0; b < nNum_globalB; b++) {

                    if (TEMPMARCS.elementAt(a * nNum_globalB + b) == (1)) {
                        arcsleft++;
                    }

                }
            }

            //test Bestarcsleft and skip rest if needed
            if (arcsleft >= bestarcsleft) {
                do {
                    y++;
                    if (y == nNum_globalB) {
                        y = 0;
                        x++;

                    }
                } while ((x < nNum_globalA) && (TEMPMARCS.elementAt(x * nNum_globalB + y) != 1)); //Correction by ASAD set value minus 1
                if (x < nNum_globalA) {

                    partsearch(x, y, TEMPMARCS);
                    TEMPMARCS.setElementAt(0, x * nNum_globalB + y);
                    partsearch(x, y, TEMPMARCS);

                } else {
                    if (arcsleft > bestarcsleft) {
                        removeTreeStructure(first);
                        first = last = new BinaryTree(-1);
                        last.equal = null;
                        last.not_equal = null;

                        while (!BESTARCS.empty()) {
                            BESTARCS.pop();
                        }

                    }
                    bestarcsleft = arcsleft;

                    if (checkMARCS(TEMPMARCS)) {
                        BESTARCS.push(TEMPMARCS);
                    }

                }
            }
        } else {
            do {
                y++;
                if (y == nNum_globalB) {
                    y = 0;
                    x++;

                }


            } while ((x < nNum_globalA) && (TEMPMARCS.elementAt(x * nNum_globalB + y) != 1)); //Correction by ASAD set value minus 1

            if (x < nNum_globalA) {

                partsearch(x, y, TEMPMARCS);
                TEMPMARCS.setElementAt(0, x * nNum_globalB + y);
                partsearch(x, y, TEMPMARCS);
            } else {
                int arcsleft = 0;
                for (int a = 0; a < nNum_globalA; a++) {
                    for (int b = 0; b < nNum_globalB; b++) {
                        if (TEMPMARCS.elementAt(a * nNum_globalB + b) == 1) {
                            arcsleft++;
                        }

                    }
                }
                if (arcsleft >= bestarcsleft) {
                    if (arcsleft > bestarcsleft) {
                        removeTreeStructure(first);
                        first = last = new BinaryTree(-1);
                        last.equal = null;
                        last.not_equal = null;
                        while (!BESTARCS.empty()) {
                            BESTARCS.pop();
                        }

                    }
                    bestarcsleft = arcsleft;

                    if (checkMARCS(TEMPMARCS)) {
                        BESTARCS.push(TEMPMARCS);
                    }

                }
            }
        }
    }

    //The function is called in function partsearch. The function is given a temporary matrix and a position (row/column)
//within this matrix. First the function sets all entries to zero, which can be exlcuded in respect to the current
//atom by atom matching. After this the function replaces all entries in the same row and column of the current
//position by zeros. Only the entry of the current position is set to one.
//Return value "count_arcsleft" counts the number of arcs, which are still in the matrix.
    private void removeRedundantArcs(int row, int column, Vector<Integer> MARCS) {

        //System.err.print("Betrachte: " + c_globalA.get(row*2+0) + c_globalA.get(row*2+1));
        //System.err.println( " und " + c_globalA.get(column*2+0) + c_globalA.get(column*2+1) );
        int G1_atom = i_globalA.get(row * 3 + 0);
        int G2_atom = i_globalA.get(row * 3 + 1);
        int G3_atom = i_globalB.get(column * 3 + 0);
        int G4_atom = i_globalB.get(column * 3 + 1);

        for (int x = 0; x < nNum_globalA; x++) {
            int row_atom1 = i_globalA.get(x * 3 + 0);
            int row_atom2 = i_globalA.get(x * 3 + 1);

            for (int y = 0; y < nNum_globalB; y++) {
                int column_atom3 = i_globalB.get(y * 3 + 0);
                int column_atom4 = i_globalB.get(y * 3 + 1);

                if (((G1_atom == row_atom1) || (G1_atom == row_atom2)) && (!(((column_atom3 == G3_atom) || (column_atom4 == G3_atom)) || ((column_atom3 == G4_atom) || (column_atom4 == G4_atom))))) {

                    MARCS.setElementAt(0, x * nNum_globalB + y);
                }

                if (((G2_atom == row_atom1) || (G2_atom == row_atom2)) && (!(((column_atom3 == G3_atom) || (column_atom4 == G3_atom)) || ((column_atom3 == G4_atom) || (column_atom4 == G4_atom))))) {
                    MARCS.setElementAt(0, x * nNum_globalB + y);
                }

                if (((G3_atom == column_atom3) || (G3_atom == column_atom4)) && (!(((row_atom1 == G1_atom) || (row_atom2 == G1_atom)) || ((row_atom1 == G2_atom) || (row_atom2 == G2_atom))))) {
                    MARCS.setElementAt(0, x * nNum_globalB + y);
                }

                if (((G4_atom == column_atom3) || (G4_atom == column_atom4)) && (!(((row_atom1 == G1_atom) || (row_atom2 == G1_atom)) || ((row_atom1 == G2_atom) || (row_atom2 == G2_atom))))) {
                    MARCS.setElementAt(0, x * nNum_globalB + y);
                }

            }
        }

        for (int v = 0; v < nNum_globalA; v++) {
            MARCS.set(v * nNum_globalB + column, 0);
        }

        for (int w = 0; w < nNum_globalB; w++) {
            MARCS.set(row * nNum_globalB + w, 0);
        }

        MARCS.set(row * nNum_globalB + column, 1);
        //System.err.println("MARCS: " + MARCS);
    }

    /* Modified function call by ASAD in Java have to check
     *
     */
    private int removeTreeStructure(BinaryTree cur_struc) {

        BinaryTree equal_struc = cur_struc.equal;
        BinaryTree not_equal_struc = cur_struc.not_equal;
        cur_struc = null;


        if (equal_struc != null) {
            removeTreeStructure(equal_struc);
        }

        if (not_equal_struc != null) {
            removeTreeStructure(not_equal_struc);
        }

        return 0;
    }

    //The function is called in function partsearch. The function is given z temporary matrix.
//The function checks whether the temporary matrix is already found by calling the function
//"verifyNodes". If the matrix already exists the function returns false which means that
//the matrix will not be stored. Otherwise the function returns true which means that the
//matrix will be stored in function partsearch.
    //private boolean checkMARCS(Vector<Integer> MARCS) {
    private boolean checkMARCS(Vector<Integer> MARCS_T) {


        Vector<Integer> posnum_list = new Vector<Integer>();
        posnum_list.setSize(nNum_globalA * nNum_globalA); /*TO DO ASAD Initi by 0;*/

        for (int i = 0; i < posnum_list.size(); i++) {
            posnum_list.set(i, 0);
        }

        int y = 0;
        int count_entries = 0;
        for (int x = 0; x < (nNum_globalA * nNum_globalB); x++) {
            if (MARCS_T.elementAt(x) == 1) {
                posnum_list.setElementAt(x, y++);
                count_entries++;
            }
        }
        boolean flag = false;

        verifyNodes(posnum_list, first, 0, count_entries);
        if (new_matrix) {
            flag = true;
        }

        return flag;

    }

    /* Modified function call by ASAD in Java have to check
     *
     */
    private boolean verifyNodes(Vector<Integer> matrix, BinaryTree cur_struc, int x, int field_length) {

//        Vector<Integer> matrix = new Vector<Integer>(matrix_org);

        if ((matrix.elementAt(x) == cur_struc.getValue()) && (x < field_length)) {
            if (cur_struc.equal != null) {
                new_matrix = false;
                verifyNodes(matrix, cur_struc.equal, x + 1, field_length);
            }

        }
        if (matrix.elementAt(x) != cur_struc.getValue()) {
            if (cur_struc.not_equal != null) {
                verifyNodes(matrix, cur_struc.not_equal, x, field_length);
            }

            if (cur_struc.not_equal == null) {
                cur_struc.not_equal = new BinaryTree(matrix.elementAt(x));
                cur_struc.not_equal.not_equal = null;
                int y = 0;


                BinaryTree last_one = cur_struc.not_equal;

                while ((y + x + 1) < field_length) {
                    last_one.equal = new BinaryTree(matrix.elementAt(y + x + 1));
                    last_one = last_one.equal;
                    last_one.not_equal = null;
                    y++;

                }


                last_one.equal = null;
                new_matrix = true;
            }

        }
        return true;
    }

    private int searchCorrespondingAtom(int mapped_atoms_size, int atom_from_other_molecule, int molecule, Vector<Integer> mapped_atoms_org) {


        Vector<Integer> mapped_atoms = new Vector<Integer>(mapped_atoms_org);

        int corresponding_atom = 0;
        for (int a = 0; a < mapped_atoms_size; a++) {
            if (molecule == 1) {
                if (mapped_atoms.elementAt(a * 2 + 0).intValue() == atom_from_other_molecule) {

                    corresponding_atom = mapped_atoms.get(a * 2 + 1);
                }

            }
            if (molecule == 2) {
                if (mapped_atoms.elementAt(a * 2 + 1).intValue() == atom_from_other_molecule) {
                    corresponding_atom = mapped_atoms.get(a * 2 + 0);
                }

            }
        }
        return corresponding_atom;
    }

    private int generateCSetBCopy(int bond_number, Vector<String> c_setB, Vector<String> c_setB_copy) {

        for (int a = 0; a < bond_number; a++) {
            c_setB_copy.add(c_setB.get(a * 4 + 0));
            c_setB_copy.add(c_setB.get(a * 4 + 1));
            c_setB_copy.add("X");
            c_setB_copy.add("X");
        }

        return 0;
    }

    private int generateCSetACopy(int bond_number, Vector<String> c_setA, Vector<String> c_setA_copy) {

        for (int a = 0; a < bond_number; a++) {
            c_setA_copy.add(c_setA.get(a * 4 + 0));
            c_setA_copy.add(c_setA.get(a * 4 + 1));
            c_setA_copy.add("X");
            c_setA_copy.add("X");
        }

        return 0;
    }

    private void startsearch() {
        Vector<Integer> FIXARCS = new Vector<Integer>(nNum_globalA * nNum_globalB);//  Initialize FIXARCS with 0

        FIXARCS.setSize(nNum_globalA * nNum_globalB);

        for (int i = 0; i < nNum_globalA * nNum_globalB; i++) {

            FIXARCS.set(i, 0);
        }

        int x = 0;
        int y = 0;

        while ((x < nNum_globalA) && (MARCS.elementAt(x * nNum_globalB + y) != 1)) {
            y++;
            if (y == nNum_globalB) {
                y = 0;
                x++;

            }

        }

        if (x == nNum_globalA) {
            y = nNum_globalB - 1;
            x = x - 1;
        }
        if (MARCS.elementAt(x * nNum_globalB + y) == 0) {
            partsearch(x, y, MARCS);
        }
        if (MARCS.elementAt(x * nNum_globalB + y) != 0) {
            partsearch(x, y, MARCS);
            MARCS.set(x * nNum_globalB + y, 0);
            partsearch(x, y, MARCS);
        }

    }

    /**
     *
     * @return
     */
    public List<List<Integer>> getMappings() {

        return this.mappings;
    }

    /**
     *
     * @return
     */
    public int getMCSSize() {

        return this.globalMCSSize;
    }
}
