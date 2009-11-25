/* Copyright (C) 2009 Syed Asad Rahman {asad@ebi.ac.uk}
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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smsd.global.BondType;
import org.openscience.cdk.smsd.helper.BinaryTree;

/**
 * @cdk.module smsd
 */
public class McGregor {

    private IAtomContainer source = null;
    private IAtomContainer target = null;
    private BinaryTree last = null;
    private BinaryTree first = null;
    private Stack<List<Integer>> bestARCS = null;
    private int bestarcsleft = 0;
    private int globalMCSSize = 0;
    private List<List<Integer>> mappings = null;
    /*This should be more or equal to all the atom types*/
    private static final String[] SignArray = {"$1", "$2", "$3", "$4", "$5", "$6", "$7", "$8", "$9", "$10", "$11", "$12",
        "$13", "$15", "$16", "$17", "$18", "$19", "$20", "$21", "$22", "$23", "$24",
        "$25", "$26", "$27", "$28", "$29", "$30", "$31", "$32", "$33", "$34", "$35", "$36",
        "$37", "$38", "$39", "$40", "$41", "$42", "$43", "$44", "$45", "$46",
        "$47", "$48", "$49", "$50", "$51", "$52", "$53", "$54", "$55"
    };
    protected boolean newMatrix = false;
    /**
     * Bond sensitive flag
     *
     */
    protected boolean bondTypeFlag = BondType.getInstance().getBondSensitiveFlag();

    /**
     * Creates a new instance of McGregor
     * @param source
     * @param target
     * @param _mappings
     */
    public McGregor(IAtomContainer source, IAtomContainer target, List<List<Integer>> _mappings) {


        this.source = source;
        this.target = target;
        this.mappings = _mappings;
        bestarcsleft = 0;

        if (_mappings.isEmpty()) {
            this.globalMCSSize = 0;
        } else {
            this.globalMCSSize = _mappings.get(0).size();
        }

        bestARCS = new Stack<List<Integer>>();
        newMatrix = false;
    }

    /**
     *
     * @param best_Mapping_size
     * @param present_Mapping
     * @throws IOException
     */
    public void startMcGregorIteration(int best_Mapping_size, Map<Integer, Integer> present_Mapping) throws IOException {

        int neighborBondnumA = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
        int setBondNumA = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
        int neighborBondNumB = 0; //number of remaining molecule B bonds after the clique search, which are neighbors of the MCS_1
        int setBondNumB = 0; //number of remaining molecule B bonds after the clique search, which aren't neighbors

        List<Integer> i_bond_neighborsA = new ArrayList<Integer>();
        List<Integer> i_bond_setA = new ArrayList<Integer>();
        List<String> c_bond_neighborsA = new ArrayList<String>();
        List<String> c_bond_setA = new ArrayList<String>();

        List<Integer> i_bond_neighborsB = new ArrayList<Integer>();
        List<Integer> i_bond_setB = new ArrayList<Integer>();
        List<String> c_bond_neighborsB = new ArrayList<String>();
        List<String> c_bond_setB = new ArrayList<String>();

        this.globalMCSSize = (best_Mapping_size / 2);
        int counter = 0;

        List<String> c_tab1_copy = McGregorChecks.generateCTabCopy(source);
        List<String> c_tab2_copy = McGregorChecks.generateCTabCopy(target);


        //find mapped atoms of both molecules and store these in mapped_atoms
        List<Integer> mapped_atoms = new ArrayList<Integer>();
//        System.out.println("\nMapped Atoms");
        for (Map.Entry<Integer, Integer> map : present_Mapping.entrySet()) {
//            System.out.println("i:" + map.getKey() + " j:" + map.getValue());
            mapped_atoms.add(map.getKey());
            mapped_atoms.add(map.getValue());
        }
        int mapping_size = present_Mapping.size();

        //find unmapped atoms of molecule A
        List<Integer> unmapped_atoms_molA = new ArrayList<Integer>();

        int unmapped_numA = 0;
        boolean atomA_is_unmapped = true;

        for (int a = 0; a < source.getAtomCount(); a++) {
            //Atomic list are only numbers from 1 to atom_number1

            for (Integer key : present_Mapping.keySet()) {
                if (key == a) {
                    atomA_is_unmapped = false;
                }
            }


            if (atomA_is_unmapped) {
                unmapped_atoms_molA.add(unmapped_numA, a);
                unmapped_numA++;
            }
            atomA_is_unmapped = true;
        }

        QueryProcessor queryProcess = new QueryProcessor
                (neighborBondNumB,
                c_bond_neighborsA,
                i_bond_neighborsA,
                c_bond_setA,
                i_bond_setA,
                setBondNumA);
        queryProcess.process(source,
                target,
                c_tab1_copy,
                c_tab2_copy,
                SignArray,
                unmapped_atoms_molA,
                mapping_size,
                mapped_atoms,
                counter);

        i_bond_setA = queryProcess.getIBondSetA();
        c_bond_setA = queryProcess.getCBondSetA();
        setBondNumA = queryProcess.getBondNumA();
        neighborBondnumA = queryProcess.getNeighborBondNumA();
        i_bond_neighborsA = queryProcess.getIBondNeighboursA();
        c_bond_neighborsA = queryProcess.getCBondNeighborsA();

        //find unmapped atoms of molecule B
        List<Integer> unmapped_atoms_molB = new ArrayList<Integer>();
        int unmapped_numB = 0;
        boolean atomB_is_unmapped = true;

        for (int a = 0; a < target.getAtomCount(); a++) {
            for (Integer value : present_Mapping.values()) {

                if (a == value) {
                    atomB_is_unmapped = false;
                }
            }
            if (atomB_is_unmapped) {
                unmapped_atoms_molB.add(unmapped_numB, a);
                unmapped_numB++;
            }
            atomB_is_unmapped = true;
        }

//        System.out.println("unmapped_atoms_molB: " + unmapped_atoms_molB.size());

        //Extract bonds which are related with unmapped atoms of molecule B.
        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
        //new special signs -> the result are two vectors: cBondNeighborsA and int_bonds_molB, which contain those
        //bonds of molecule B, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule A


        TargetProcessor targetProcess = new TargetProcessor(
                c_bond_neighborsB,
                i_bond_setB,
                c_bond_setB,
                neighborBondNumB,
                setBondNumB,
                neighborBondnumA,
                i_bond_neighborsA,
                c_bond_neighborsA);

        targetProcess.process(target, unmapped_atoms_molB,
                mapping_size,
                i_bond_neighborsB,
                mapped_atoms,
                counter,
                c_tab1_copy,
                c_tab2_copy,
                SignArray);

        setBondNumB = targetProcess.getBondNumB();
        neighborBondNumB = targetProcess.getNeighborBondNumB();
        c_bond_neighborsB = targetProcess.getCBondNeighborsB();
        i_bond_setB = targetProcess.getIBondSetB();
        c_bond_setB = targetProcess.getCBondSetB();

        boolean dummy = false;

        iterator(dummy, present_Mapping.size(), mapped_atoms, neighborBondnumA, neighborBondNumB, i_bond_neighborsA, i_bond_neighborsB, c_bond_neighborsA, c_bond_neighborsB, setBondNumA, setBondNumB, i_bond_setA, i_bond_setB, c_bond_setA, c_bond_setB);

    }

    /**
     *
     * @param best_Mapping_size
     * @param clique_vector
     * @param comp_graph_nodes
     * @throws IOException
     */
    public void startMcGregorIteration(int best_Mapping_size, List<Integer> clique_vector, List<Integer> comp_graph_nodes) throws IOException {


        int neighborBondnumA = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
        int setBondNumA = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
        int neighborBondNumB = 0; //number of remaining molecule B bonds after the clique search, which are neighbors of the MCS_1
        int setBondNumB = 0; //number of remaining molecule B bonds after the clique search, which aren't neighbors

        List<Integer> i_bond_neighborsA = new ArrayList<Integer>();
        List<Integer> i_bond_setA = new ArrayList<Integer>();
        List<String> c_bond_neighborsA = new ArrayList<String>();
        List<String> c_bond_setA = new ArrayList<String>();

        List<Integer> i_bond_neighborsB = new ArrayList<Integer>();
        List<Integer> i_bond_setB = new ArrayList<Integer>();
        List<String> c_bond_neighborsB = new ArrayList<String>();
        List<String> c_bond_setB = new ArrayList<String>();

        this.globalMCSSize = (best_Mapping_size / 2);
        int counter = 0;

        List<String> c_tab1_copy = McGregorChecks.generateCTabCopy(source);
        List<String> c_tab2_copy = McGregorChecks.generateCTabCopy(target);


        //find mapped atoms of both molecules and store these in mapped_atoms
        List<Integer> mapped_atoms = new ArrayList<Integer>();

        int mapped_atom_number = 0;
//

        i_bond_neighborsA = new ArrayList<Integer>();

        int clique_siz = clique_vector.size();
        int vec_size = comp_graph_nodes.size();


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
        List<Integer> unmapped_atoms_molA = new ArrayList<Integer>();

        int unmapped_numA = 0;
        boolean atomA_is_unmapped = true;

//        System.out.println("Mapped Atoms: " + mapped_atoms);

        for (int a = 0; a < source.getAtomCount(); a++) {
            //Atomic list are only numbers from 1 to atom_number1

            for (int b = 0; b < clique_siz; b++) {
                //the number of nodes == number of assigned pairs
                if (mapped_atoms.get(b * 2) == a) {
                    atomA_is_unmapped = false;
                }
            }
            if (atomA_is_unmapped == true) {
                unmapped_atoms_molA.add(unmapped_numA++, a);
            }
            atomA_is_unmapped = true;
        }
        //Extract bonds which are related with unmapped atoms of molecule A.
        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
        //new special signs -> the result are two vectors: cBondNeighborsA and int_bonds_molA, which contain those
        //bonds of molecule A, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule B

       QueryProcessor queryProcess = new QueryProcessor
                (neighborBondNumB,
                c_bond_neighborsA,
                i_bond_neighborsA,
                c_bond_setA,
                i_bond_setA,
                setBondNumA);
        queryProcess.process(source,
                target,
                c_tab1_copy,
                c_tab2_copy,
                SignArray,
                unmapped_atoms_molA,
                clique_siz,
                mapped_atoms,
                counter);

        i_bond_setA = queryProcess.getIBondSetA();
        c_bond_setA = queryProcess.getCBondSetA();
        setBondNumA = queryProcess.getBondNumA();
        neighborBondnumA = queryProcess.getNeighborBondNumA();
        i_bond_neighborsA = queryProcess.getIBondNeighboursA();
        c_bond_neighborsA = queryProcess.getCBondNeighborsA();
        //find unmapped atoms of molecule B
        List<Integer> unmapped_atoms_molB = new ArrayList<Integer>();
        int unmapped_numB = 0;
        boolean atomB_is_unmapped = true;

//        System.out.println("neighborBondnumA After:" + neighborBondnumA);

        for (int a = 0; a < target.getAtomCount(); a++) {
            for (int b = 0; b < clique_siz; b++) {
                if (a == mapped_atoms.get(b * 2 + 1)) {
                    atomB_is_unmapped = false;
                }
            }
            if (atomB_is_unmapped == true) {
                unmapped_atoms_molB.add(unmapped_numB++, a);
            }
            atomB_is_unmapped = true;
        }


        //Extract bonds which are related with unmapped atoms of molecule B.
        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
        //new special signs -> the result are two vectors: cBondNeighborsA and int_bonds_molB, which contain those
        //bonds of molecule B, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule A



        TargetProcessor targetProcess = new TargetProcessor(
                c_bond_neighborsB,
                i_bond_setB,
                c_bond_setB,
                neighborBondNumB,
                setBondNumB,
                neighborBondnumA,
                i_bond_neighborsA,
                c_bond_neighborsA);


        targetProcess.process(target,
                unmapped_atoms_molB,
                clique_siz,
                i_bond_neighborsB,
                mapped_atoms,
                counter,
                c_tab1_copy,
                c_tab2_copy,
                SignArray);

        setBondNumB = targetProcess.getBondNumB();
        neighborBondNumB = targetProcess.getNeighborBondNumB();
        c_bond_neighborsB = targetProcess.getCBondNeighborsB();
        i_bond_setB = targetProcess.getIBondSetB();
        c_bond_setB = targetProcess.getCBondSetB();

        boolean dummy = false;

        iterator(dummy, mapped_atom_number, mapped_atoms, neighborBondnumA, neighborBondNumB, i_bond_neighborsA, i_bond_neighborsB, c_bond_neighborsA, c_bond_neighborsB, setBondNumA, setBondNumB, i_bond_setA, i_bond_setB, c_bond_setA, c_bond_setB);

    }

    private int iterator(boolean mappingCheckFlag,
            int mappedAtomCount,
            List<Integer> mapped_atoms,
            int neighborBondNumA,
            int neighborBondNumB,
            List<Integer> iBondNeighborAtomsA,
            List<Integer> iBondNeighborAtomsB,
            List<String> cBondNeighborsA,
            List<String> cBondNeighborsB,
            int setNumA, int setNumB,
            List<Integer> iBondSetA,
            List<Integer> iBondSetB,
            List<String> cBondSetA,
            List<String> cBondSetB) throws IOException {

        boolean moreMappingPossible = McGregorChecks.isFurtherMappingPossible(
                source,
                target,
                neighborBondNumA,
                neighborBondNumB,
                iBondNeighborAtomsA,
                iBondNeighborAtomsB,
                cBondNeighborsA,
                cBondNeighborsB);

        if (neighborBondNumA == 0 || neighborBondNumB == 0 || mappingCheckFlag || !moreMappingPossible) {
            try {
                if (mappedAtomCount >= globalMCSSize) {
//                    System.out.println("Hello-1");
                    if (mappedAtomCount > globalMCSSize) {
//                        System.out.println("Hello-2");
                        this.globalMCSSize = mappedAtomCount;
//                        System.out.println("best_MAPPING_size: " + globalMCSSize);
                        mappings.clear();
                    }
                    mappings.add(mapped_atoms);
//                    System.out.println("mappings " + mappings);
                }
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            return 0;
        }
        List<Integer> modifiedARCS = new ArrayList<Integer>();
        int size = neighborBondNumA * neighborBondNumB;
        for (int i = 0; i < size; i++) {
            modifiedARCS.add(i, 0);
        }

        modifiedARCS = McGregorChecks.setArcs(
                source,
                target,
                neighborBondNumA,
                neighborBondNumB,
                iBondNeighborAtomsA,
                iBondNeighborAtomsB,
                cBondNeighborsA,
                cBondNeighborsB,
                modifiedARCS);

        first = last = new BinaryTree(-1);
        last.equal = null;
        last.not_equal = null;
        bestarcsleft = 0;

        startsearch(neighborBondNumA, neighborBondNumB, iBondNeighborAtomsA, iBondNeighborAtomsB, modifiedARCS);
        Stack<List<Integer>> BESTARCS_copy = new Stack<List<Integer>>();

        BESTARCS_copy.addAll(bestARCS);
        while (!bestARCS.empty()) {
            bestARCS.pop();
        }
        while (!BESTARCS_copy.empty()) {

            List<Integer> MARCS_vector = new ArrayList<Integer>(BESTARCS_copy.peek());
            List<Integer> localMAPPING = findMcGregorMapping(MARCS_vector, mappedAtomCount, mapped_atoms, neighborBondNumA, iBondNeighborAtomsA, neighborBondNumB, iBondNeighborAtomsB);

            int localMappingSize = localMAPPING.size() / 2;
            boolean no_further_MAPPINGS = false;
            if (mappedAtomCount == localMappingSize) {
                no_further_MAPPINGS = true;
            }


            int newNeighborNumA = 0; //instead of neighborBondnumA
            int newNeighborNumB = 0; //instead of neighborBondNumB

            List<Integer> localIBondNeighborsA = new ArrayList<Integer>(); //instead of iBondNeighborAtomsA
            List<Integer> localIBondNeighborsB = new ArrayList<Integer>(); //instead of iBondNeighborAtomsB
            List<String> localCBondNeighborsA = new ArrayList<String>(); //instead of cBondNeighborsA
            List<String> localCBondNeighborsB = new ArrayList<String>(); //instead of cBondNeighborsB

            //new values for setNumA + setNumB
            //new arrays for iBondSetA + iBondSetB + cBondSetB + cBondSetB

            int localBondNumA = 0; //instead of setNumA
            int localBondNumB = 0; //instead of setNumB

            List<Integer> localIBondSetA = new ArrayList<Integer>(); //instead of iBondSetA
            List<Integer> localIBondSetB = new ArrayList<Integer>(); //instead of iBondSetB
            List<String> localCBondSetA = new ArrayList<String>(); //instead of cBondSetA
            List<String> localCBondSetB = new ArrayList<String>(); //instead of cBondSetB

            List<String> cTab1Copy = McGregorChecks.generateCSetCopy(setNumA, cBondSetA);
            List<String> cTab2Copy = McGregorChecks.generateCSetCopy(setNumB, cBondSetB);

            //find unmapped atoms of molecule A
            List<Integer> unmapped_atoms_molA = new ArrayList<Integer>();
            int unmapped_numA = 0;
            boolean atomA_is_unmapped = true;

            for (int a = 0; a < source.getAtomCount(); a++) {
                for (int b = 0; b < localMappingSize; b++) {
                    if (a == localMAPPING.get(b * 2 + 0)) {
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

            int counter = 0;

            QueryProcessor queryProcess = new QueryProcessor(
                    newNeighborNumA,
                    localCBondNeighborsA,
                    localIBondNeighborsA,
                    localCBondSetA,
                    localIBondSetA,
                    localBondNumA);

            queryProcess.process(
                    setNumA,
                    iBondSetA,
                    iBondSetB,
                    unmapped_atoms_molA,
                    localMappingSize,
                    localMAPPING,
                    cTab1Copy,
                    cTab2Copy,
                    SignArray,
                    counter,
                    setNumB);

               
            localIBondSetA = queryProcess.getIBondSetA();
            localCBondSetA = queryProcess.getCBondSetA();
            localBondNumA = queryProcess.getBondNumA();
            newNeighborNumA = queryProcess.getNeighborBondNumA();
            localIBondNeighborsA = queryProcess.getIBondNeighboursA();
            localCBondNeighborsA = queryProcess.getCBondNeighborsA();
            //find unmapped atoms of molecule B

            List<Integer> unmapped_atoms_molB = new ArrayList<Integer>();
            int unmapped_numB = 0;
            boolean atomB_is_unmapped = true;

            for (int a = 0; a < target.getAtomCount(); a++) {
                for (int b = 0; b < localMappingSize; b++) {
                    if (a == localMAPPING.get(b * 2 + 1)) {
                        atomB_is_unmapped = false;
                    }
                }
                if (atomB_is_unmapped) {
                    unmapped_atoms_molB.add(a);
                    unmapped_numB++;
                }
                atomB_is_unmapped = true;
            }


//          The special signs must be transfered to the corresponding atoms of molecule A

            TargetProcessor targetP = new TargetProcessor(
                    localCBondNeighborsB,
                    iBondSetB,
                    localCBondSetB,
                    newNeighborNumB,
                    localBondNumB,
                    newNeighborNumA,
                    localIBondNeighborsA,
                    localCBondNeighborsA);

            targetP.process(setNumB, unmapped_atoms_molA,
                    localMappingSize,
                    localIBondNeighborsB,
                    localMAPPING,
                    counter,
                    cTab1Copy,
                    cTab2Copy,
                    SignArray);

            localBondNumB = targetP.getBondNumB();
            newNeighborNumB = targetP.getNeighborBondNumB();
            localCBondNeighborsB = targetP.getCBondNeighborsB();
            localIBondSetB = targetP.getIBondSetB();
            localCBondSetB = targetP.getCBondSetB();

//             System.out.println("Mapped Atoms before Iterator2: " + mapped_atoms);
            iterator(no_further_MAPPINGS, localMappingSize, localMAPPING, newNeighborNumA, newNeighborNumB, localIBondNeighborsA, localIBondNeighborsB, localCBondNeighborsA, localCBondNeighborsB,
                    localBondNumA, localBondNumB, localIBondSetA, localIBondSetB, localCBondSetA, localCBondSetB);
            BESTARCS_copy.pop();
        }
        //System.out.println("Mapped Atoms before iterator Over: " + mapped_atoms);
        return 0;
    }

    private List<Integer> findMcGregorMapping(List<Integer> MARCS, int mapped_atoms_num, List<Integer> current_MAPPING_org, int bondnum_A, List<Integer> i_bonds_A_org, int bondnum_B, List<Integer> i_bonds_B_org) {

        List<Integer> currentMapping = new ArrayList<Integer>(current_MAPPING_org);
        List<Integer> i_bonds_A = new ArrayList<Integer>(i_bonds_A_org);
        List<Integer> i_bonds_B = new ArrayList<Integer>(i_bonds_B_org);
        List<Integer> additional_mapping = new ArrayList<Integer>();



        for (int x = 0; x < bondnum_A; x++) {
            for (int y = 0; y < bondnum_B; y++) {


                if (MARCS.get(x * bondnum_B + y) == 1) {


                    int Atom1_moleculeA = i_bonds_A.get(x * 3 + 0);
                    int Atom2_moleculeA = i_bonds_A.get(x * 3 + 1);
                    int Atom1_moleculeB = i_bonds_B.get(y * 3 + 0);
                    int Atom2_moleculeB = i_bonds_B.get(y * 3 + 1);

                    IAtom R1_A = source.getAtom(Atom1_moleculeA);
                    IAtom R2_A = source.getAtom(Atom2_moleculeA);
                    IBond ReactantBond = source.getBond(R1_A, R2_A);

                    IAtom P1_B = target.getAtom(Atom1_moleculeB);
                    IAtom P2_B = target.getAtom(Atom2_moleculeB);
                    IBond ProductBond = target.getBond(P1_B, P2_B);

//                  Bond Order Check Introduced by Asad

                    boolean bMatch = McGregorChecks.bondMatch(ReactantBond, ProductBond);

                    if ((bondTypeFlag && bMatch) || (!bondTypeFlag)) {

                        for (int z = 0; z < mapped_atoms_num; z++) {

                            int Mapped_Atom_1 = currentMapping.get(z * 2 + 0);
                            int Mapped_Atom_2 = currentMapping.get(z * 2 + 1);

                            if ((Mapped_Atom_1 == Atom1_moleculeA) && (Mapped_Atom_2 == Atom1_moleculeB)) {
                                additional_mapping.add(Atom2_moleculeA);
                                additional_mapping.add(Atom2_moleculeB);
                            }

                            if ((Mapped_Atom_1 == Atom1_moleculeA) && (Mapped_Atom_2 == Atom2_moleculeB)) {
                                additional_mapping.add(Atom2_moleculeA);
                                additional_mapping.add(Atom1_moleculeB);
                            }

                            if ((Mapped_Atom_1 == Atom2_moleculeA) && (Mapped_Atom_2 == Atom1_moleculeB)) {
                                additional_mapping.add(Atom1_moleculeA);
                                additional_mapping.add(Atom2_moleculeB);
                            }

                            if ((Mapped_Atom_1 == Atom2_moleculeA) && (Mapped_Atom_2 == Atom2_moleculeB)) {
                                additional_mapping.add(Atom1_moleculeA);
                                additional_mapping.add(Atom1_moleculeB);
                            }
                        }//for loop
                    }
                }
            }
        }


        int additionalMappingSize = additional_mapping.size();

        //add McGregorBondTypeInSensitive mapping to the Clique mapping
        for (int a = 0; a < additionalMappingSize; a = a + 2) {
            currentMapping.add(additional_mapping.get(a + 0));
            currentMapping.add(additional_mapping.get(a + 1));
        }
//        remove recurring mappings from currentMapping

        List<Integer> unique_MAPPING = McGregorChecks.removeRecurringMappings(currentMapping);

        return unique_MAPPING;
    }

    private List<Integer> partsearch(int xstart, int ystart, List<Integer> TempArcs_Org, int nNumGlobalA, int nNumGlobalB, List<Integer> iGlobalA, List<Integer> iGlobalB) {
        int xIndex = xstart;
        int yIndex = ystart;

        List<Integer> TempArcs = new ArrayList<Integer>(TempArcs_Org);

        if (TempArcs.get(xstart * nNumGlobalB + ystart) == 1) {

            removeRedundantArcs(xstart, ystart, TempArcs, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
            int arcsleft = 0;

            for (int a = 0; a < nNumGlobalA; a++) {
                for (int b = 0; b < nNumGlobalB; b++) {

                    if (TempArcs.get(a * nNumGlobalB + b) == (1)) {
                        arcsleft++;
                    }
                }
            }

            //test Bestarcsleft and skip rest if needed
            if (arcsleft >= bestarcsleft) {
                do {
                    yIndex++;
                    if (yIndex == nNumGlobalB) {
                        yIndex = 0;
                        xIndex++;

                    }
                } while ((xIndex < nNumGlobalA) && (TempArcs.get(xIndex * nNumGlobalB + yIndex) != 1)); //Correction by ASAD set value minus 1
                if (xIndex < nNumGlobalA) {

                    partsearch(xIndex, yIndex, TempArcs, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
                    TempArcs.set(xIndex * nNumGlobalB + yIndex, 0);
                    partsearch(xIndex, yIndex, TempArcs, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);

                } else {
                    if (arcsleft > bestarcsleft) {
                        McGregorChecks.removeTreeStructure(first);
                        first = last = new BinaryTree(-1);
                        last.equal = null;
                        last.not_equal = null;

                        while (!bestARCS.empty()) {
                            bestARCS.pop();
                        }
                    }
                    bestarcsleft = arcsleft;

                    if (checkMARCS(TempArcs, nNumGlobalA, nNumGlobalB)) {
                        bestARCS.push(TempArcs);
                    }
                }
            }
        } else {
            do {
                yIndex++;
                if (yIndex == nNumGlobalB) {
                    yIndex = 0;
                    xIndex++;
                }

            } while ((xIndex < nNumGlobalA) && (TempArcs.get(xIndex * nNumGlobalB + yIndex) != 1)); //Correction by ASAD set value minus 1

            if (xIndex < nNumGlobalA) {

                partsearch(xIndex, yIndex, TempArcs, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
                TempArcs.set(xIndex * nNumGlobalB + yIndex, 0);
                partsearch(xIndex, yIndex, TempArcs, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
            } else {
                int arcsleft = 0;
                for (int a = 0; a <
                        nNumGlobalA; a++) {
                    for (int b = 0; b <
                            nNumGlobalB; b++) {
                        if (TempArcs.get(a * nNumGlobalB + b) == 1) {
                            arcsleft++;
                        }

                    }
                }
                if (arcsleft >= bestarcsleft) {
                    if (arcsleft > bestarcsleft) {
                        McGregorChecks.removeTreeStructure(first);
                        first = last = new BinaryTree(-1);
                        last.equal = null;
                        last.not_equal = null;
                        while (!bestARCS.empty()) {
                            bestARCS.pop();
                        }
                    }
                    bestarcsleft = arcsleft;

                    if (checkMARCS(TempArcs, nNumGlobalA, nNumGlobalB)) {
                        bestARCS.push(TempArcs);
                    }

                }
            }
        }
        return TempArcs;
    }

//The function is called in function partsearch. The function is given a temporary matrix and a position (row/column)
//within this matrix. First the function sets all entries to zero, which can be exlcuded in respect to the current
//atom by atom matching. After this the function replaces all entries in the same row and column of the current
//position by zeros. Only the entry of the current position is set to one.
//Return value "count_arcsleft" counts the number of arcs, which are still in the matrix.
    private void removeRedundantArcs(int row, int column, List<Integer> MARCS, int nNumGlobalA, int nNumGlobalB,
            List<Integer> iGlobalA, List<Integer> iGlobalB) {

        int G1_atom = iGlobalA.get(row * 3 + 0);
        int G2_atom = iGlobalA.get(row * 3 + 1);
        int G3_atom = iGlobalB.get(column * 3 + 0);
        int G4_atom = iGlobalB.get(column * 3 + 1);

        for (int x = 0; x < nNumGlobalA; x++) {
            int row_atom1 = iGlobalA.get(x * 3 + 0);
            int row_atom2 = iGlobalA.get(x * 3 + 1);

            for (int y = 0; y < nNumGlobalB; y++) {
                int column_atom3 = iGlobalB.get(y * 3 + 0);
                int column_atom4 = iGlobalB.get(y * 3 + 1);

                if (McGregorChecks.cases(G1_atom, G2_atom, G3_atom, G4_atom, row_atom1, row_atom2, column_atom3, column_atom4)) {
                    MARCS.set(x * nNumGlobalB + y, 0);
                }

            }
        }

        for (int v = 0; v < nNumGlobalA; v++) {
            MARCS.set(v * nNumGlobalB + column, 0);
        }

        for (int w = 0; w < nNumGlobalB; w++) {
            MARCS.set(row * nNumGlobalB + w, 0);
        }

        MARCS.set(row * nNumGlobalB + column, 1);
    }

//The function is called in function partsearch. The function is given z temporary matrix.
//The function checks whether the temporary matrix is already found by calling the function
//"verifyNodes". If the matrix already exists the function returns false which means that
//the matrix will not be stored. Otherwise the function returns true which means that the
//matrix will be stored in function partsearch.
    private boolean checkMARCS(List<Integer> MARCS_T, int nNumGlobalA, int nNumGlobalB) {

        int size = nNumGlobalA * nNumGlobalA;
        List<Integer> posnum_list = new ArrayList<Integer>(size);

        for (int i = 0; i < posnum_list.size(); i++) {
            posnum_list.add(i, 0);
        }

        int yCounter = 0;
        int count_entries = 0;
        for (int x = 0; x < (nNumGlobalA * nNumGlobalB); x++) {
            if (MARCS_T.get(x) == 1) {
                posnum_list.add(yCounter++, x);
                count_entries++;
            }
        }
        boolean flag = false;

        verifyNodes(posnum_list, first, 0, count_entries);
        if (newMatrix) {
            flag = true;
        }

        return flag;

    }

    private boolean verifyNodes(List<Integer> matrix, BinaryTree currentStructure, int index, int fieldLength) {

        if (((matrix.get(index) == currentStructure.getValue()) && (index < fieldLength)) && (currentStructure.equal != null)) {
            newMatrix = false;
            verifyNodes(matrix, currentStructure.equal, index + 1, fieldLength);
        }
        if (matrix.get(index) != currentStructure.getValue()) {
            if (currentStructure.not_equal != null) {
                verifyNodes(matrix, currentStructure.not_equal, index, fieldLength);
            }

            if (currentStructure.not_equal == null) {
                currentStructure.not_equal = new BinaryTree(matrix.get(index));
                currentStructure.not_equal.not_equal = null;
                int yIndex = 0;


                BinaryTree last_one = currentStructure.not_equal;

                while ((yIndex + index + 1) < fieldLength) {
                    last_one.equal = new BinaryTree(matrix.get(yIndex + index + 1));
                    last_one = last_one.equal;
                    last_one.not_equal = null;
                    yIndex++;

                }
                last_one.equal = null;
                newMatrix = true;
            }

        }
        return true;
    }

    private void startsearch(int nNumGlobalA, int nNumGlobalB, List<Integer> iGlobalA, List<Integer> iGlobalB, List<Integer> modifiedARCS) {
        int size = nNumGlobalA * nNumGlobalB;
        List<Integer> FIXARCS = new ArrayList<Integer>(size);//  Initialize FIXARCS with 0
        for (int i = 0; i < size; i++) {
            FIXARCS.add(i, 0);
        }

        int xIndex = 0;
        int yIndex = 0;

        while ((xIndex < nNumGlobalA) && (modifiedARCS.get(xIndex * nNumGlobalB + yIndex) != 1)) {
            yIndex++;
            if (yIndex == nNumGlobalB) {
                yIndex = 0;
                xIndex++;
            }
        }

        if (xIndex == nNumGlobalA) {
            yIndex = nNumGlobalB - 1;
            xIndex = xIndex - 1;
        }

        if (modifiedARCS.get(xIndex * nNumGlobalB + yIndex) == 0) {
            partsearch(xIndex, yIndex, modifiedARCS, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
        }

        if (modifiedARCS.get(xIndex * nNumGlobalB + yIndex) != 0) {
            partsearch(xIndex, yIndex, modifiedARCS, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
            modifiedARCS.set(xIndex * nNumGlobalB + yIndex, 0);
            partsearch(xIndex, yIndex, modifiedARCS, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
        }

    }

    /**
     * @return
     */
    public List<List<Integer>> getMappings() {

        return this.mappings;
    }

    /**
     * @return
     */
    public int getMCSSize() {

        return this.globalMCSSize;
    }
}
