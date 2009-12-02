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
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smsd.global.BondType;
import org.openscience.cdk.smsd.helper.BinaryTree;

/**
 * @cdk.module smsd
 */
//public class McGregor {
//
//    private IAtomContainer source = null;
//    private IAtomContainer target = null;
//    private BinaryTree last = null;
//    private BinaryTree first = null;
//    private Stack<List<Integer>> bestARCS = null;
//    private int bestarcsleft = 0;
//    private int globalMCSSize = 0;
//    private List<List<Integer>> mappings = null;
//    /*This should be more or equal to all the atom types*/
//    private static final String[] SignArray = {"$1", "$2", "$3", "$4", "$5", "$6", "$7", "$8", "$9", "$10", "$11", "$12",
//        "$13", "$15", "$16", "$17", "$18", "$19", "$20", "$21", "$22", "$23", "$24",
//        "$25", "$26", "$27", "$28", "$29", "$30", "$31", "$32", "$33", "$34", "$35", "$36",
//        "$37", "$38", "$39", "$40", "$41", "$42", "$43", "$44", "$45", "$46",
//        "$47", "$48", "$49", "$50", "$51", "$52", "$53", "$54", "$55"
//    };
//    protected boolean newMatrix = false;
//    /**
//     * Bond sensitive flag
//     *
//     */
//    protected boolean bondTypeFlag = BondType.getInstance().getBondSensitiveFlag();
//
//    /**
//     * Creates a new instance of McGregor
//     * @param source
//     * @param target
//     * @param _mappings
//     */
//    public McGregor(IAtomContainer source, IAtomContainer target, List<List<Integer>> _mappings) {
//
//
//        this.source = source;
//        this.target = target;
//        this.mappings = _mappings;
//        bestarcsleft = 0;
//
//        if (_mappings.isEmpty()) {
//            this.globalMCSSize = 0;
//        } else {
//            this.globalMCSSize = _mappings.get(0).size();
//        }
//
//        bestARCS = new Stack<List<Integer>>();
//        newMatrix = false;
//    }
//
//    /**
//     *
//     * @param best_Mapping_size
//     * @param present_Mapping
//     * @throws IOException
//     */
//    public void startMcGregorIteration(int best_Mapping_size, Map<Integer, Integer> present_Mapping) throws IOException {
//
//        int neighborBondnumA = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
//        int setBondNumA = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
//        int neighborBondNumB = 0; //number of remaining molecule B bonds after the clique search, which are neighbors of the MCS_1
//        int setBondNumB = 0; //number of remaining molecule B bonds after the clique search, which aren't neighbors
//
//        List<Integer> i_bond_neighborsA = new ArrayList<Integer>();
//        List<Integer> i_bond_setA = new ArrayList<Integer>();
//        List<String> c_bond_neighborsA = new ArrayList<String>();
//        List<String> c_bond_setA = new ArrayList<String>();
//
//        List<Integer> i_bond_neighborsB = new ArrayList<Integer>();
//        List<Integer> i_bond_setB = new ArrayList<Integer>();
//        List<String> c_bond_neighborsB = new ArrayList<String>();
//        List<String> c_bond_setB = new ArrayList<String>();
//
//        this.globalMCSSize = (best_Mapping_size / 2);
//        int counter = 0;
//
//        List<String> c_tab1_copy = McGregorChecks.generateCTabCopy(source);
//        List<String> c_tab2_copy = McGregorChecks.generateCTabCopy(target);
//
//
//        //find mapped atoms of both molecules and store these in mapped_atoms
//        List<Integer> mapped_atoms = new ArrayList<Integer>();
////        System.out.println("\nMapped Atoms");
//        for (Map.Entry<Integer, Integer> map : present_Mapping.entrySet()) {
////            System.out.println("i:" + map.getKey() + " j:" + map.getValue());
//            mapped_atoms.add(map.getKey());
//            mapped_atoms.add(map.getValue());
//        }
//        int mapping_size = present_Mapping.size();
//
//        //find unmapped atoms of molecule A
//        List<Integer> unmapped_atoms_molA = new ArrayList<Integer>();
//
//        int unmapped_numA = 0;
//        boolean atomA_is_unmapped = true;
//
//        for (int a = 0; a < source.getAtomCount(); a++) {
//            //Atomic list are only numbers from 1 to atom_number1
//
//            for (Integer key : present_Mapping.keySet()) {
//                if (key == a) {
//                    atomA_is_unmapped = false;
//                }
//            }
//
//
//            if (atomA_is_unmapped) {
//                unmapped_atoms_molA.add(unmapped_numA, a);
//                unmapped_numA++;
//            }
//            atomA_is_unmapped = true;
//        }
//
//        QueryProcessor queryProcess = new QueryProcessor
//                (neighborBondNumB,
//                c_bond_neighborsA,
//                i_bond_neighborsA,
//                c_bond_setA,
//                i_bond_setA,
//                setBondNumA);
//        queryProcess.process(source,
//                target,
//                c_tab1_copy,
//                c_tab2_copy,
//                SignArray,
//                unmapped_atoms_molA,
//                mapping_size,
//                mapped_atoms,
//                counter);
//
//        i_bond_setA = queryProcess.getIBondSetA();
//        c_bond_setA = queryProcess.getCBondSetA();
//        setBondNumA = queryProcess.getBondNumA();
//        neighborBondnumA = queryProcess.getNeighborBondNumA();
//        i_bond_neighborsA = queryProcess.getIBondNeighboursA();
//        c_bond_neighborsA = queryProcess.getCBondNeighborsA();
//
//        //find unmapped atoms of molecule B
//        List<Integer> unmapped_atoms_molB = new ArrayList<Integer>();
//        int unmapped_numB = 0;
//        boolean atomB_is_unmapped = true;
//
//        for (int a = 0; a < target.getAtomCount(); a++) {
//            for (Integer value : present_Mapping.values()) {
//
//                if (a == value) {
//                    atomB_is_unmapped = false;
//                }
//            }
//            if (atomB_is_unmapped) {
//                unmapped_atoms_molB.add(unmapped_numB, a);
//                unmapped_numB++;
//            }
//            atomB_is_unmapped = true;
//        }
//
////        System.out.println("unmapped_atoms_molB: " + unmapped_atoms_molB.size());
//
//        //Extract bonds which are related with unmapped atoms of molecule B.
//        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
//        //new special signs -> the result are two vectors: cBondNeighborsA and int_bonds_molB, which contain those
//        //bonds of molecule B, which are relevant for the McGregorBondTypeInSensitive algorithm.
//        //The special signs must be transfered to the corresponding atoms of molecule A
//
//
//        TargetProcessor targetProcess = new TargetProcessor(
//                c_bond_neighborsB,
//                i_bond_setB,
//                c_bond_setB,
//                neighborBondNumB,
//                setBondNumB,
//                neighborBondnumA,
//                i_bond_neighborsA,
//                c_bond_neighborsA);
//
//        targetProcess.process(target, unmapped_atoms_molB,
//                mapping_size,
//                i_bond_neighborsB,
//                mapped_atoms,
//                counter,
//                c_tab1_copy,
//                c_tab2_copy,
//                SignArray);
//
//        setBondNumB = targetProcess.getBondNumB();
//        neighborBondNumB = targetProcess.getNeighborBondNumB();
//        c_bond_neighborsB = targetProcess.getCBondNeighborsB();
//        i_bond_setB = targetProcess.getIBondSetB();
//        c_bond_setB = targetProcess.getCBondSetB();
//
//        boolean dummy = false;
//
//        iterator(dummy, present_Mapping.size(), mapped_atoms, neighborBondnumA, neighborBondNumB, i_bond_neighborsA, i_bond_neighborsB, c_bond_neighborsA, c_bond_neighborsB, setBondNumA, setBondNumB, i_bond_setA, i_bond_setB, c_bond_setA, c_bond_setB);
//
//    }
//
//    /**
//     *
//     * @param best_Mapping_size
//     * @param clique_vector
//     * @param comp_graph_nodes
//     * @throws IOException
//     */
//    public void startMcGregorIteration(int best_Mapping_size, List<Integer> clique_vector, List<Integer> comp_graph_nodes) throws IOException {
//
//
//        int neighborBondnumA = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
//        int setBondNumA = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
//        int neighborBondNumB = 0; //number of remaining molecule B bonds after the clique search, which are neighbors of the MCS_1
//        int setBondNumB = 0; //number of remaining molecule B bonds after the clique search, which aren't neighbors
//
//        List<Integer> i_bond_neighborsA = new ArrayList<Integer>();
//        List<Integer> i_bond_setA = new ArrayList<Integer>();
//        List<String> c_bond_neighborsA = new ArrayList<String>();
//        List<String> c_bond_setA = new ArrayList<String>();
//
//        List<Integer> i_bond_neighborsB = new ArrayList<Integer>();
//        List<Integer> i_bond_setB = new ArrayList<Integer>();
//        List<String> c_bond_neighborsB = new ArrayList<String>();
//        List<String> c_bond_setB = new ArrayList<String>();
//
//        this.globalMCSSize = (best_Mapping_size / 2);
//        int counter = 0;
//
//        List<String> c_tab1_copy = McGregorChecks.generateCTabCopy(source);
//        List<String> c_tab2_copy = McGregorChecks.generateCTabCopy(target);
//
//
//        //find mapped atoms of both molecules and store these in mapped_atoms
//        List<Integer> mapped_atoms = new ArrayList<Integer>();
//
//        int mapped_atom_number = 0;
////
//
//        i_bond_neighborsA = new ArrayList<Integer>();
//
//        int clique_siz = clique_vector.size();
//        int vec_size = comp_graph_nodes.size();
//
//
//        int cliqueNumber = 0;
//        for (int a = 0; a < clique_siz; a++) {
//            //go through all clique nodes
//            cliqueNumber = clique_vector.get(a);
//            for (int b = 0; b < vec_size; b = b + 3) {
//                //go through all nodes in the compatibility graph
//                if (cliqueNumber == comp_graph_nodes.get(b + 2)) {
//                    mapped_atoms.add(comp_graph_nodes.get(b));
//                    mapped_atoms.add(comp_graph_nodes.get(b + 1));
//                    mapped_atom_number++;
//                }
//            }
//        }
//
//
//        //find unmapped atoms of molecule A
//        List<Integer> unmapped_atoms_molA = new ArrayList<Integer>();
//
//        int unmapped_numA = 0;
//        boolean atomA_is_unmapped = true;
//
////        System.out.println("Mapped Atoms: " + mapped_atoms);
//
//        for (int a = 0; a < source.getAtomCount(); a++) {
//            //Atomic list are only numbers from 1 to atom_number1
//
//            for (int b = 0; b < clique_siz; b++) {
//                //the number of nodes == number of assigned pairs
//                if (mapped_atoms.get(b * 2) == a) {
//                    atomA_is_unmapped = false;
//                }
//            }
//            if (atomA_is_unmapped == true) {
//                unmapped_atoms_molA.add(unmapped_numA++, a);
//            }
//            atomA_is_unmapped = true;
//        }
//        //Extract bonds which are related with unmapped atoms of molecule A.
//        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
//        //new special signs -> the result are two vectors: cBondNeighborsA and int_bonds_molA, which contain those
//        //bonds of molecule A, which are relevant for the McGregorBondTypeInSensitive algorithm.
//        //The special signs must be transfered to the corresponding atoms of molecule B
//
//       QueryProcessor queryProcess = new QueryProcessor
//                (neighborBondNumB,
//                c_bond_neighborsA,
//                i_bond_neighborsA,
//                c_bond_setA,
//                i_bond_setA,
//                setBondNumA);
//        queryProcess.process(source,
//                target,
//                c_tab1_copy,
//                c_tab2_copy,
//                SignArray,
//                unmapped_atoms_molA,
//                clique_siz,
//                mapped_atoms,
//                counter);
//
//        i_bond_setA = queryProcess.getIBondSetA();
//        c_bond_setA = queryProcess.getCBondSetA();
//        setBondNumA = queryProcess.getBondNumA();
//        neighborBondnumA = queryProcess.getNeighborBondNumA();
//        i_bond_neighborsA = queryProcess.getIBondNeighboursA();
//        c_bond_neighborsA = queryProcess.getCBondNeighborsA();
//        //find unmapped atoms of molecule B
//        List<Integer> unmapped_atoms_molB = new ArrayList<Integer>();
//        int unmapped_numB = 0;
//        boolean atomB_is_unmapped = true;
//
////        System.out.println("neighborBondnumA After:" + neighborBondnumA);
//
//        for (int a = 0; a < target.getAtomCount(); a++) {
//            for (int b = 0; b < clique_siz; b++) {
//                if (a == mapped_atoms.get(b * 2 + 1)) {
//                    atomB_is_unmapped = false;
//                }
//            }
//            if (atomB_is_unmapped == true) {
//                unmapped_atoms_molB.add(unmapped_numB++, a);
//            }
//            atomB_is_unmapped = true;
//        }
//
//
//        //Extract bonds which are related with unmapped atoms of molecule B.
//        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
//        //new special signs -> the result are two vectors: cBondNeighborsA and int_bonds_molB, which contain those
//        //bonds of molecule B, which are relevant for the McGregorBondTypeInSensitive algorithm.
//        //The special signs must be transfered to the corresponding atoms of molecule A
//
//
//
//        TargetProcessor targetProcess = new TargetProcessor(
//                c_bond_neighborsB,
//                i_bond_setB,
//                c_bond_setB,
//                neighborBondNumB,
//                setBondNumB,
//                neighborBondnumA,
//                i_bond_neighborsA,
//                c_bond_neighborsA);
//
//
//        targetProcess.process(target,
//                unmapped_atoms_molB,
//                clique_siz,
//                i_bond_neighborsB,
//                mapped_atoms,
//                counter,
//                c_tab1_copy,
//                c_tab2_copy,
//                SignArray);
//
//        setBondNumB = targetProcess.getBondNumB();
//        neighborBondNumB = targetProcess.getNeighborBondNumB();
//        c_bond_neighborsB = targetProcess.getCBondNeighborsB();
//        i_bond_setB = targetProcess.getIBondSetB();
//        c_bond_setB = targetProcess.getCBondSetB();
//
//        boolean dummy = false;
//
//        iterator(dummy, mapped_atom_number, mapped_atoms, neighborBondnumA, neighborBondNumB, i_bond_neighborsA, i_bond_neighborsB, c_bond_neighborsA, c_bond_neighborsB, setBondNumA, setBondNumB, i_bond_setA, i_bond_setB, c_bond_setA, c_bond_setB);
//
//    }
//
//    private int iterator(boolean mappingCheckFlag,
//            int mappedAtomCount,
//            List<Integer> mapped_atoms,
//            int neighborBondNumA,
//            int neighborBondNumB,
//            List<Integer> iBondNeighborAtomsA,
//            List<Integer> iBondNeighborAtomsB,
//            List<String> cBondNeighborsA,
//            List<String> cBondNeighborsB,
//            int setNumA, int setNumB,
//            List<Integer> iBondSetA,
//            List<Integer> iBondSetB,
//            List<String> cBondSetA,
//            List<String> cBondSetB) throws IOException {
//
//        boolean moreMappingPossible = McGregorChecks.isFurtherMappingPossible(
//                source,
//                target,
//                neighborBondNumA,
//                neighborBondNumB,
//                iBondNeighborAtomsA,
//                iBondNeighborAtomsB,
//                cBondNeighborsA,
//                cBondNeighborsB);
//
//        if (neighborBondNumA == 0 || neighborBondNumB == 0 || mappingCheckFlag || !moreMappingPossible) {
//            try {
//                if (mappedAtomCount >= globalMCSSize) {
////                    System.out.println("Hello-1");
//                    if (mappedAtomCount > globalMCSSize) {
////                        System.out.println("Hello-2");
//                        this.globalMCSSize = mappedAtomCount;
////                        System.out.println("best_MAPPING_size: " + globalMCSSize);
//                        mappings.clear();
//                    }
//                    mappings.add(mapped_atoms);
////                    System.out.println("mappings " + mappings);
//                }
//            } catch (Exception ex) {
//                ex.printStackTrace();
//            }
//            return 0;
//        }
//        List<Integer> modifiedARCS = new ArrayList<Integer>();
//        int size = neighborBondNumA * neighborBondNumB;
//        for (int i = 0; i < size; i++) {
//            modifiedARCS.add(i, 0);
//        }
//
//        modifiedARCS = McGregorChecks.setArcs(
//                source,
//                target,
//                neighborBondNumA,
//                neighborBondNumB,
//                iBondNeighborAtomsA,
//                iBondNeighborAtomsB,
//                cBondNeighborsA,
//                cBondNeighborsB,
//                modifiedARCS);
//
//        first = last = new BinaryTree(-1);
//        last.equal = null;
//        last.not_equal = null;
//        bestarcsleft = 0;
//
//        startsearch(neighborBondNumA, neighborBondNumB, iBondNeighborAtomsA, iBondNeighborAtomsB, modifiedARCS);
//        Stack<List<Integer>> BESTARCS_copy = new Stack<List<Integer>>();
//
//        BESTARCS_copy.addAll(bestARCS);
//        while (!bestARCS.empty()) {
//            bestARCS.pop();
//        }
//        while (!BESTARCS_copy.empty()) {
//
//            List<Integer> MARCS_vector = new ArrayList<Integer>(BESTARCS_copy.peek());
//            List<Integer> localMAPPING = findMcGregorMapping(MARCS_vector, mappedAtomCount, mapped_atoms, neighborBondNumA, iBondNeighborAtomsA, neighborBondNumB, iBondNeighborAtomsB);
//
//            int localMappingSize = localMAPPING.size() / 2;
//            boolean no_further_MAPPINGS = false;
//            if (mappedAtomCount == localMappingSize) {
//                no_further_MAPPINGS = true;
//            }
//
//
//            int newNeighborNumA = 0; //instead of neighborBondnumA
//            int newNeighborNumB = 0; //instead of neighborBondNumB
//
//            List<Integer> localIBondNeighborsA = new ArrayList<Integer>(); //instead of iBondNeighborAtomsA
//            List<Integer> localIBondNeighborsB = new ArrayList<Integer>(); //instead of iBondNeighborAtomsB
//            List<String> localCBondNeighborsA = new ArrayList<String>(); //instead of cBondNeighborsA
//            List<String> localCBondNeighborsB = new ArrayList<String>(); //instead of cBondNeighborsB
//
//            //new values for setNumA + setNumB
//            //new arrays for iBondSetA + iBondSetB + cBondSetB + cBondSetB
//
//            int localBondNumA = 0; //instead of setNumA
//            int localBondNumB = 0; //instead of setNumB
//
//            List<Integer> localIBondSetA = new ArrayList<Integer>(); //instead of iBondSetA
//            List<Integer> localIBondSetB = new ArrayList<Integer>(); //instead of iBondSetB
//            List<String> localCBondSetA = new ArrayList<String>(); //instead of cBondSetA
//            List<String> localCBondSetB = new ArrayList<String>(); //instead of cBondSetB
//
//            List<String> cTab1Copy = McGregorChecks.generateCSetCopy(setNumA, cBondSetA);
//            List<String> cTab2Copy = McGregorChecks.generateCSetCopy(setNumB, cBondSetB);
//
//            //find unmapped atoms of molecule A
//            List<Integer> unmapped_atoms_molA = new ArrayList<Integer>();
//            int unmapped_numA = 0;
//            boolean atomA_is_unmapped = true;
//
//            for (int a = 0; a < source.getAtomCount(); a++) {
//                for (int b = 0; b < localMappingSize; b++) {
//                    if (a == localMAPPING.get(b * 2 + 0)) {
//                        atomA_is_unmapped = false;
//                    }
//                }
//                if (atomA_is_unmapped) {
//                    unmapped_atoms_molA.add(a);
//                    unmapped_numA++;
//                }
//                atomA_is_unmapped = true;
//            }
//
//
//            //The special signs must be transfered to the corresponding atoms of molecule B
//
//            int counter = 0;
//
//            QueryProcessor queryProcess = new QueryProcessor(
//                    newNeighborNumA,
//                    localCBondNeighborsA,
//                    localIBondNeighborsA,
//                    localCBondSetA,
//                    localIBondSetA,
//                    localBondNumA);
//
//            queryProcess.process(
//                    setNumA,
//                    iBondSetA,
//                    iBondSetB,
//                    unmapped_atoms_molA,
//                    localMappingSize,
//                    localMAPPING,
//                    cTab1Copy,
//                    cTab2Copy,
//                    SignArray,
//                    counter,
//                    setNumB);
//
//
//            localIBondSetA = queryProcess.getIBondSetA();
//            localCBondSetA = queryProcess.getCBondSetA();
//            localBondNumA = queryProcess.getBondNumA();
//            newNeighborNumA = queryProcess.getNeighborBondNumA();
//            localIBondNeighborsA = queryProcess.getIBondNeighboursA();
//            localCBondNeighborsA = queryProcess.getCBondNeighborsA();
//            //find unmapped atoms of molecule B
//
//            List<Integer> unmapped_atoms_molB = new ArrayList<Integer>();
//            int unmapped_numB = 0;
//            boolean atomB_is_unmapped = true;
//
//            for (int a = 0; a < target.getAtomCount(); a++) {
//                for (int b = 0; b < localMappingSize; b++) {
//                    if (a == localMAPPING.get(b * 2 + 1)) {
//                        atomB_is_unmapped = false;
//                    }
//                }
//                if (atomB_is_unmapped) {
//                    unmapped_atoms_molB.add(a);
//                    unmapped_numB++;
//                }
//                atomB_is_unmapped = true;
//            }
//
//
////          The special signs must be transfered to the corresponding atoms of molecule A
//
//            TargetProcessor targetP = new TargetProcessor(
//                    localCBondNeighborsB,
//                    iBondSetB,
//                    localCBondSetB,
//                    newNeighborNumB,
//                    localBondNumB,
//                    newNeighborNumA,
//                    localIBondNeighborsA,
//                    localCBondNeighborsA);
//
//            targetP.process(setNumB, unmapped_atoms_molA,
//                    localMappingSize,
//                    localIBondNeighborsB,
//                    localMAPPING,
//                    counter,
//                    cTab1Copy,
//                    cTab2Copy,
//                    SignArray);
//
//            localBondNumB = targetP.getBondNumB();
//            newNeighborNumB = targetP.getNeighborBondNumB();
//            localCBondNeighborsB = targetP.getCBondNeighborsB();
//            localIBondSetB = targetP.getIBondSetB();
//            localCBondSetB = targetP.getCBondSetB();
//
////             System.out.println("Mapped Atoms before Iterator2: " + mapped_atoms);
//            iterator(no_further_MAPPINGS, localMappingSize, localMAPPING, newNeighborNumA, newNeighborNumB, localIBondNeighborsA, localIBondNeighborsB, localCBondNeighborsA, localCBondNeighborsB,
//                    localBondNumA, localBondNumB, localIBondSetA, localIBondSetB, localCBondSetA, localCBondSetB);
//            BESTARCS_copy.pop();
//        }
//        //System.out.println("Mapped Atoms before iterator Over: " + mapped_atoms);
//        return 0;
//    }
//
//    private List<Integer> findMcGregorMapping(List<Integer> MARCS, int mapped_atoms_num, List<Integer> current_MAPPING_org, int bondnum_A, List<Integer> i_bonds_A_org, int bondnum_B, List<Integer> i_bonds_B_org) {
//
//        List<Integer> currentMapping = new ArrayList<Integer>(current_MAPPING_org);
//        List<Integer> i_bonds_A = new ArrayList<Integer>(i_bonds_A_org);
//        List<Integer> i_bonds_B = new ArrayList<Integer>(i_bonds_B_org);
//        List<Integer> additional_mapping = new ArrayList<Integer>();
//
//
//
//        for (int x = 0; x < bondnum_A; x++) {
//            for (int y = 0; y < bondnum_B; y++) {
//
//
//                if (MARCS.get(x * bondnum_B + y) == 1) {
//
//
//                    int Atom1_moleculeA = i_bonds_A.get(x * 3 + 0);
//                    int Atom2_moleculeA = i_bonds_A.get(x * 3 + 1);
//                    int Atom1_moleculeB = i_bonds_B.get(y * 3 + 0);
//                    int Atom2_moleculeB = i_bonds_B.get(y * 3 + 1);
//
//                    IAtom R1_A = source.getAtom(Atom1_moleculeA);
//                    IAtom R2_A = source.getAtom(Atom2_moleculeA);
//                    IBond ReactantBond = source.getBond(R1_A, R2_A);
//
//                    IAtom P1_B = target.getAtom(Atom1_moleculeB);
//                    IAtom P2_B = target.getAtom(Atom2_moleculeB);
//                    IBond ProductBond = target.getBond(P1_B, P2_B);
//
////                  Bond Order Check Introduced by Asad
//
//                    boolean bMatch = McGregorChecks.bondMatch(ReactantBond, ProductBond);
//
//                    if ((bondTypeFlag && bMatch) || (!bondTypeFlag)) {
//
//                        for (int z = 0; z < mapped_atoms_num; z++) {
//
//                            int Mapped_Atom_1 = currentMapping.get(z * 2 + 0);
//                            int Mapped_Atom_2 = currentMapping.get(z * 2 + 1);
//
//                            if ((Mapped_Atom_1 == Atom1_moleculeA) && (Mapped_Atom_2 == Atom1_moleculeB)) {
//                                additional_mapping.add(Atom2_moleculeA);
//                                additional_mapping.add(Atom2_moleculeB);
//                            }
//
//                            if ((Mapped_Atom_1 == Atom1_moleculeA) && (Mapped_Atom_2 == Atom2_moleculeB)) {
//                                additional_mapping.add(Atom2_moleculeA);
//                                additional_mapping.add(Atom1_moleculeB);
//                            }
//
//                            if ((Mapped_Atom_1 == Atom2_moleculeA) && (Mapped_Atom_2 == Atom1_moleculeB)) {
//                                additional_mapping.add(Atom1_moleculeA);
//                                additional_mapping.add(Atom2_moleculeB);
//                            }
//
//                            if ((Mapped_Atom_1 == Atom2_moleculeA) && (Mapped_Atom_2 == Atom2_moleculeB)) {
//                                additional_mapping.add(Atom1_moleculeA);
//                                additional_mapping.add(Atom1_moleculeB);
//                            }
//                        }//for loop
//                    }
//                }
//            }
//        }
//
//
//        int additionalMappingSize = additional_mapping.size();
//
//        //add McGregorBondTypeInSensitive mapping to the Clique mapping
//        for (int a = 0; a < additionalMappingSize; a = a + 2) {
//            currentMapping.add(additional_mapping.get(a + 0));
//            currentMapping.add(additional_mapping.get(a + 1));
//        }
////        remove recurring mappings from currentMapping
//
//        List<Integer> unique_MAPPING = McGregorChecks.removeRecurringMappings(currentMapping);
//
//        return unique_MAPPING;
//    }
//
//    private List<Integer> partsearch(int xstart, int ystart, List<Integer> TempArcs_Org, int nNumGlobalA, int nNumGlobalB, List<Integer> iGlobalA, List<Integer> iGlobalB) {
//        int xIndex = xstart;
//        int yIndex = ystart;
//
//        List<Integer> TempArcs = new ArrayList<Integer>(TempArcs_Org);
//
//        if (TempArcs.get(xstart * nNumGlobalB + ystart) == 1) {
//
//            removeRedundantArcs(xstart, ystart, TempArcs, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
//            int arcsleft = 0;
//
//            for (int a = 0; a < nNumGlobalA; a++) {
//                for (int b = 0; b < nNumGlobalB; b++) {
//
//                    if (TempArcs.get(a * nNumGlobalB + b) == (1)) {
//                        arcsleft++;
//                    }
//                }
//            }
//
//            //test Bestarcsleft and skip rest if needed
//            if (arcsleft >= bestarcsleft) {
//                do {
//                    yIndex++;
//                    if (yIndex == nNumGlobalB) {
//                        yIndex = 0;
//                        xIndex++;
//
//                    }
//                } while ((xIndex < nNumGlobalA) && (TempArcs.get(xIndex * nNumGlobalB + yIndex) != 1)); //Correction by ASAD set value minus 1
//                if (xIndex < nNumGlobalA) {
//
//                    partsearch(xIndex, yIndex, TempArcs, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
//                    TempArcs.set(xIndex * nNumGlobalB + yIndex, 0);
//                    partsearch(xIndex, yIndex, TempArcs, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
//
//                } else {
//                    if (arcsleft > bestarcsleft) {
//                        McGregorChecks.removeTreeStructure(first);
//                        first = last = new BinaryTree(-1);
//                        last.equal = null;
//                        last.not_equal = null;
//
//                        while (!bestARCS.empty()) {
//                            bestARCS.pop();
//                        }
//                    }
//                    bestarcsleft = arcsleft;
//
//                    if (checkMARCS(TempArcs, nNumGlobalA, nNumGlobalB)) {
//                        bestARCS.push(TempArcs);
//                    }
//                }
//            }
//        } else {
//            do {
//                yIndex++;
//                if (yIndex == nNumGlobalB) {
//                    yIndex = 0;
//                    xIndex++;
//                }
//
//            } while ((xIndex < nNumGlobalA) && (TempArcs.get(xIndex * nNumGlobalB + yIndex) != 1)); //Correction by ASAD set value minus 1
//
//            if (xIndex < nNumGlobalA) {
//
//                partsearch(xIndex, yIndex, TempArcs, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
//                TempArcs.set(xIndex * nNumGlobalB + yIndex, 0);
//                partsearch(xIndex, yIndex, TempArcs, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
//            } else {
//                int arcsleft = 0;
//                for (int a = 0; a <
//                        nNumGlobalA; a++) {
//                    for (int b = 0; b <
//                            nNumGlobalB; b++) {
//                        if (TempArcs.get(a * nNumGlobalB + b) == 1) {
//                            arcsleft++;
//                        }
//
//                    }
//                }
//                if (arcsleft >= bestarcsleft) {
//                    if (arcsleft > bestarcsleft) {
//                        McGregorChecks.removeTreeStructure(first);
//                        first = last = new BinaryTree(-1);
//                        last.equal = null;
//                        last.not_equal = null;
//                        while (!bestARCS.empty()) {
//                            bestARCS.pop();
//                        }
//                    }
//                    bestarcsleft = arcsleft;
//
//                    if (checkMARCS(TempArcs, nNumGlobalA, nNumGlobalB)) {
//                        bestARCS.push(TempArcs);
//                    }
//
//                }
//            }
//        }
//        return TempArcs;
//    }
//
////The function is called in function partsearch. The function is given a temporary matrix and a position (row/column)
////within this matrix. First the function sets all entries to zero, which can be exlcuded in respect to the current
////atom by atom matching. After this the function replaces all entries in the same row and column of the current
////position by zeros. Only the entry of the current position is set to one.
////Return value "count_arcsleft" counts the number of arcs, which are still in the matrix.
//    private void removeRedundantArcs(int row, int column, List<Integer> MARCS, int nNumGlobalA, int nNumGlobalB,
//            List<Integer> iGlobalA, List<Integer> iGlobalB) {
//
//        int G1_atom = iGlobalA.get(row * 3 + 0);
//        int G2_atom = iGlobalA.get(row * 3 + 1);
//        int G3_atom = iGlobalB.get(column * 3 + 0);
//        int G4_atom = iGlobalB.get(column * 3 + 1);
//
//        for (int x = 0; x < nNumGlobalA; x++) {
//            int row_atom1 = iGlobalA.get(x * 3 + 0);
//            int row_atom2 = iGlobalA.get(x * 3 + 1);
//
//            for (int y = 0; y < nNumGlobalB; y++) {
//                int column_atom3 = iGlobalB.get(y * 3 + 0);
//                int column_atom4 = iGlobalB.get(y * 3 + 1);
//
//                if (McGregorChecks.cases(G1_atom, G2_atom, G3_atom, G4_atom, row_atom1, row_atom2, column_atom3, column_atom4)) {
//                    MARCS.set(x * nNumGlobalB + y, 0);
//                }
//
//            }
//        }
//
//        for (int v = 0; v < nNumGlobalA; v++) {
//            MARCS.set(v * nNumGlobalB + column, 0);
//        }
//
//        for (int w = 0; w < nNumGlobalB; w++) {
//            MARCS.set(row * nNumGlobalB + w, 0);
//        }
//
//        MARCS.set(row * nNumGlobalB + column, 1);
//    }
//
////The function is called in function partsearch. The function is given z temporary matrix.
////The function checks whether the temporary matrix is already found by calling the function
////"verifyNodes". If the matrix already exists the function returns false which means that
////the matrix will not be stored. Otherwise the function returns true which means that the
////matrix will be stored in function partsearch.
//    private boolean checkMARCS(List<Integer> MARCS_T, int nNumGlobalA, int nNumGlobalB) {
//
//        int size = nNumGlobalA * nNumGlobalA;
//        List<Integer> posnum_list = new ArrayList<Integer>(size);
//
//        for (int i = 0; i < posnum_list.size(); i++) {
//            posnum_list.add(i, 0);
//        }
//
//        int yCounter = 0;
//        int count_entries = 0;
//        for (int x = 0; x < (nNumGlobalA * nNumGlobalB); x++) {
//            if (MARCS_T.get(x) == 1) {
//                posnum_list.add(yCounter++, x);
//                count_entries++;
//            }
//        }
//        boolean flag = false;
//
//        verifyNodes(posnum_list, first, 0, count_entries);
//        if (newMatrix) {
//            flag = true;
//        }
//
//        return flag;
//
//    }
//
//    private boolean verifyNodes(List<Integer> matrix, BinaryTree currentStructure, int index, int fieldLength) {
//
//        if (((matrix.get(index) == currentStructure.getValue()) && (index < fieldLength)) && (currentStructure.equal != null)) {
//            newMatrix = false;
//            verifyNodes(matrix, currentStructure.equal, index + 1, fieldLength);
//        }
//        if (matrix.get(index) != currentStructure.getValue()) {
//            if (currentStructure.not_equal != null) {
//                verifyNodes(matrix, currentStructure.not_equal, index, fieldLength);
//            }
//
//            if (currentStructure.not_equal == null) {
//                currentStructure.not_equal = new BinaryTree(matrix.get(index));
//                currentStructure.not_equal.not_equal = null;
//                int yIndex = 0;
//
//
//                BinaryTree last_one = currentStructure.not_equal;
//
//                while ((yIndex + index + 1) < fieldLength) {
//                    last_one.equal = new BinaryTree(matrix.get(yIndex + index + 1));
//                    last_one = last_one.equal;
//                    last_one.not_equal = null;
//                    yIndex++;
//
//                }
//                last_one.equal = null;
//                newMatrix = true;
//            }
//
//        }
//        return true;
//    }
//
//    private void startsearch(int nNumGlobalA, int nNumGlobalB, List<Integer> iGlobalA, List<Integer> iGlobalB, List<Integer> modifiedARCS) {
//        int size = nNumGlobalA * nNumGlobalB;
//        List<Integer> FIXARCS = new ArrayList<Integer>(size);//  Initialize FIXARCS with 0
//        for (int i = 0; i < size; i++) {
//            FIXARCS.add(i, 0);
//        }
//
//        int xIndex = 0;
//        int yIndex = 0;
//
//        while ((xIndex < nNumGlobalA) && (modifiedARCS.get(xIndex * nNumGlobalB + yIndex) != 1)) {
//            yIndex++;
//            if (yIndex == nNumGlobalB) {
//                yIndex = 0;
//                xIndex++;
//            }
//        }
//
//        if (xIndex == nNumGlobalA) {
//            yIndex = nNumGlobalB - 1;
//            xIndex = xIndex - 1;
//        }
//
//        if (modifiedARCS.get(xIndex * nNumGlobalB + yIndex) == 0) {
//            partsearch(xIndex, yIndex, modifiedARCS, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
//        }
//
//        if (modifiedARCS.get(xIndex * nNumGlobalB + yIndex) != 0) {
//            partsearch(xIndex, yIndex, modifiedARCS, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
//            modifiedARCS.set(xIndex * nNumGlobalB + yIndex, 0);
//            partsearch(xIndex, yIndex, modifiedARCS, nNumGlobalA, nNumGlobalB, iGlobalA, iGlobalB);
//        }
//
//    }
//
//    /**
//     * @return
//     */
//    public List<List<Integer>> getMappings() {
//
//        return this.mappings;
//    }
//
//    /**
//     * @return
//     */
//    public int getMCSSize() {
//
//        return this.globalMCSSize;
//    }
//}


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

/**
 * @cdk.module smsd
 */
public class McGregor {

    private IAtomContainer source = null;
    private IAtomContainer target = null;
    private BinaryTree last = null;
    private BinaryTree first = null;
    private Stack<List<Integer>> bestARCS = null;
    private List<Integer> modifiedARCS = null;
    private List<Integer> i_globalA = null;
    private List<Integer> i_globalB = null;
    private List<String> c_globalA = null;
    private List<String> c_globalB = null;
    private List<String> c_tab1_copy = null;
    private List<String> c_tab2_copy = null;
    private List<Integer> i_bond_neighborsA = null;
    private List<String> c_bond_neighborsA = null;
    private int nNumGlobalA = 0;
    private int nNumGlobalB = 0;
    private int bestarcsleft = 0;
    private int globalMCSSize = 0;
    private List<List<Integer>> mappings = null;
    private int neighborBondnumA = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
    private int setBondNumA = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
    private int neighborBondNumB = 0; //number of remaining molecule B bonds after the clique search, which are neighbors of the MCS_1
    private int setBondNumB = 0; //number of remaining molecule B bonds after the clique search, which aren't neighbors
    /*This should be more or equal to all the atom types*/
    private static final String[] SignArray = {"$1", "$2", "$3", "$4", "$5", "$6", "$7", "$8", "$9", "$10", "$11", "$12",
        "$13", "$15", "$16", "$17", "$18", "$19", "$20", "$21", "$22", "$23", "$24",
        "$25", "$26", "$27", "$28", "$29", "$30", "$31", "$32", "$33", "$34", "$35", "$36",
        "$37", "$38", "$39", "$40", "$41", "$42", "$43", "$44", "$45", "$46",
        "$47", "$48", "$49", "$50", "$51", "$52", "$53", "$54", "$55"
    };
    protected boolean newMatrix = false;
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
        this.nNumGlobalA = 0;
        this.nNumGlobalB = 0;
        bestarcsleft = 0;

        if (_mappings.isEmpty()) {
            this.globalMCSSize = 0;
        } else {
            this.globalMCSSize = _mappings.get(0).size();
        }
        modifiedARCS = new ArrayList<Integer>();

        bestARCS = new Stack<List<Integer>>();

        //Initialization of global vectors
        i_globalA = new ArrayList<Integer>();
        i_globalB = new ArrayList<Integer>();
        c_globalA = new ArrayList<String>();
        c_globalB = new ArrayList<String>();

        c_tab1_copy = new ArrayList<String>();
        c_tab2_copy = new ArrayList<String>();

        newMatrix = false;
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
        List<Integer> mapped_atoms = new ArrayList<Integer>();
//        System.out.println("\nMapped Atoms");
        for (Map.Entry<Integer, Integer> map : present_Mapping.entrySet()) {
//            System.out.println("i:" + map.getKey() + " j:" + map.getValue());
            mapped_atoms.add(map.getKey());
            mapped_atoms.add(map.getValue());
        }
        int mapping_size = present_Mapping.size();


        i_bond_neighborsA = new ArrayList<Integer>();
        List<Integer> i_bond_setA = new ArrayList<Integer>();
        c_bond_neighborsA = new ArrayList<String>();
        List<String> c_bond_setA = new ArrayList<String>();

        List<Integer> i_bond_neighborsB = new ArrayList<Integer>();
        List<Integer> i_bond_setB = new ArrayList<Integer>();
        List<String> c_bond_neighborsB = new ArrayList<String>();
        List<String> c_bond_setB = new ArrayList<String>();

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

        int counter = 0;

        QueryProcessor queryProcess = new QueryProcessor(source, target, c_tab1_copy, c_tab2_copy, SignArray, neighborBondnumA, setBondNumA);


        queryProcess.process(unmapped_atoms_molA, mapping_size,
                i_bond_neighborsA, i_bond_setA,
                c_bond_neighborsA, c_bond_setA,
                mapped_atoms, counter);

        this.c_tab1_copy = queryProcess.getCTab1();
        this.c_tab2_copy = queryProcess.getCTab2();
        this.setBondNumA = queryProcess.getBondNumA();
        this.neighborBondnumA = queryProcess.getNeighborBondNumA();
        this.i_bond_neighborsA = queryProcess.getIBondNeighboursA();
        this.c_bond_neighborsA = queryProcess.getCBondNeighborsA();

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
        //new special signs -> the result are two vectors: c_bond_neighborsA and int_bonds_molB, which contain those
        //bonds of molecule B, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule A


        TargetProcessor targetProcess = new TargetProcessor(target, c_tab1_copy, c_tab2_copy, SignArray, neighborBondNumB, setBondNumB, neighborBondnumA, i_bond_neighborsA, c_bond_neighborsA);


        targetProcess.process(unmapped_atoms_molB,
                mapping_size,
                i_bond_neighborsB,
                i_bond_setB,
                c_bond_neighborsB,
                c_bond_setB,
                mapped_atoms,
                counter);

        this.c_tab1_copy = targetProcess.getCTab1();
        this.c_tab2_copy = targetProcess.getCTab2();
        this.setBondNumB = targetProcess.getBondNumB();
        this.neighborBondNumB = targetProcess.getNeighborBondNumB();

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
        this.globalMCSSize = (best_Mapping_size / 2);
        int SR_count = 0;


        c_tab1_copy.clear();
        generateCTab1Copy();


        c_tab2_copy.clear();
        generateCTab2Copy();

        //find mapped atoms of both molecules and store these in mapped_atoms
        List<Integer> mapped_atoms = new ArrayList<Integer>();

        int mapped_atom_number = 0;
//

        i_bond_neighborsA = new ArrayList<Integer>();
        List<Integer> i_bond_setA = new ArrayList<Integer>();
        c_bond_neighborsA = new ArrayList<String>();
        List<String> c_bond_setA = new ArrayList<String>();

        List<Integer> i_bond_neighborsB = new ArrayList<Integer>();
        List<Integer> i_bond_setB = new ArrayList<Integer>();
        List<String> c_bond_neighborsB = new ArrayList<String>();
        List<String> c_bond_setB = new ArrayList<String>();



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
        //new special signs -> the result are two vectors: c_bond_neighborsA and int_bonds_molA, which contain those
        //bonds of molecule A, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule B

        QueryProcessor queryProcess = new QueryProcessor(source, target, c_tab1_copy, c_tab2_copy, SignArray, neighborBondnumA, setBondNumA);


        queryProcess.process(unmapped_atoms_molA, clique_siz,
                i_bond_neighborsA, i_bond_setA,
                c_bond_neighborsA, c_bond_setA,
                mapped_atoms, SR_count);

        this.c_tab1_copy = queryProcess.getCTab1();
        this.c_tab2_copy = queryProcess.getCTab2();
        this.setBondNumA = queryProcess.getBondNumA();
        this.neighborBondnumA = queryProcess.getNeighborBondNumA();
        this.i_bond_neighborsA = queryProcess.getIBondNeighboursA();
        this.c_bond_neighborsA = queryProcess.getCBondNeighborsA();

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
        //new special signs -> the result are two vectors: c_bond_neighborsA and int_bonds_molB, which contain those
        //bonds of molecule B, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule A



        TargetProcessor targetProcess = new TargetProcessor(target, c_tab1_copy, c_tab2_copy, SignArray, neighborBondNumB, setBondNumB, neighborBondnumA, i_bond_neighborsA, c_bond_neighborsA);


        targetProcess.process(unmapped_atoms_molB,
                clique_siz,
                i_bond_neighborsB,
                i_bond_setB,
                c_bond_neighborsB,
                c_bond_setB,
                mapped_atoms,
                SR_count);

        this.c_tab1_copy = targetProcess.getCTab1();
        this.c_tab2_copy = targetProcess.getCTab2();
        this.setBondNumB = targetProcess.getBondNumB();
        this.neighborBondNumB = targetProcess.getNeighborBondNumB();

        boolean dummy = false;

        iterator(dummy, mapped_atom_number, mapped_atoms, neighborBondnumA, neighborBondNumB, i_bond_neighborsA, i_bond_neighborsB, c_bond_neighborsA, c_bond_neighborsB, setBondNumA, setBondNumB, i_bond_setA, i_bond_setB, c_bond_setA, c_bond_setB);

    }

    private int iterator(boolean mappingCheckFlag,
            int mapped_atom_num,
            List<Integer> mapped_atoms_org,
            int neighbor_bondnum_A,
            int neighbor_bondnum_B,
            List<Integer> i_bond_neighbor_atoms_A,
            List<Integer> i_bond_neighbor_atoms_B,
            List<String> c_bond_neighborsA,
            List<String> c_bond_neighborsB,
            int set_num_A, int set_num_B,
            List<Integer> i_bond_setA,
            List<Integer> i_bond_setB,
            List<String> c_bond_setA,
            List<String> c_bond_setB) throws IOException {

        List<Integer> mapped_atoms = new ArrayList<Integer>(mapped_atoms_org);

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
                IAtom R1_A = source.getAtom(Index_I);
                IAtom R2_A = source.getAtom(Index_IPlus1);
                IBond ReactantBond = source.getBond(R1_A, R2_A);

                //Get the atom index from the i_bond neighbor vactor
                int Index_J = i_bond_neighbor_atoms_B.get(column * 3 + 0);
                int Index_JPlus1 = i_bond_neighbor_atoms_B.get(column * 3 + 1);

                //Get the atoms
                IAtom P1_B = target.getAtom(Index_J);
                IAtom P2_B = target.getAtom(Index_JPlus1);
                IBond ProductBond = target.getBond(P1_B, P2_B);


                if (bondTypeFlag && bondMatch(ReactantBond, ProductBond)) {
                    if ((G1A.compareToIgnoreCase(G1B) == 0 && G2A.compareToIgnoreCase(G2B) == 0) || (G1A.compareToIgnoreCase(G2B) == 0 && G2A.compareToIgnoreCase(G1B) == 0)) {

                        no_further_mapping_possible = false;
                    }
//
                } else if ((!bondTypeFlag) &&
                        (((G1A.compareToIgnoreCase(G1B) == 0) && (G2A.compareToIgnoreCase(G2B) == 0)) || ((G1A.compareToIgnoreCase(G2B) == 0) && (G2A.compareToIgnoreCase(G1B) == 0)))) {
                    no_further_mapping_possible = false;
                }
            }
        }
        if (neighbor_bondnum_A == 0 || neighbor_bondnum_B == 0 || mappingCheckFlag || no_further_mapping_possible) {
            try {
                if (mapped_atom_num >= globalMCSSize) {
//                    System.out.println("Hello-1");
                    if (mapped_atom_num > globalMCSSize) {
//                        System.out.println("Hello-2");
                        this.globalMCSSize = mapped_atom_num;
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

        i_globalA.clear();
        i_globalB.clear();
        c_globalA.clear();
        c_globalB.clear();
        //redefining of global vectors and variables
        nNumGlobalA = neighbor_bondnum_A; //N global variable defined
        nNumGlobalB = neighbor_bondnum_B; //N global variable defined

        i_globalA.addAll(i_bond_neighbor_atoms_A);
        i_globalB.addAll(i_bond_neighbor_atoms_B);
        c_globalA.addAll(c_bond_neighborsA);
        c_globalB.addAll(c_bond_neighborsB);

        modifiedARCS.clear();
        int size = neighbor_bondnum_A * neighbor_bondnum_B;
//        modifiedARCS.setSize(size);
        for (int i = 0; i < size; i++) {
            modifiedARCS.add(i, 0);
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

                    IAtom R1_A = source.getAtom(Index_I);
                    IAtom R2_A = source.getAtom(Index_IPlus1);
                    IBond ReactantBond = source.getBond(R1_A, R2_A);

                    int Index_J = i_bond_neighbor_atoms_B.get(column * 3 + 0);
                    int Index_JPlus1 = i_bond_neighbor_atoms_B.get(column * 3 + 1);

                    IAtom P1_B = target.getAtom(Index_J);
                    IAtom P2_B = target.getAtom(Index_JPlus1);
                    IBond ProductBond = target.getBond(P1_B, P2_B);
                    if (bondTypeFlag && bondMatch(ReactantBond, ProductBond)) {
                        modifiedARCS.set(row * neighbor_bondnum_B + column, 1);
                    } else if (!bondTypeFlag) {
                        modifiedARCS.set(row * neighbor_bondnum_B + column, 1);
                    }
                }


            }
        }
        first = last = new BinaryTree(-1);
        last.equal = null;
        last.not_equal = null;
        bestarcsleft = 0;

        startsearch();
        Stack<List<Integer>> BESTARCS_copy = new Stack<List<Integer>>();

        BESTARCS_copy.addAll(bestARCS);
        while (!bestARCS.empty()) {
            bestARCS.pop();
        }
        while (!BESTARCS_copy.empty()) {

            List<Integer> MARCS_vector = new ArrayList<Integer>(BESTARCS_copy.peek());
            List<Integer> new_MAPPING = findMcGregorMapping(MARCS_vector, mapped_atom_num, mapped_atoms, neighbor_bondnum_A, i_bond_neighbor_atoms_A, neighbor_bondnum_B, i_bond_neighbor_atoms_B);

            int new_MAPPING_size = new_MAPPING.size() / 2;
            boolean no_further_MAPPINGS = false;
            if (mapped_atom_num == new_MAPPING_size) {
                no_further_MAPPINGS = true;
            }


            int new_neighbor_numA = 0; //instead of neighborBondnumA
            int new_neighbor_numB = 0; //instead of neighborBondNumB

            List<Integer> new_i_neighborsA = new ArrayList<Integer>(); //instead of i_bond_neighbor_atoms_A
            List<Integer> new_i_neighborsB = new ArrayList<Integer>(); //instead of i_bond_neighbor_atoms_B
            List<String> new_c_neighborsA = new ArrayList<String>(); //instead of c_bond_neighborsA
            List<String> new_c_neighborsB = new ArrayList<String>(); //instead of c_bond_neighborsB

            new_i_neighborsA.clear();
            new_i_neighborsB.clear();
            new_c_neighborsA.clear();
            new_c_neighborsB.clear();


            //new values for set_num_A + set_num_B
            //new arrays for i_bond_setA + i_bond_setB + c_bond_setB + c_bond_setB

            setBondNumA = 0; //instead of set_num_A
            setBondNumB = 0; //instead of set_num_B

            List<Integer> new_i_bond_setA = new ArrayList<Integer>(); //instead of i_bond_setA
            List<Integer> new_i_bond_setB = new ArrayList<Integer>(); //instead of i_bond_setB
            List<String> new_c_bond_setA = new ArrayList<String>(); //instead of c_bond_setA
            List<String> new_c_bond_setB = new ArrayList<String>(); //instead of c_bond_setB
            List<String> c_setB_copy = new ArrayList<String>();
            List<String> c_setA_copy = new ArrayList<String>();

            generateCSetACopy(set_num_A, c_bond_setA, c_setA_copy);
            generateCSetBCopy(set_num_B, c_bond_setB, c_setB_copy);

            //find unmapped atoms of molecule A
            List<Integer> unmapped_atoms_molA = new ArrayList<Integer>();
            int unmapped_numA = 0;
            boolean atomA_is_unmapped = true;

            for (int a = 0; a < source.getAtomCount(); a++) {
                for (int b = 0; b < new_MAPPING_size; b++) {
                    if (a == new_MAPPING.get(b * 2 + 0)) {
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
            boolean bond_considered = false;
            boolean normal_bond = true;
            for (int a = 0; a < set_num_A; a++) {

                int _elementAt_a = i_bond_setA.get(a * 3 + 0).intValue();
                for (int b = 0; b < unmapped_numA; b++) {
                    Integer unMappedAtomIndex = unmapped_atoms_molA.get(b);
                    if (unMappedAtomIndex == _elementAt_a) {
                        for (int c = 0; c < new_MAPPING_size; c++) {

                            if (new_MAPPING.get(c * 2 + 0).equals(i_bond_setA.get(a * 3 + 1))) {

                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 0));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 1));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 2));
                                new_c_neighborsA.add(c_setA_copy.get(a * 4 + 0));
                                if (c_setA_copy.get(a * 4 + 3).compareToIgnoreCase("X") == 0) {

                                    new_c_neighborsA.add(SignArray[counter]);
                                    new_c_neighborsA.add("X");
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    changeCharBonds(i_bond_setA.get(a * 3 + 1), SignArray[counter], set_num_A, i_bond_setA, c_setA_copy);
                                    int cor_atom = searchCorrespondingAtom(new_MAPPING_size, i_bond_setA.get(a * 3 + 1), 1, new_MAPPING);
                                    changeCharBonds(cor_atom, SignArray[counter], set_num_B, i_bond_setB, c_setB_copy);
                                    counter++;

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
                            setBondNumA++;

                        }
                        normal_bond = true;
                        bond_considered = true;
                    }
                    if (unMappedAtomIndex.equals(i_bond_setA.get(a * 3 + 1))) {
                        for (int c = 0; c < new_MAPPING_size; c++) {

                            if (new_MAPPING.get(c * 2 + 0).equals(i_bond_setA.get(a * 3 + 0))) {

                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 0));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 1));
                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 2));
                                if (c_setA_copy.get(a * 4 + 2).compareToIgnoreCase("X") == 0) {

                                    new_c_neighborsA.add(SignArray[counter]);
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 1));
                                    new_c_neighborsA.add(c_setA_copy.get(a * 4 + 0));
                                    new_c_neighborsA.add("X");
                                    changeCharBonds(i_bond_setA.get(a * 3 + 0), SignArray[counter], set_num_A, i_bond_setA, c_setA_copy);
                                    int cor_atom = searchCorrespondingAtom(new_MAPPING_size, i_bond_setA.get(a * 3 + 0), 1, new_MAPPING);
                                    changeCharBonds(cor_atom, SignArray[counter], set_num_B, i_bond_setB, c_setB_copy);
                                    counter++;

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
                            setBondNumA++;
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

            List<Integer> unmapped_atoms_molB = new ArrayList<Integer>();
            int unmapped_numB = 0;
            boolean atomB_is_unmapped = true;

            for (int a = 0; a < target.getAtomCount(); a++) {
                for (int b = 0; b < new_MAPPING_size; b++) {
                    if (a == new_MAPPING.get(b * 2 + 1)) {
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
                    if (unmapped_atoms_molB.get(b).equals(i_bond_setB.get(a * 3 + 0))) {
                        for (int c = 0; c < new_MAPPING_size; c++) {
                            if (new_MAPPING.get(c * 2 + 1).equals(i_bond_setB.get(a * 3 + 1))) {
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 0));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 1));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 2));
                                if (c_setB_copy.get(a * 4 + 3).compareToIgnoreCase("X") == 0) {
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add(SignArray[counter]);
                                    new_c_neighborsB.add("X");
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    changeCharBonds(i_bond_setB.get(a * 3 + 1), SignArray[counter], set_num_B, i_bond_setB, c_setB_copy);
                                    int cor_atom = searchCorrespondingAtom(new_MAPPING_size, i_bond_setB.get(a * 3 + 1), 2, new_MAPPING);
                                    changeCharBonds(cor_atom, SignArray[counter], new_neighbor_numA, new_i_neighborsA, new_c_neighborsA);
                                    counter++;

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
                            setBondNumB++;
                        }


                        normal_bond = true;
                        bond_considered = true;
                    }
                    if (unmapped_atoms_molB.get(b).equals(i_bond_setB.get(a * 3 + 1))) {
                        for (int c = 0; c < new_MAPPING_size; c++) {

                            if (new_MAPPING.get(c * 2 + 1).equals(i_bond_setB.get(a * 3 + 0))) {

                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 0));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 1));
                                new_i_neighborsB.add(i_bond_setB.get(a * 3 + 2));

                                if (c_setB_copy.get(a * 4 + 2).compareToIgnoreCase("X") == 0) {
                                    new_c_neighborsB.add(SignArray[counter]);
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 1));
                                    new_c_neighborsB.add(c_setB_copy.get(a * 4 + 0));
                                    new_c_neighborsB.add("X");
                                    changeCharBonds(i_bond_setB.get(a * 3 + 0), SignArray[counter], set_num_B, i_bond_setB, c_setB_copy);
                                    int cor_atom = searchCorrespondingAtom(new_MAPPING_size, i_bond_setB.get(a * 3 + 0), 2, new_MAPPING);
                                    changeCharBonds(cor_atom, SignArray[counter], new_neighbor_numA, new_i_neighborsA, new_c_neighborsA);
                                    counter++;
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
                            setBondNumB++;
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
                    setBondNumA, setBondNumB, new_i_bond_setA, new_i_bond_setB, new_c_bond_setA, new_c_bond_setB);
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
        IAtomContainer reactant = source;
        for (int a = 0; a < reactant.getBondCount(); a++) {
            String AtomI = reactant.getBond(a).getAtom(0).getSymbol();
            String AtomJ = reactant.getBond(a).getAtom(1).getSymbol();
            c_tab1_copy.add(AtomI);
            c_tab1_copy.add(AtomJ);
            c_tab1_copy.add("X");
            c_tab1_copy.add("X");
        }

        return 0;
    }

    private int generateCTab2Copy() throws IOException {
        IAtomContainer product = target;
        for (int a = 0; a < product.getBondCount(); a++) {
            String AtomI = product.getBond(a).getAtom(0).getSymbol();
            String AtomJ = product.getBond(a).getAtom(1).getSymbol();
            c_tab2_copy.add(AtomI);
            c_tab2_copy.add(AtomJ);
            c_tab2_copy.add("X");
            c_tab2_copy.add("X");
        }

        return 0;

    }

    private int changeCharBonds(int corresponding_atom, String new_symbol, int neighbor_bondnum, List<Integer> i_bond_neighbors, List<String> c_bond_neighbors) {

        for (int a = 0; a < neighbor_bondnum; a++) {
            if ((i_bond_neighbors.get(a * 3 + 0) == (corresponding_atom)) && (c_bond_neighbors.get(a * 4 + 2).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 2, c_bond_neighbors.get(a * 4 + 0));
                c_bond_neighbors.set(a * 4 + 0, new_symbol);
            }

            if ((i_bond_neighbors.get(a * 3 + 1) == (corresponding_atom)) && (c_bond_neighbors.get(a * 4 + 3).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 3, c_bond_neighbors.get(a * 4 + 1));
                c_bond_neighbors.set(a * 4 + 1, new_symbol);
            }

        }

        return 0;
    }

    /**
     *
     * @param ReactantBond
     * @param ProductBond
     * @return
     */
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

                    boolean bMatch = bondMatch(ReactantBond, ProductBond);

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

        List<Integer> unique_MAPPING = removeRecurringMappings(currentMapping);

        return unique_MAPPING;
    }

//Function compaires a structure array with itself. Sometimes a mapping occurs several times within the array.
//The function eliminates these recurring mappings. Function is called in function best_solution.
//The function is called by itself as long as the last list element is processed.
    private List<Integer> removeRecurringMappings(List<Integer> atom_mapping) {


        boolean exist = true;
        List<Integer> temp_map = new ArrayList<Integer>();
        int temp_counter = 0;
        int atom_mapping_size = atom_mapping.size();
        for (int x = 0; x < atom_mapping_size; x = x + 2) {
            int atom = atom_mapping.get(x);
            for (int y = x + 2; y < atom_mapping_size; y = y + 2) {
                if (atom == atom_mapping.get(y)) {
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

    private void partsearch(int xstart, int ystart, List<Integer> TEMPMARCS_ORG) {
        int xIndex = xstart;
        int yIndex = ystart;

        List<Integer> TEMPMARCS = new ArrayList<Integer>(TEMPMARCS_ORG);

        if (TEMPMARCS.get(xstart * nNumGlobalB + ystart) == 1) {

            removeRedundantArcs(xstart, ystart, TEMPMARCS);
            int arcsleft = 0;

            for (int a = 0; a < nNumGlobalA; a++) {
                for (int b = 0; b < nNumGlobalB; b++) {

                    if (TEMPMARCS.get(a * nNumGlobalB + b) == (1)) {
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
                } while ((xIndex < nNumGlobalA) && (TEMPMARCS.get(xIndex * nNumGlobalB + yIndex) != 1)); //Correction by ASAD set value minus 1
                if (xIndex < nNumGlobalA) {

                    partsearch(xIndex, yIndex, TEMPMARCS);
                    TEMPMARCS.set(xIndex * nNumGlobalB + yIndex, 0);
                    partsearch(xIndex, yIndex, TEMPMARCS);

                } else {
                    if (arcsleft > bestarcsleft) {
                        removeTreeStructure(first);
                        first = last = new BinaryTree(-1);
                        last.equal = null;
                        last.not_equal = null;

                        while (!bestARCS.empty()) {
                            bestARCS.pop();
                        }
                    }
                    bestarcsleft = arcsleft;

                    if (checkMARCS(TEMPMARCS)) {
                        bestARCS.push(TEMPMARCS);
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

            } while ((xIndex < nNumGlobalA) && (TEMPMARCS.get(xIndex * nNumGlobalB + yIndex) != 1)); //Correction by ASAD set value minus 1

            if (xIndex < nNumGlobalA) {

                partsearch(xIndex, yIndex, TEMPMARCS);
                TEMPMARCS.set(xIndex * nNumGlobalB + yIndex, 0);
                partsearch(xIndex, yIndex, TEMPMARCS);
            } else {
                int arcsleft = 0;
                for (int a = 0; a <
                        nNumGlobalA; a++) {
                    for (int b = 0; b <
                            nNumGlobalB; b++) {
                        if (TEMPMARCS.get(a * nNumGlobalB + b) == 1) {
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
                        while (!bestARCS.empty()) {
                            bestARCS.pop();
                        }
                    }
                    bestarcsleft = arcsleft;

                    if (checkMARCS(TEMPMARCS)) {
                        bestARCS.push(TEMPMARCS);
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
    private void removeRedundantArcs(int row, int column, List<Integer> MARCS) {

        int G1_atom = i_globalA.get(row * 3 + 0);
        int G2_atom = i_globalA.get(row * 3 + 1);
        int G3_atom = i_globalB.get(column * 3 + 0);
        int G4_atom = i_globalB.get(column * 3 + 1);

        for (int x = 0; x < nNumGlobalA; x++) {
            int row_atom1 = i_globalA.get(x * 3 + 0);
            int row_atom2 = i_globalA.get(x * 3 + 1);

            for (int y = 0; y < nNumGlobalB; y++) {
                int column_atom3 = i_globalB.get(y * 3 + 0);
                int column_atom4 = i_globalB.get(y * 3 + 1);

                if (cases(G1_atom, G2_atom, G3_atom, G4_atom, row_atom1, row_atom2, column_atom3, column_atom4)) {
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

    private boolean case1(int G1_atom, int G3_atom, int G4_atom, int row_atom1, int row_atom2, int column_atom3, int column_atom4) {
        if (((G1_atom == row_atom1) || (G1_atom == row_atom2)) &&
                (!(((column_atom3 == G3_atom) || (column_atom4 == G3_atom)) || ((column_atom3 == G4_atom) || (column_atom4 == G4_atom))))) {
            return true;
        }
        return false;
    }

    private boolean case2(int G2_atom, int G3_atom, int G4_atom, int row_atom1, int row_atom2, int column_atom3, int column_atom4) {
        if (((G2_atom == row_atom1) ||
                (G2_atom == row_atom2)) &&
                (!(((column_atom3 == G3_atom) || (column_atom4 == G3_atom)) || ((column_atom3 == G4_atom) || (column_atom4 == G4_atom))))) {
            return true;
        }
        return false;
    }

    private boolean case3(int G1_atom, int G3_atom, int G2_atom, int row_atom1, int row_atom2, int column_atom3, int column_atom4) {
        if (((G3_atom == column_atom3) || (G3_atom == column_atom4)) &&
                (!(((row_atom1 == G1_atom) || (row_atom2 == G1_atom)) || ((row_atom1 == G2_atom) || (row_atom2 == G2_atom))))) {
            return true;
        }
        return false;
    }

    private boolean case4(int G1_atom, int G2_atom, int G4_atom, int row_atom1, int row_atom2, int column_atom3, int column_atom4) {
        if (((G4_atom == column_atom3) || (G4_atom == column_atom4)) &&
                (!(((row_atom1 == G1_atom) || (row_atom2 == G1_atom)) || ((row_atom1 == G2_atom) || (row_atom2 == G2_atom))))) {
            return true;
        }
        return false;
    }

    private boolean cases(int G1_atom, int G2_atom, int G3_atom, int G4_atom, int row_atom1, int row_atom2, int column_atom3, int column_atom4) {
        if (case1(G1_atom, G3_atom, G4_atom, row_atom1, row_atom2, column_atom3, column_atom4) || case2(G2_atom, G3_atom, G4_atom, row_atom1, row_atom2, column_atom3, column_atom4) || case3(G1_atom, G3_atom, G2_atom, row_atom1, row_atom2, column_atom3, column_atom4) || case4(G1_atom, G2_atom, G4_atom, row_atom1, row_atom2, column_atom3, column_atom4)) {
            return true;
        }
        return false;
    }

    /*
     * Modified function call by ASAD in Java have to check
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
    private boolean checkMARCS(List<Integer> MARCS_T) {

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

    private int searchCorrespondingAtom(int mappedAtomsSize, int atomFromOtherMolecule, int molecule, List<Integer> mapped_atoms_org) {


        List<Integer> mapped_atoms = new ArrayList<Integer>(mapped_atoms_org);

        int corresponding_atom = 0;
        for (int a = 0; a < mappedAtomsSize; a++) {
            if ((molecule == 1) && (mapped_atoms.get(a * 2 + 0).intValue() == atomFromOtherMolecule)) {
                corresponding_atom = mapped_atoms.get(a * 2 + 1);
            }
            if ((molecule == 2) && (mapped_atoms.get(a * 2 + 1).intValue() == atomFromOtherMolecule)) {
                corresponding_atom = mapped_atoms.get(a * 2 + 0);
            }
        }
        return corresponding_atom;
    }

    private int generateCSetBCopy(int bond_number, List<String> c_setB, List<String> c_setB_copy) {

        for (int a = 0; a < bond_number; a++) {
            c_setB_copy.add(c_setB.get(a * 4 + 0));
            c_setB_copy.add(c_setB.get(a * 4 + 1));
            c_setB_copy.add("X");
            c_setB_copy.add("X");
        }

        return 0;
    }

    private int generateCSetACopy(int bond_number, List<String> c_setA, List<String> c_setA_copy) {

        for (int a = 0; a < bond_number; a++) {
            c_setA_copy.add(c_setA.get(a * 4 + 0));
            c_setA_copy.add(c_setA.get(a * 4 + 1));
            c_setA_copy.add("X");
            c_setA_copy.add("X");
        }

        return 0;
    }

    private void startsearch() {
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
            partsearch(xIndex, yIndex, modifiedARCS);
        }

        if (modifiedARCS.get(xIndex * nNumGlobalB + yIndex) != 0) {
            partsearch(xIndex, yIndex, modifiedARCS);
            modifiedARCS.set(xIndex * nNumGlobalB + yIndex, 0);
            partsearch(xIndex, yIndex, modifiedARCS);
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