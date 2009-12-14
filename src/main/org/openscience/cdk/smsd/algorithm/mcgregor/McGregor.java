
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
    private List<Integer> modifiedARCS = null;
    private int bestarcsleft = 0;
    private int globalMCSSize = 0;
    private List<List<Integer>> mappings = null;
    /*This should be more or equal to all the atom types*/
    private static final String[] SignArray = {
        "$1", "$2", "$3", "$4", "$5", "$6", "$7", "$8", "$9", "$10", "$11", "$12",
        "$13", "$15", "$16", "$17", "$18", "$19", "$20", "$21", "$22", "$23", "$24",
        "$25", "$26", "$27", "$28", "$29", "$30", "$31", "$32", "$33", "$34", "$35",
        "$36", "$37", "$38", "$39", "$40", "$41", "$42", "$43", "$44", "$45", "$46",
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
        bestarcsleft = 0;

        if (_mappings.isEmpty()) {
            this.globalMCSSize = 0;
        } else {
            this.globalMCSSize = _mappings.get(0).size();
        }
        modifiedARCS = new ArrayList<Integer>();
        bestARCS = new Stack<List<Integer>>();
        newMatrix = false;
    }

    /**
     *
     * @param largestMappingSize
     * @param present_Mapping
     * @throws IOException
     */
    public void startMcGregorIteration(int largestMappingSize, Map<Integer, Integer> present_Mapping) throws IOException {

        this.globalMCSSize = (largestMappingSize / 2);
        List<String> c_tab1_copy = McGregorChecks.generateCTabCopy(source);

        List<String> c_tab2_copy = McGregorChecks.generateCTabCopy(target);


        //find mapped atoms of both molecules and store these in mappedAtoms
        List<Integer> mapped_atoms = new ArrayList<Integer>();
//        System.out.println("\nMapped Atoms");
        for (Map.Entry<Integer, Integer> map : present_Mapping.entrySet()) {
//            System.out.println("i:" + map.getKey() + " j:" + map.getValue());
            mapped_atoms.add(map.getKey());
            mapped_atoms.add(map.getValue());
        }
        int mapping_size = present_Mapping.size();


        List<Integer> i_bond_neighborsA = new ArrayList<Integer>();
        List<String> c_bond_neighborsA = new ArrayList<String>();

        List<Integer> i_bond_setA = new ArrayList<Integer>();
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
        int gSetBondNumA = 0;
        int gSetBondNumB = 0;
        int gNeighborBondnumA = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
        int gNeighborBondNumB = 0; //number of remaining molecule B bonds after the clique search, which are neighbors of the MCS_1



        QueryProcessor queryProcess = new QueryProcessor(
                c_tab1_copy,
                c_tab2_copy,
                SignArray,
                gNeighborBondnumA,
                gSetBondNumA,
                i_bond_neighborsA,
                c_bond_neighborsA);


        queryProcess.process(
                source,
                target,
                unmapped_atoms_molA,
                mapping_size,
                i_bond_setA,
                c_bond_setA,
                mapped_atoms,
                counter);

        c_tab1_copy = queryProcess.getCTab1();
        c_tab2_copy = queryProcess.getCTab2();
        gSetBondNumA = queryProcess.getBondNumA();
        gNeighborBondnumA = queryProcess.getNeighborBondNumA();
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
                c_tab1_copy,
                c_tab2_copy,
                SignArray,
                gNeighborBondNumB,
                gSetBondNumB,
                i_bond_neighborsB,
                c_bond_neighborsB,
                gNeighborBondnumA,
                i_bond_neighborsA,
                c_bond_neighborsA);


        targetProcess.process(
                target,
                unmapped_atoms_molB,
                mapping_size,
                i_bond_setB,
                c_bond_setB,
                mapped_atoms,
                counter);

        c_tab1_copy = targetProcess.getCTab1();
        c_tab2_copy = targetProcess.getCTab2();
        gSetBondNumB = targetProcess.getBondNumB();
        gNeighborBondNumB = targetProcess.getNeighborBondNumB();
        i_bond_neighborsB = targetProcess.getIBondNeighboursB();
        c_bond_neighborsB = targetProcess.getCBondNeighborsB();

        boolean dummy = false;

        iterator(dummy,
                present_Mapping.size(),
                mapped_atoms,
                gNeighborBondnumA,
                gNeighborBondNumB,
                i_bond_neighborsA,
                i_bond_neighborsB,
                c_bond_neighborsA,
                c_bond_neighborsB,
                gSetBondNumA,
                gSetBondNumB,
                i_bond_setA,
                i_bond_setB,
                c_bond_setA,
                c_bond_setB);

    }

    /**
     *
     * @param largestMappingSize
     * @param clique_vector
     * @param comp_graph_nodes
     * @throws IOException
     */
    public void startMcGregorIteration(int largestMappingSize, List<Integer> clique_vector, List<Integer> comp_graph_nodes) throws IOException {
        this.globalMCSSize = (largestMappingSize / 2);
        List<String> c_tab1_copy = McGregorChecks.generateCTabCopy(source);

        List<String> c_tab2_copy = McGregorChecks.generateCTabCopy(target);

        //find mapped atoms of both molecules and store these in mappedAtoms
        List<Integer> mapped_atoms = new ArrayList<Integer>();

        int mappedAtomCount = 0;


        List<Integer> iBondNeighborAtomsA = new ArrayList<Integer>();
        List<String> cBondNeighborsA = new ArrayList<String>();

        List<Integer> i_bond_setA = new ArrayList<Integer>();
        List<String> c_bond_setA = new ArrayList<String>();

        List<Integer> iBondNeighborAtomsB = new ArrayList<Integer>();
        List<Integer> i_bond_setB = new ArrayList<Integer>();
        List<String> cBondNeighborsB = new ArrayList<String>();
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
                    mappedAtomCount++;
                }
            }
        }


        //find unmapped atoms of molecule A
        List<Integer> unmapped_atoms_molA = new ArrayList<Integer>();

        int unmapped_numA = 0;
        boolean atomA_is_unmapped = true;

//        System.out.println("Mapped Atoms: " + mappedAtoms);

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

        int counter = 0;
        int setNumA = 0;
        int setNumB = 0;
        int localNeighborBondnumA = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
        int localNeighborBondNumB = 0; //number of remaining molecule B bonds after the clique search, which are neighbors of the MCS_1

        //Extract bonds which are related with unmapped atoms of molecule A.
        //In case that unmapped atoms are connected with already mapped atoms, the mapped atoms are labelled with
        //new special signs -> the result are two vectors: cBondNeighborsA and int_bonds_molA, which contain those
        //bonds of molecule A, which are relevant for the McGregorBondTypeInSensitive algorithm.
        //The special signs must be transfered to the corresponding atoms of molecule B

        QueryProcessor queryProcess = new QueryProcessor(
                c_tab1_copy,
                c_tab2_copy,
                SignArray,
                localNeighborBondnumA,
                setNumA,
                iBondNeighborAtomsA,
                cBondNeighborsA);


        queryProcess.process(
                source,
                target,
                unmapped_atoms_molA,
                clique_siz,
                i_bond_setA,
                c_bond_setA,
                mapped_atoms,
                counter);

        c_tab1_copy = queryProcess.getCTab1();
        c_tab2_copy = queryProcess.getCTab2();
        setNumA = queryProcess.getBondNumA();
        localNeighborBondnumA = queryProcess.getNeighborBondNumA();
        iBondNeighborAtomsA = queryProcess.getIBondNeighboursA();
        cBondNeighborsA = queryProcess.getCBondNeighborsA();

        //find unmapped atoms of molecule B
        List<Integer> unmapped_atoms_molB = new ArrayList<Integer>();
        int unmapped_numB = 0;
        boolean atomB_is_unmapped = true;

//        System.out.println("localNeighborBondnumA After:" + localNeighborBondnumA);

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
                c_tab1_copy,
                c_tab2_copy,
                SignArray,
                localNeighborBondNumB,
                setNumB,
                iBondNeighborAtomsB,
                cBondNeighborsB,
                localNeighborBondnumA,
                iBondNeighborAtomsA,
                cBondNeighborsA);


        targetProcess.process(
                target,
                unmapped_atoms_molB,
                clique_siz,
                i_bond_setB,
                c_bond_setB,
                mapped_atoms,
                counter);

        c_tab1_copy = targetProcess.getCTab1();
        c_tab2_copy = targetProcess.getCTab2();
        setNumB = targetProcess.getBondNumB();
        localNeighborBondNumB = targetProcess.getNeighborBondNumB();
        iBondNeighborAtomsB = targetProcess.getIBondNeighboursB();
        cBondNeighborsB = targetProcess.getCBondNeighborsB();

        boolean dummy = false;

        iterator(dummy,
                mappedAtomCount,
                mapped_atoms,
                localNeighborBondnumA,
                localNeighborBondNumB,
                iBondNeighborAtomsA,
                iBondNeighborAtomsB,
                cBondNeighborsA,
                cBondNeighborsB,
                setNumA,
                setNumB,
                i_bond_setA,
                i_bond_setB,
                c_bond_setA,
                c_bond_setB);

    }

    private int iterator(boolean mappingCheckFlag,
            int mappedAtomCount,
            List<Integer> mappedAtomsOrg,
            int neighborBondNumA,
            int neighborBondNumB,
            List<Integer> iBondNeighborAtomsA,
            List<Integer> iBondNeighborAtomsB,
            List<String> cBondNeighborsA,
            List<String> cBondNeighborsB,
            int setNumA, int setNumB,
            List<Integer> i_bond_setA,
            List<Integer> i_bond_setB,
            List<String> c_bond_setA,
            List<String> c_bond_setB) throws IOException {

        List<Integer> mappedAtoms = new ArrayList<Integer>(mappedAtomsOrg);

//        //check possible mappings:
        boolean furtherMappingFlag = McGregorChecks.isFurtherMappingPossible(source, target, neighborBondNumA, neighborBondNumB, iBondNeighborAtomsA, iBondNeighborAtomsB, cBondNeighborsA, cBondNeighborsB);

        if (neighborBondNumA == 0 || neighborBondNumB == 0 || mappingCheckFlag || !furtherMappingFlag) {
            setFinalMappings(mappedAtoms, mappedAtomCount);
            return 0;
        }

        modifiedARCS.clear();
        int size = neighborBondNumA * neighborBondNumB;
        for (int i = 0; i < size; i++) {
            modifiedARCS.add(i, 0);
        }
        setModifedArcs(neighborBondNumA,
                neighborBondNumB,
                iBondNeighborAtomsA,
                iBondNeighborAtomsB,
                cBondNeighborsA,
                cBondNeighborsB);
        first = last = new BinaryTree(-1);
        last.equal = null;
        last.not_equal = null;
        bestarcsleft = 0;

        startsearch(iBondNeighborAtomsA, iBondNeighborAtomsB, neighborBondNumA, neighborBondNumB);
        Stack<List<Integer>> BESTARCS_copy = new Stack<List<Integer>>();

        BESTARCS_copy.addAll(bestARCS);
        while (!bestARCS.empty()) {
            bestARCS.pop();
        }
        searchAndExtendMappings(BESTARCS_copy,
                mappedAtoms,
                mappedAtomCount,
                neighborBondNumA,
                neighborBondNumB,
                iBondNeighborAtomsA,
                iBondNeighborAtomsB,
                setNumA,
                setNumB,
                i_bond_setA,
                i_bond_setB,
                c_bond_setA,
                c_bond_setB);

        //System.out.println("In the iterator Termination");
        //System.out.println("============+++++++++==============");
        //System.out.println("Mapped Atoms before iterator Over: " + mappedAtoms);
        return 0;
    }

    private List<Integer> findMcGregorMapping(List<Integer> MARCS, int mapped_atoms_num, List<Integer> current_MAPPING_org, int bondnum_A, List<Integer> i_bonds_A, int bondnum_B, List<Integer> i_bonds_B) {

        List<Integer> currentMapping = new ArrayList<Integer>(current_MAPPING_org);
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

    private void setModifedArcs(int neighborBondNumA,
            int neighborBondNumB,
            List<Integer> iBondNeighborAtomsA,
            List<Integer> iBondNeighborAtomsB,
            List<String> cBondNeighborsA,
            List<String> cBondNeighborsB) {
        for (int row = 0; row < neighborBondNumA; row++) {
            for (int column = 0; column < neighborBondNumB; column++) {

                String G1A = cBondNeighborsA.get(row * 4 + 0);
                String G2A = cBondNeighborsA.get(row * 4 + 1);
                String G1B = cBondNeighborsB.get(column * 4 + 0);
                String G2B = cBondNeighborsB.get(column * 4 + 1);


                if ((G1A.compareToIgnoreCase(G1B) == 0 && G2A.compareToIgnoreCase(G2B) == 0) || (G1A.compareToIgnoreCase(G2B) == 0 && G2A.compareToIgnoreCase(G1B) == 0)) {


                    int Index_I = iBondNeighborAtomsA.get(row * 3 + 0);
                    int Index_IPlus1 = iBondNeighborAtomsA.get(row * 3 + 1);

                    IAtom R1_A = source.getAtom(Index_I);
                    IAtom R2_A = source.getAtom(Index_IPlus1);
                    IBond ReactantBond = source.getBond(R1_A, R2_A);

                    int Index_J = iBondNeighborAtomsB.get(column * 3 + 0);
                    int Index_JPlus1 = iBondNeighborAtomsB.get(column * 3 + 1);

                    IAtom P1_B = target.getAtom(Index_J);
                    IAtom P2_B = target.getAtom(Index_JPlus1);
                    IBond ProductBond = target.getBond(P1_B, P2_B);
                    if ((bondTypeFlag && McGregorChecks.bondMatch(ReactantBond, ProductBond)) || (!bondTypeFlag)) {
                        modifiedARCS.set(row * neighborBondNumB + column, 1);
                    }
                }
            }
        }
    }

    private void partsearch(int xstart, int ystart, List<Integer> TEMPMARCS_ORG, List<Integer> iBondNeighborAtomsA, List<Integer> iBondNeighborAtomsB, int neighborBondNumA, int neighborBondNumB) {
        int xIndex = xstart;
        int yIndex = ystart;

        List<Integer> TEMPMARCS = new ArrayList<Integer>(TEMPMARCS_ORG);

        if (TEMPMARCS.get(xstart * neighborBondNumB + ystart) == 1) {

            McGregorChecks.removeRedundantArcs(xstart, ystart, TEMPMARCS, iBondNeighborAtomsA, iBondNeighborAtomsB, neighborBondNumA, neighborBondNumB);
            int arcsleft = McGregorChecks.countArcsLeft(TEMPMARCS, neighborBondNumA, neighborBondNumB);

            //test Best arcs left and skip rest if needed
            if (arcsleft >= bestarcsleft) {
                setArcs(xIndex, yIndex, arcsleft, TEMPMARCS, iBondNeighborAtomsA, iBondNeighborAtomsB, neighborBondNumA, neighborBondNumB);
            }
        } else {
            do {
                yIndex++;
                if (yIndex == neighborBondNumB) {
                    yIndex = 0;
                    xIndex++;
                }

            } while ((xIndex < neighborBondNumA) && (TEMPMARCS.get(xIndex * neighborBondNumB + yIndex) != 1)); //Correction by ASAD set value minus 1

            if (xIndex < neighborBondNumA) {

                partsearch(xIndex, yIndex, TEMPMARCS, iBondNeighborAtomsA, iBondNeighborAtomsB, neighborBondNumA, neighborBondNumB);
                TEMPMARCS.set(xIndex * neighborBondNumB + yIndex, 0);
                partsearch(xIndex, yIndex, TEMPMARCS, iBondNeighborAtomsA, iBondNeighborAtomsB, neighborBondNumA, neighborBondNumB);
            } else {
                int arcsleft = McGregorChecks.countArcsLeft(TEMPMARCS, neighborBondNumA, neighborBondNumB);
                if (arcsleft >= bestarcsleft) {
                    popBestArcs(arcsleft);

                    if (checkMARCS(TEMPMARCS, neighborBondNumA, neighborBondNumB)) {
                        bestARCS.push(TEMPMARCS);
                    }

                }
            }
        }
    }

//The function is called in function partsearch. The function is given z temporary matrix.
//The function checks whether the temporary matrix is already found by calling the function
//"verifyNodes". If the matrix already exists the function returns false which means that
//the matrix will not be stored. Otherwise the function returns true which means that the
//matrix will be stored in function partsearch.
    private boolean checkMARCS(List<Integer> MARCS_T, int neighborBondNumA, int neighborBondNumB) {

        int size = neighborBondNumA * neighborBondNumA;
        List<Integer> posnum_list = new ArrayList<Integer>(size);

        for (int i = 0; i < posnum_list.size(); i++) {
            posnum_list.add(i, 0);
        }

        int yCounter = 0;
        int count_entries = 0;
        for (int x = 0; x < (neighborBondNumA * neighborBondNumB); x++) {
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
        if (index < fieldLength) {
            if (matrix.get(index) == currentStructure.getValue() && currentStructure.equal != null) {
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
        }
        return true;
    }

    private void startsearch(List<Integer> iBondNeighborAtomsA, List<Integer> iBondNeighborAtomsB, int neighborBondNumA, int neighborBondNumB) {
        int size = neighborBondNumA * neighborBondNumB;
        List<Integer> FIXARCS = new ArrayList<Integer>(size);//  Initialize FIXARCS with 0
        for (int i = 0; i < size; i++) {
            FIXARCS.add(i, 0);
        }

        int xIndex = 0;
        int yIndex = 0;

        while ((xIndex < neighborBondNumA) && (modifiedARCS.get(xIndex * neighborBondNumB + yIndex) != 1)) {
            yIndex++;
            if (yIndex == neighborBondNumB) {
                yIndex = 0;
                xIndex++;
            }
        }

        if (xIndex == neighborBondNumA) {
            yIndex = neighborBondNumB - 1;
            xIndex = xIndex - 1;
        }

        if (modifiedARCS.get(xIndex * neighborBondNumB + yIndex) == 0) {
            partsearch(xIndex, yIndex, modifiedARCS, iBondNeighborAtomsA, iBondNeighborAtomsB, neighborBondNumA, neighborBondNumB);
        }

        if (modifiedARCS.get(xIndex * neighborBondNumB + yIndex) != 0) {
            partsearch(xIndex, yIndex, modifiedARCS, iBondNeighborAtomsA, iBondNeighborAtomsB, neighborBondNumA, neighborBondNumB);
            modifiedARCS.set(xIndex * neighborBondNumB + yIndex, 0);
            partsearch(xIndex, yIndex, modifiedARCS, iBondNeighborAtomsA, iBondNeighborAtomsB, neighborBondNumA, neighborBondNumB);
        }

    }

    /**
     * @return mappings
     */
    public List<List<Integer>> getMappings() {

        return this.mappings;
    }

    /**
     * @return MCS size
     */
    public int getMCSSize() {

        return this.globalMCSSize;
    }

    private void setFinalMappings(List<Integer> mapped_atoms, int mappedAtomCount) {
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
    }

    private void searchAndExtendMappings(
            Stack<List<Integer>> BESTARCS_copy,
            List<Integer> mappedAtoms,
            int mappedAtomCount,
            int neighborBondNumA,
            int neighborBondNumB,
            List<Integer> iBondNeighborAtomsA,
            List<Integer> iBondNeighborAtomsB,
            int setNumA, int setNumB,
            List<Integer> i_bond_setA,
            List<Integer> i_bond_setB,
            List<String> c_bond_setA,
            List<String> c_bond_setB) throws IOException {

        while (!BESTARCS_copy.empty()) {

            List<Integer> MARCS_vector = new ArrayList<Integer>(BESTARCS_copy.peek());
            List<Integer> new_Mapping = findMcGregorMapping(MARCS_vector, mappedAtomCount, mappedAtoms, neighborBondNumA, iBondNeighborAtomsA, neighborBondNumB, iBondNeighborAtomsB);

            int newMapingSize = new_Mapping.size() / 2;
            boolean no_further_MAPPINGS = false;
            if (mappedAtomCount == newMapingSize) {
                no_further_MAPPINGS = true;
            }

            List<Integer> new_i_neighborsA = new ArrayList<Integer>(); //instead of iBondNeighborAtomsA
            List<Integer> new_i_neighborsB = new ArrayList<Integer>(); //instead of iBondNeighborAtomsB
            List<String> new_c_neighborsA = new ArrayList<String>(); //instead of cBondNeighborsA
            List<String> new_c_neighborsB = new ArrayList<String>(); //instead of cBondNeighborsB
            List<Integer> new_i_bond_setA = new ArrayList<Integer>(); //instead of i_bond_setA
            List<Integer> new_i_bond_setB = new ArrayList<Integer>(); //instead of i_bond_setB
            List<String> new_c_bond_setA = new ArrayList<String>(); //instead of c_bond_setA
            List<String> new_c_bond_setB = new ArrayList<String>(); //instead of c_bond_setB
            //new values for setNumA + setNumB
            //new arrays for i_bond_setA + i_bond_setB + c_bond_setB + c_bond_setB

            List<String> c_setA_copy = McGregorChecks.generateCSetCopy(setNumA, c_bond_setA);
            List<String> c_setB_copy = McGregorChecks.generateCSetCopy(setNumB, c_bond_setB);

            //find unmapped atoms of molecule A
            List<Integer> unmapped_atoms_molA = new ArrayList<Integer>();
            int unmapped_numA = 0;
            boolean atomA_is_unmapped = true;

            for (int a = 0; a < source.getAtomCount(); a++) {
                for (int b = 0; b < newMapingSize; b++) {
                    if (a == new_Mapping.get(b * 2 + 0)) {
                        atomA_is_unmapped = false;
                    }

                }
                if (atomA_is_unmapped) {
                    unmapped_atoms_molA.add(unmapped_numA++, a);
                }
                atomA_is_unmapped = true;
            }


            //The special signs must be transfered to the corresponding atoms of molecule B

            int counter = 0;
            //number of remaining molecule A bonds after the clique search, which aren't neighbors
            int newSetBondNumA = 0; //instead of setNumA
            int newNeighborNumA = 0; //instead of localNeighborBondnumA


            QueryProcessor queryProcess =
                    new QueryProcessor(
                    c_setA_copy,
                    c_setB_copy,
                    SignArray,
                    newNeighborNumA,
                    newSetBondNumA,
                    new_i_neighborsA,
                    new_c_neighborsA);


            queryProcess.process(
                    setNumA,
                    setNumB,
                    i_bond_setA,
                    i_bond_setB,
                    unmapped_atoms_molA,
                    newMapingSize,
                    new_i_bond_setA,
                    new_c_bond_setA,
                    new_Mapping,
                    counter);

            c_setA_copy = queryProcess.getCTab1();
            c_setB_copy = queryProcess.getCTab2();
            newSetBondNumA = queryProcess.getBondNumA();
            newNeighborNumA = queryProcess.getNeighborBondNumA();
            new_i_neighborsA = queryProcess.getIBondNeighboursA();
            new_c_neighborsA = queryProcess.getCBondNeighborsA();

            //find unmapped atoms of molecule B

            List<Integer> unmapped_atoms_molB = new ArrayList<Integer>();
            int unmapped_numB = 0;
            boolean atomB_is_unmapped = true;

            for (int a = 0; a < target.getAtomCount(); a++) {
                for (int b = 0; b < newMapingSize; b++) {
                    if (a == new_Mapping.get(b * 2 + 1)) {
                        atomB_is_unmapped = false;
                    }
                }
                if (atomB_is_unmapped) {
                    unmapped_atoms_molB.add(unmapped_numB++, a);
                }
                atomB_is_unmapped = true;
            }


            //number of remaining molecule B bonds after the clique search, which aren't neighbors
            int newSetBondNumB = 0; //instead of setNumB
            int newNeighborNumB = 0; //instead of localNeighborBondNumB

            TargetProcessor targetProcess = new TargetProcessor(
                    c_setA_copy,
                    c_setB_copy,
                    SignArray,
                    newNeighborNumB,
                    newSetBondNumB,
                    new_i_neighborsB,
                    new_c_neighborsB,
                    newNeighborNumA,
                    new_i_neighborsA,
                    new_c_neighborsA);


            targetProcess.process(
                    setNumB,
                    unmapped_atoms_molB,
                    newMapingSize,
                    i_bond_setB,
                    c_bond_setB,
                    new_Mapping,
                    counter,
                    new_i_bond_setB,
                    new_c_bond_setB);

            c_setA_copy = targetProcess.getCTab1();
            c_setB_copy = targetProcess.getCTab2();
            newSetBondNumB = targetProcess.getBondNumB();
            newNeighborNumB = targetProcess.getNeighborBondNumB();
            new_i_neighborsB = targetProcess.getIBondNeighboursB();
            new_c_neighborsB = targetProcess.getCBondNeighborsB();

//             System.out.println("Mapped Atoms before Iterator2: " + mappedAtoms);
            iterator(no_further_MAPPINGS,
                    newMapingSize,
                    new_Mapping,
                    newNeighborNumA,
                    newNeighborNumB,
                    new_i_neighborsA,
                    new_i_neighborsB,
                    new_c_neighborsA,
                    new_c_neighborsB,
                    newSetBondNumA,
                    newSetBondNumB,
                    new_i_bond_setA,
                    new_i_bond_setB,
                    new_c_bond_setA,
                    new_c_bond_setB);

            BESTARCS_copy.pop();
//            System.out.println("End of the iterator!!!!");
        }
    }

    private void setArcs(int xIndex, int yIndex, int arcsleft, List<Integer> TEMPMARCS, List<Integer> iBondNeighborAtomsA, List<Integer> iBondNeighborAtomsB, int neighborBondNumA, int neighborBondNumB) {
        do {
            yIndex++;
            if (yIndex == neighborBondNumB) {
                yIndex = 0;
                xIndex++;

            }
        } while ((xIndex < neighborBondNumA) && (TEMPMARCS.get(xIndex * neighborBondNumB + yIndex) != 1)); //Correction by ASAD set value minus 1
        if (xIndex < neighborBondNumA) {

            partsearch(xIndex, yIndex, TEMPMARCS, iBondNeighborAtomsA, iBondNeighborAtomsB, neighborBondNumA, neighborBondNumB);
            TEMPMARCS.set(xIndex * neighborBondNumB + yIndex, 0);
            partsearch(xIndex, yIndex, TEMPMARCS, iBondNeighborAtomsA, iBondNeighborAtomsB, neighborBondNumA, neighborBondNumB);

        } else {
            popBestArcs(arcsleft);


            if (checkMARCS(TEMPMARCS, neighborBondNumA, neighborBondNumB)) {
                bestARCS.push(TEMPMARCS);
            }
        }
    }

    private void popBestArcs(int arcsleft) {
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
    }
}
