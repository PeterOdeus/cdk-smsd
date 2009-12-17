
/* Copyright (C) 2005-2006 Markus Leber
 *               2006-2009 Syed Asad Rahman {asad@ebi.atomContainer.uk}
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

import java.util.ArrayList;
import java.util.List;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 * @cdk.module smsd
 */
public class QueryProcessor {

    private List<String> c_tab1_copy;
    private List<String> c_tab2_copy;
    private String[] SignArray;
    private int neighborBondNumA = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
    private int setBondNumA = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
    private List<Integer> iBondNeighborsA;
    private List<String> cBondNeighborsA;

    /**
     * @param c_tab1_copy
     * @param c_tab2_copy
     * @param SignArray
     * @param neighbor_bondnum_A
     * @param set_bondnum_A
     * @param i_bond_neighborsA
     * @param c_bond_neighborsA
     */
    protected QueryProcessor(
            List<String> c_tab1_copy,
            List<String> c_tab2_copy,
            String[] SignArray,
            int neighbor_bondnum_A,
            int set_bondnum_A,
            List<Integer> i_bond_neighborsA,
            List<String> c_bond_neighborsA) {

        this.c_tab1_copy = c_tab1_copy;
        this.c_tab2_copy = c_tab2_copy;
        this.SignArray = SignArray;
        this.neighborBondNumA = neighbor_bondnum_A;
        this.setBondNumA = set_bondnum_A;
        this.iBondNeighborsA = i_bond_neighborsA;
        this.cBondNeighborsA = c_bond_neighborsA;
    }

    /**
     *
     * @param query
     * @param target 
     * @param unmapped_atoms_molA
     * @param mappingSize
     * @param i_bond_setA
     * @param c_bond_setA
     * @param mapped_atoms
     * @param counter
     */
    protected void process(
            IAtomContainer query,
            IAtomContainer target,
            List<Integer> unmapped_atoms_molA,
            int mappingSize,
            List<Integer> i_bond_setA,
            List<String> c_bond_setA,
            List<Integer> mapped_atoms,
            int counter) {

        int unmapped_numA = unmapped_atoms_molA.size();
        boolean bond_considered = false;
        boolean normal_bond = true;

//        System.out.println("\n" + c_tab1_copy + "\n");


        for (int atomIndex = 0; atomIndex < query.getBondCount(); atomIndex++) {


            Integer indexI = query.getAtomNumber(query.getBond(atomIndex).getAtom(0));
            Integer indexJ = query.getAtomNumber(query.getBond(atomIndex).getAtom(1));
            Integer order = query.getBond(atomIndex).getOrder().ordinal() + 1;

//            System.out.println(AtomI + "= , =" + AtomJ );
            for (Integer unMappedAtomIndex = 0; unMappedAtomIndex < unmapped_numA; unMappedAtomIndex++) {

                if (unmapped_atoms_molA.get(unMappedAtomIndex).equals(indexI)) {
                    normal_bond = unMappedAtomsEqualsIndexI(query, target, mappingSize, atomIndex, counter, mapped_atoms, indexI, indexJ, order);
                    bond_considered = true;
                } else //Does a ungemaptes atom at second position in the connection occur?
                if (unmapped_atoms_molA.get(unMappedAtomIndex).equals(indexJ)) {
                    normal_bond = unMappedAtomsEqualsIndexJ(query, target, mappingSize, atomIndex, counter, mapped_atoms, indexI, indexJ, order);
                    bond_considered = true;
                }
                if (normal_bond && bond_considered) {
                    markNormalBonds(atomIndex, i_bond_setA, c_bond_setA, indexI, indexJ, order);
                    normal_bond = true;
                    break;
                }
            }
            bond_considered = false;
        }

        /*******************************************************************************///
    }

    /**
     * 
     * @param setNumA
     * @param setNumB
     * @param i_bond_setA
     * @param i_bond_setB
     * @param unmapped_atoms_molA
     * @param newMapingSize
     * @param new_i_bond_setA
     * @param new_c_bond_setA
     * @param new_Mapping
     * @param counter
     */
    protected void process(
            int setNumA,
            int setNumB,
            List<Integer> i_bond_setA,
            List<Integer> i_bond_setB,
            List<Integer> unmapped_atoms_molA,
            int newMapingSize,
            List<Integer> new_i_bond_setA,
            List<String> new_c_bond_setA,
            List<Integer> new_Mapping,
            int counter) {


        boolean bond_considered = false;
        boolean normal_bond = true;


        for (int atomIndex = 0; atomIndex < setNumA; atomIndex++) {
            Integer indexI = i_bond_setA.get(atomIndex * 3 + 0);
            Integer indexJ = i_bond_setA.get(atomIndex * 3 + 1);
            Integer order = i_bond_setA.get(atomIndex * 3 + 2);

            for (Integer unMappedAtomIndex : unmapped_atoms_molA) {
                if (unMappedAtomIndex.equals(indexI)) {
                    normal_bond = UnMappedAtomsEqualsIndexI(setNumA, setNumB, i_bond_setA, i_bond_setB, newMapingSize, atomIndex, counter, new_Mapping, indexI, indexJ, order);
                    bond_considered = true;
                } else if (unMappedAtomIndex.equals(indexJ)) {
                    normal_bond = UnMappedAtomsEqualsIndexJ(setNumA, setNumB, i_bond_setA, i_bond_setB, newMapingSize, atomIndex, counter, new_Mapping, indexI, indexJ, order);
                    bond_considered = true;
                }

                if (normal_bond && bond_considered) {
                    markNormalBonds(atomIndex, new_i_bond_setA, new_c_bond_setA, indexI, indexJ, order);
                    normal_bond = true;
                    break;
                }
            }
            bond_considered = false;
        }
    }

    private int searchCorrespondingAtom(int mapped_atoms_size, int atom_from_other_molecule, int molecule, List<Integer> mapped_atoms_org) {


        List<Integer> mapped_atoms = new ArrayList<Integer>(mapped_atoms_org);

        int corresponding_atom = 0;
        for (int a = 0; a < mapped_atoms_size; a++) {
            if ((molecule == 1)
                    && (mapped_atoms.get(a * 2 + 0).intValue() == atom_from_other_molecule)) {

                corresponding_atom = mapped_atoms.get(a * 2 + 1);


            }
            if ((molecule == 2)
                    && (mapped_atoms.get(a * 2 + 1).intValue() == atom_from_other_molecule)) {
                corresponding_atom = mapped_atoms.get(a * 2 + 0);
            }
        }
        return corresponding_atom;
    }

    private int changeCharBonds(int corresponding_atom, String new_symbol, int neighbor_bondnum, IAtomContainer atomContainer, List<String> c_bond_neighbors) {
        for (int a = 0; a < neighbor_bondnum; a++) {
            IBond bond = atomContainer.getBond(a);
            if ((atomContainer.getAtomNumber(bond.getAtom(0)) == corresponding_atom) && (c_bond_neighbors.get(a * 4 + 2).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 2, c_bond_neighbors.get(a * 4 + 0));
                c_bond_neighbors.set(a * 4 + 0, new_symbol);
            }

            if ((atomContainer.getAtomNumber(bond.getAtom(1)) == corresponding_atom) && (c_bond_neighbors.get(a * 4 + 3).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 3, c_bond_neighbors.get(a * 4 + 1));
                c_bond_neighbors.set(a * 4 + 1, new_symbol);
            }

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

    private void markNormalBonds(int atomIndex,
            List<Integer> i_bond_setA,
            List<String> c_bond_setA,
            Integer indexI,
            Integer indexJ,
            Integer order) {
        i_bond_setA.add(indexI);
        i_bond_setA.add(indexJ);
        i_bond_setA.add(order);
        c_bond_setA.add(c_tab1_copy.get(atomIndex * 4 + 0));
        c_bond_setA.add(c_tab1_copy.get(atomIndex * 4 + 1));
        c_bond_setA.add("X");
        c_bond_setA.add("X");
        setBondNumA++;
    }

    private void step1(int atomIndex, int counter) {
        cBondNeighborsA.add(c_tab1_copy.get(atomIndex * 4 + 0));
        cBondNeighborsA.add(SignArray[counter]);
        cBondNeighborsA.add("X");
        cBondNeighborsA.add(c_tab1_copy.get(atomIndex * 4 + 1));
    }

    private void step2(int atomIndex) {
        cBondNeighborsA.add(c_tab1_copy.get(atomIndex * 4 + 0));
        cBondNeighborsA.add(c_tab1_copy.get(atomIndex * 4 + 1));
        cBondNeighborsA.add("X");
        cBondNeighborsA.add(c_tab1_copy.get(atomIndex * 4 + 3));
    }

    private void step3(int atomIndex, int counter) {
        cBondNeighborsA.add(SignArray[counter]);
        cBondNeighborsA.add(c_tab1_copy.get(atomIndex * 4 + 1));
        cBondNeighborsA.add(c_tab1_copy.get(atomIndex * 4 + 0));
        cBondNeighborsA.add("X");
    }

    private void step4(int atomIndex) {
        cBondNeighborsA.add(c_tab1_copy.get(atomIndex * 4 + 0));
        cBondNeighborsA.add(c_tab1_copy.get(atomIndex * 4 + 1));
        cBondNeighborsA.add(c_tab1_copy.get(atomIndex * 4 + 2));
        cBondNeighborsA.add("X");
    }

    private boolean unMappedAtomsEqualsIndexI(
            IAtomContainer query,
            IAtomContainer target,
            int mappingSize,
            int atomIndex,
            int counter,
            List<Integer> mapped_atoms,
            Integer indexI,
            Integer indexJ,
            Integer order) {
        boolean normal_bond = true;
        for (int c = 0; c < mappingSize; c++) {

            if (mapped_atoms.get(c * 2).equals(indexJ)) {
                setBondNeighbors(indexI, indexJ, order);
                if (c_tab1_copy.get(atomIndex * 4 + 3).compareToIgnoreCase("X") == 0) {

                    step1(atomIndex, counter);
                    changeCharBonds(indexJ, SignArray[counter], query.getBondCount(), query, c_tab1_copy);

                    int cor_atom = searchCorrespondingAtom(mappingSize, indexJ, 1, mapped_atoms);
                    changeCharBonds(cor_atom, SignArray[counter], target.getBondCount(), target, c_tab2_copy);
                    counter++;
                } else {
                    step2(atomIndex);
                }
                normal_bond = false;
                neighborBondNumA++;
            }
        }
        return normal_bond;
    }

    private boolean unMappedAtomsEqualsIndexJ(
            IAtomContainer query,
            IAtomContainer target,
            int mappingSize,
            int atomIndex,
            int counter,
            List<Integer> mapped_atoms,
            Integer indexI,
            Integer indexJ,
            Integer order) {

        boolean normal_bond = true;
        for (int c = 0; c < mappingSize; c++) {


            if (mapped_atoms.get(c * 2 + 0).equals(indexI)) {
                setBondNeighbors(indexI, indexJ, order);
                if (c_tab1_copy.get(atomIndex * 4 + 2).compareToIgnoreCase("X") == 0) {
                    step3(atomIndex, counter);
                    changeCharBonds(indexI, SignArray[counter], query.getBondCount(), query, c_tab1_copy);

                    int cor_atom = searchCorrespondingAtom(mappingSize, indexI, 1, mapped_atoms);
                    changeCharBonds(cor_atom, SignArray[counter], target.getBondCount(), target, c_tab2_copy);
                    counter++;
                } else {
                    step4(atomIndex);
                }
                normal_bond = false;
                neighborBondNumA++;
                //System.out.println("Neighbor");
                //System.out.println(neighborBondNumA);
            }
        }
        return normal_bond;
    }

    private boolean UnMappedAtomsEqualsIndexI(
            int setNumA,
            int setNumB,
            List<Integer> i_bond_setA,
            List<Integer> i_bond_setB,
            int newMapingSize,
            int atomIndex,
            int counter,
            List<Integer> new_Mapping,
            Integer indexI,
            Integer indexJ,
            Integer order) {
        boolean normal_bond = true;
        for (int c = 0; c < newMapingSize; c++) {

            if (new_Mapping.get(c * 2 + 0).equals(indexJ)) {

                setBondNeighbors(indexI, indexJ, order);
                if (c_tab1_copy.get(atomIndex * 4 + 3).compareToIgnoreCase("X") == 0) {
                    step1(atomIndex, counter);
                    changeCharBonds(indexJ, SignArray[counter], setNumA, i_bond_setA, c_tab1_copy);
                    int cor_atom = McGregorChecks.searchCorrespondingAtom(newMapingSize, indexJ, 1, new_Mapping);
                    changeCharBonds(cor_atom, SignArray[counter], setNumB, i_bond_setB, c_tab2_copy);
                    counter++;

                } else {
                    step2(atomIndex);
                }
                normal_bond = false;
                neighborBondNumA++;

            }
        }
        return normal_bond;
    }

    private boolean UnMappedAtomsEqualsIndexJ(
            int setNumA,
            int setNumB,
            List<Integer> i_bond_setA,
            List<Integer> i_bond_setB,
            int newMappingSize,
            int atomIndex,
            int counter, List<Integer> new_Mapping,
            Integer indexI,
            Integer indexJ,
            Integer order) {
        boolean normal_bond = true;
        for (int c = 0; c < newMappingSize; c++) {

            if (new_Mapping.get(c * 2 + 0).equals(indexI)) {

                setBondNeighbors(indexI, indexJ, order);
                if (c_tab1_copy.get(atomIndex * 4 + 2).compareToIgnoreCase("X") == 0) {
                    step3(atomIndex, counter);
                    changeCharBonds(indexI, SignArray[counter], setNumA, i_bond_setA, c_tab1_copy);
                    int cor_atom = McGregorChecks.searchCorrespondingAtom(newMappingSize, indexI, 1, new_Mapping);
                    changeCharBonds(cor_atom, SignArray[counter], setNumB, i_bond_setB, c_tab2_copy);
                    counter++;
                } else {
                    step4(atomIndex);
                }

                normal_bond = false;
                neighborBondNumA++;

            }
        }
        return normal_bond;
    }

    private void setBondNeighbors(Integer indexI,
            Integer indexJ,
            Integer order) {
        iBondNeighborsA.add(indexI);
        iBondNeighborsA.add(indexJ);
        iBondNeighborsA.add(order);
    }

    /**
     *
     * @return
     */
    protected List<String> getCTab1() {
        return this.c_tab1_copy;
    }

    /**
     *
     * @return
     */
    protected List<String> getCTab2() {
        return this.c_tab2_copy;
    }

    /**
     *
     * @return number of remaining molecule A bonds after the clique search,
     * which are neighbors of the MCS
     *
     */
    protected int getNeighborBondNumA() {
        return this.neighborBondNumA;
    }

    /**
     *
     * @return number of remaining molecule A bonds after the clique search,
     * which aren't neighbors
     */
    protected int getBondNumA() {
        return this.setBondNumA;
    }

    List<Integer> getIBondNeighboursA() {
        return this.iBondNeighborsA;
    }

    List<String> getCBondNeighborsA() {
        return this.cBondNeighborsA;
    }
}
