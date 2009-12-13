
/* Copyright (C) 2005-2006  Markus Leber
 *               2006-2009  Syed Asad Rahman {asad@ebi.atomContainer.uk}
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

import java.util.List;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 * @cdk.module smsd
 */
public class TargetProcessor {

    private List<String> c_tab1_copy;
    private List<String> c_tab2_copy;
    private String[] SignArray;
    //number of remaining molecule A bonds after the clique search, which are
    //neighbors of the MCS
    private int neighborBondNumB = 0;
    //number of remaining molecule A bonds after the clique search, which aren't
    //neighbors
    private int setBondNumB = 0;
    private List<Integer> iBondNeighborsB;
    private List<String> cBondNeighborsB;

    /**
     *
     * @param c_tab1_copy
     * @param c_tab2_copy
     * @param SignArray
     * @param neighbor_bondnum_B
     * @param set_bondnum_B
     * @param i_bond_neighborsB
     * @param c_bond_neighborsB
     */
    protected TargetProcessor(
            List<String> c_tab1_copy,
            List<String> c_tab2_copy,
            String[] SignArray,
            int neighbor_bondnum_B,
            int set_bondnum_B,
            List<Integer> i_bond_neighborsB,
            List<String> c_bond_neighborsB) {

        this.c_tab1_copy = c_tab1_copy;
        this.c_tab2_copy = c_tab2_copy;
        this.SignArray = SignArray;
        this.neighborBondNumB = neighbor_bondnum_B;
        this.setBondNumB = set_bondnum_B;
        this.iBondNeighborsB = i_bond_neighborsB;
        this.cBondNeighborsB = c_bond_neighborsB;
    }

    protected void process(
            IAtomContainer target,
            List<Integer> unmapped_atoms_molB,
            int mappingSize,
            List<Integer> i_bond_setB,
            List<String> c_bond_setB,
            List<Integer> mapped_atoms,
            int counter,
            int neighborBondNumA,
            List<Integer> i_bond_neighborsA,
            List<String> c_bond_neighborsA) {


        int unmapped_numB = unmapped_atoms_molB.size();
        boolean bond_considered = false;
        boolean normal_bond = true;


        for (int a = 0; a < target.getBondCount(); a++) {

            Integer indexI = target.getAtomNumber(target.getBond(a).getAtom(0));
            Integer indexJ = target.getAtomNumber(target.getBond(a).getAtom(1));
            Integer order = target.getBond(a).getOrder().ordinal() + 1;

            for (int b = 0; b < unmapped_numB; b++) {
                if (unmapped_atoms_molB.get(b).equals(indexI)) {
                    normal_bond = unMappedAtomsEqualsIndexI(target, neighborBondNumA, mappingSize, a, c_bond_neighborsA, i_bond_neighborsA, counter, mapped_atoms, indexI, indexJ, order);
                    bond_considered = true;
                } else if (unmapped_atoms_molB.get(b) == indexJ) {
                    normal_bond = unMappedAtomsEqualsIndexJ(target, neighborBondNumA, mappingSize, a, c_bond_neighborsA, i_bond_neighborsA, counter, mapped_atoms, indexI, indexJ, order);
                    bond_considered = true;
                }

                if (normal_bond && bond_considered) {
                    markNormalBonds(a, i_bond_setB, c_bond_setB, indexI, indexJ, order);
                    normal_bond = true;
                    break;
                }

            }
            bond_considered = false;
        }

    }

    /**
     *
     * @param setNumB
     * @param unmapped_atoms_molB
     * @param newMapingSize
     * @param i_bond_setB
     * @param c_bond_setB
     * @param new_Mapping
     * @param counter
     * @param new_i_bond_setB
     * @param new_c_bond_setB
     * @param newNeighborNumA
     * @param new_i_neighborsA
     * @param new_c_neighborsA
     */
    protected void process(
            int setNumB,
            List<Integer> unmapped_atoms_molB,
            int newMapingSize,
            List<Integer> i_bond_setB,
            List<String> c_bond_setB,
            List<Integer> new_Mapping,
            int counter,
            List<Integer> new_i_bond_setB,
            List<String> new_c_bond_setB,
            int newNeighborNumA,
            List<Integer> new_i_neighborsA,
            List<String> new_c_neighborsA) {

        //The special signs must be transfered to the corresponding atoms of molecule A

        boolean bond_considered = false;
        boolean normal_bond = true;
        for (int a = 0; a < setNumB; a++) {

            Integer indexI = i_bond_setB.get(a * 3 + 0);
            Integer indexJ = i_bond_setB.get(a * 3 + 1);
            Integer order = i_bond_setB.get(a * 3 + 2);

            for (Integer unMappedAtomIndex : unmapped_atoms_molB) {
                if (unMappedAtomIndex.equals(indexI)) {
                    normal_bond = unMappedAtomsEqualsIndexI(setNumB, i_bond_setB, newNeighborNumA, newMapingSize, a, new_c_neighborsA, new_i_neighborsA, counter, new_Mapping, indexI, indexJ, order);
                    bond_considered = true;
                } else if (unMappedAtomIndex.equals(indexJ)) {
                    normal_bond = unMappedAtomsEqualsIndexJ(setNumB, i_bond_setB, newNeighborNumA, newMapingSize, a, new_c_neighborsA, new_i_neighborsA, counter, new_Mapping, indexI, indexJ, order);
                    bond_considered = true;
                }
                if (normal_bond && bond_considered) {
                    markNormalBonds(a, new_i_bond_setB, new_c_bond_setB, indexI, indexJ, order);
                    normal_bond = true;
                    break;
                }

            }
            bond_considered = false;
        }
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

    private boolean unMappedAtomsEqualsIndexI(
            IAtomContainer target,
            int neighborBondNumA,
            int mappingSize,
            int a,
            List<String> c_bond_neighborsA,
            List<Integer> i_bond_neighborsA,
            int counter,
            List<Integer> mapped_atoms,
            Integer indexI,
            Integer indexJ,
            Integer order) {
        boolean normal_bond = true;
        for (int c = 0; c < mappingSize; c++) {
            if (mapped_atoms.get(c * 2 + 1).equals(indexJ)) {
                setIBondNeighbors(indexI, indexJ, order);
                if (c_tab2_copy.get(a * 4 + 3).compareToIgnoreCase("X") == 0) {
                    step1(a, counter);
                    changeCharBonds(indexJ, SignArray[counter], target.getBondCount(), target, c_tab2_copy);
                    int cor_atom = McGregorChecks.searchCorrespondingAtom(mappingSize, indexJ, 2, mapped_atoms);
                    //Commented by Asad
                    changeCharBonds(cor_atom, SignArray[counter], neighborBondNumA, i_bond_neighborsA, c_bond_neighborsA);
//                                changeCharBonds(cor_atom, SignArray[counter], query.getBondCount(), query, c_tab1_copy);
                    counter++;
                } else {
                    step2(a);
                }
                normal_bond = false;
                neighborBondNumB++;
            }
        }
        return normal_bond;
    }

    private boolean unMappedAtomsEqualsIndexJ(
            IAtomContainer target,
            int neighborBondNumA,
            int mappingSize,
            int a,
            List<String> c_bond_neighborsA,
            List<Integer> i_bond_neighborsA,
            int counter, List<Integer> mapped_atoms,
            Integer indexI,
            Integer indexJ,
            Integer order) {
        boolean normal_bond = true;
        for (int c = 0; c < mappingSize; c++) {
            if (mapped_atoms.get(c * 2 + 1).equals(indexI)) {
                setIBondNeighbors(indexI, indexJ, order);
                if (c_tab2_copy.get(a * 4 + 2).compareToIgnoreCase("X") == 0) {
                    step3(a, counter);
                    changeCharBonds(indexI, SignArray[counter], target.getBondCount(), target, c_tab2_copy);
                    int cor_atom = McGregorChecks.searchCorrespondingAtom(mappingSize, indexI, 2, mapped_atoms);
                    changeCharBonds(cor_atom, SignArray[counter], neighborBondNumA, i_bond_neighborsA, c_bond_neighborsA);
//                                changeCharBonds(cor_atom, SignArray[counter], query.getBondCount(), query, c_tab1_copy);
                    counter++;
                } else {
                    step4(a);
                }
                normal_bond = false;
                neighborBondNumB++;
            }
        }

        return normal_bond;
    }

    private boolean unMappedAtomsEqualsIndexI(
            int setNumB,
            List<Integer> i_bond_setB,
            int newNeighborNumA,
            int newMappingSize,
            int a,
            List<String> new_c_neighborsA,
            List<Integer> new_i_neighborsA,
            int counter,
            List<Integer> new_Mapping,
            Integer indexI,
            Integer indexJ,
            Integer order) {
        boolean normal_bond = true;
        for (int c = 0; c < newMappingSize; c++) {
            if (new_Mapping.get(c * 2 + 1).equals(indexJ)) {
                setIBondNeighbors(indexI, indexJ, order);
                if (c_tab2_copy.get(a * 4 + 3).compareToIgnoreCase("X") == 0) {
                    step1(a, counter);
                    changeCharBonds(indexJ, SignArray[counter], setNumB, i_bond_setB, c_tab2_copy);
                    int cor_atom = McGregorChecks.searchCorrespondingAtom(newMappingSize, indexJ, 2, new_Mapping);
                    changeCharBonds(cor_atom, SignArray[counter], newNeighborNumA, new_i_neighborsA, new_c_neighborsA);
                    counter++;

                } else {
                    step2(a);
                }

                normal_bond = false;
                neighborBondNumB++;

            }
        }
        return normal_bond;
    }

    private boolean unMappedAtomsEqualsIndexJ(
            int setNumB,
            List<Integer> i_bond_setB,
            int newNeighborNumA,
            int newMappingSize,
            int a,
            List<String> new_c_neighborsA,
            List<Integer> new_i_neighborsA,
            int counter, List<Integer> new_Mapping,
            Integer indexI,
            Integer indexJ,
            Integer order) {
        boolean normal_bond = true;
        for (int c = 0; c < newMappingSize; c++) {
            if (new_Mapping.get(c * 2 + 1).equals(indexI)) {
                setIBondNeighbors(indexI, indexJ, order);

                if (c_tab2_copy.get(a * 4 + 2).compareToIgnoreCase("X") == 0) {

                    step3(a, counter);
                    changeCharBonds(indexI, SignArray[counter], setNumB, i_bond_setB, c_tab2_copy);
                    int cor_atom = McGregorChecks.searchCorrespondingAtom(newMappingSize, indexI, 2, new_Mapping);
                    changeCharBonds(cor_atom, SignArray[counter], newNeighborNumA, new_i_neighborsA, new_c_neighborsA);
                    counter++;
                } else {
                    step4(a);
                }

                normal_bond = false;
                neighborBondNumB++;

            }


        }

        return normal_bond;
    }

    private void markNormalBonds(
            int a,
            List<Integer> i_bond_setB,
            List<String> c_bond_setB,
            Integer indexI,
            Integer indexJ,
            Integer order) {
        i_bond_setB.add(indexI);
        i_bond_setB.add(indexJ);
        i_bond_setB.add(order);
        c_bond_setB.add(c_tab2_copy.get(a * 4 + 0));
        c_bond_setB.add(c_tab2_copy.get(a * 4 + 1));
        c_bond_setB.add("X");
        c_bond_setB.add("X");
        setBondNumB++;
    }

    private void setIBondNeighbors(Integer indexI,
            Integer indexJ,
            Integer order) {
        iBondNeighborsB.add(indexI);
        iBondNeighborsB.add(indexJ);
        iBondNeighborsB.add(order);
    }

    private void step1(int a, int counter) {
        cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
        cBondNeighborsB.add(SignArray[counter]);
        cBondNeighborsB.add("X");
        cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
    }

    private void step2(int a) {
        cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
        cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
        cBondNeighborsB.add("X");
        cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 3));
    }

    private void step3(int a, int counter) {
        cBondNeighborsB.add(SignArray[counter]);
        cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
        cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
        cBondNeighborsB.add("X");
    }

    private void step4(int a) {
        cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
        cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
        cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 2));
        cBondNeighborsB.add("X");
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
    protected int getNeighborBondNumB() {
        return this.neighborBondNumB;
    }

    /**
     *
     * @return number of remaining molecule A bonds after the clique search,
     * which aren't neighbors
     */
    protected int getBondNumB() {
        return this.setBondNumB;
    }

    List<Integer> getIBondNeighboursB() {
        return this.iBondNeighborsB;
    }

    List<String> getCBondNeighborsB() {
        return this.cBondNeighborsB;
    }
}
