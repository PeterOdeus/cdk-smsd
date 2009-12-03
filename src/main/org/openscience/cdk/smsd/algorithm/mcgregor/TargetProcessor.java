
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

import java.util.ArrayList;
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
    private int neighborBondNumB = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
    private int setBondNumB = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
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
            String AtomI = target.getBond(a).getAtom(0).getSymbol();
            String AtomJ = target.getBond(a).getAtom(1).getSymbol();
            Integer order = target.getBond(a).getOrder().ordinal() + 1;

            for (int b = 0; b < unmapped_numB; b++) {
                if (unmapped_atoms_molB.get(b).equals(indexI)) {
                    for (int c = 0; c < mappingSize; c++) {
                        if (mapped_atoms.get(c * 2 + 1).equals(indexJ)) {
                            iBondNeighborsB.add(indexI);
                            iBondNeighborsB.add(indexJ);
                            iBondNeighborsB.add(order);
                            if (c_tab2_copy.get(a * 4 + 3).compareToIgnoreCase("X") == 0) {
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                cBondNeighborsB.add(SignArray[counter]);
                                cBondNeighborsB.add("X");
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                changeCharBonds(indexJ, SignArray[counter], target.getBondCount(), target, c_tab2_copy);
                                int cor_atom = searchCorrespondingAtom(mappingSize, indexJ, 2, mapped_atoms);
                                //Commented by Asad
                                changeCharBonds(cor_atom, SignArray[counter], neighborBondNumA, i_bond_neighborsA, c_bond_neighborsA);
//                                changeCharBonds(cor_atom, SignArray[counter], query.getBondCount(), query, c_tab1_copy);
                                counter++;
                            } else {
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                cBondNeighborsB.add("X");
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 3));
                            }
                            normal_bond = false;
                            neighborBondNumB++;
                        }
                    }
                    if (normal_bond) {
                        i_bond_setB.add(indexI);
                        i_bond_setB.add(indexJ);
                        i_bond_setB.add(order);
                        c_bond_setB.add(c_tab2_copy.get(a * 4 + 0));
                        c_bond_setB.add(c_tab2_copy.get(a * 4 + 1));
                        c_bond_setB.add("X");
                        c_bond_setB.add("X");
                        setBondNumB++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                if (unmapped_atoms_molB.get(b) == indexJ) {
                    for (int c = 0; c < mappingSize; c++) {
                        if (mapped_atoms.get(c * 2 + 1).equals(indexI)) {
                            iBondNeighborsB.add(indexI);
                            iBondNeighborsB.add(indexJ);
                            iBondNeighborsB.add(order);
                            if (c_tab2_copy.get(a * 4 + 2).compareToIgnoreCase("X") == 0) {
                                cBondNeighborsB.add(SignArray[counter]);
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                cBondNeighborsB.add("X");
                                changeCharBonds(indexI, SignArray[counter], target.getBondCount(), target, c_tab2_copy);
                                int cor_atom = searchCorrespondingAtom(mappingSize, indexI, 2, mapped_atoms);
                                changeCharBonds(cor_atom, SignArray[counter], neighborBondNumA, i_bond_neighborsA, c_bond_neighborsA);
//                                changeCharBonds(cor_atom, SignArray[counter], query.getBondCount(), query, c_tab1_copy);
                                counter++;
                            } else {
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 2));
                                cBondNeighborsB.add("X");
                            }
                            normal_bond = false;
                            neighborBondNumB++;
                        }
                    }
                    if (normal_bond) {

//                        System.out.println("\n\n normal_bond:" + normal_bond);

                        i_bond_setB.add(indexI);
                        i_bond_setB.add(indexJ);
                        i_bond_setB.add(order);
                        c_bond_setB.add(AtomI);
                        c_bond_setB.add(AtomJ);
                        c_bond_setB.add("X");
                        c_bond_setB.add("X");
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
            for (Integer unMappedAtomIndex : unmapped_atoms_molB) {
                if (unMappedAtomIndex.equals(i_bond_setB.get(a * 3 + 0))) {
                    for (int c = 0; c < newMapingSize; c++) {
                        if (new_Mapping.get(c * 2 + 1).equals(i_bond_setB.get(a * 3 + 1))) {
                            iBondNeighborsB.add(i_bond_setB.get(a * 3 + 0));
                            iBondNeighborsB.add(i_bond_setB.get(a * 3 + 1));
                            iBondNeighborsB.add(i_bond_setB.get(a * 3 + 2));
                            if (c_tab2_copy.get(a * 4 + 3).compareToIgnoreCase("X") == 0) {
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                cBondNeighborsB.add(SignArray[counter]);
                                cBondNeighborsB.add("X");
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                changeCharBonds(i_bond_setB.get(a * 3 + 1), SignArray[counter], setNumB, i_bond_setB, c_tab2_copy);
                                int cor_atom = McGregorChecks.searchCorrespondingAtom(newMapingSize, i_bond_setB.get(a * 3 + 1), 2, new_Mapping);
                                changeCharBonds(cor_atom, SignArray[counter], newNeighborNumA, new_i_neighborsA, new_c_neighborsA);
                                counter++;

                            } else {
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                cBondNeighborsB.add("X");
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 3));
                            }

                            normal_bond = false;
                            neighborBondNumB++;

                        }
                    }
                    if (normal_bond) {

                        new_i_bond_setB.add(i_bond_setB.get(a * 3 + 0));
                        new_i_bond_setB.add(i_bond_setB.get(a * 3 + 1));
                        new_i_bond_setB.add(i_bond_setB.get(a * 3 + 2));
                        new_c_bond_setB.add(c_tab2_copy.get(a * 4 + 0));
                        new_c_bond_setB.add(c_tab2_copy.get(a * 4 + 1));
                        new_c_bond_setB.add("X");
                        new_c_bond_setB.add("X");
                        setBondNumB++;
                    }


                    normal_bond = true;
                    bond_considered = true;
                }
                if (unMappedAtomIndex.equals(i_bond_setB.get(a * 3 + 1))) {
                    for (int c = 0; c < newMapingSize; c++) {

                        if (new_Mapping.get(c * 2 + 1).equals(i_bond_setB.get(a * 3 + 0))) {

                            iBondNeighborsB.add(i_bond_setB.get(a * 3 + 0));
                            iBondNeighborsB.add(i_bond_setB.get(a * 3 + 1));
                            iBondNeighborsB.add(i_bond_setB.get(a * 3 + 2));

                            if (c_tab2_copy.get(a * 4 + 2).compareToIgnoreCase("X") == 0) {

                                cBondNeighborsB.add(SignArray[counter]);
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                cBondNeighborsB.add("X");
                                changeCharBonds(i_bond_setB.get(a * 3 + 0), SignArray[counter], setNumB, i_bond_setB, c_tab2_copy);
                                int cor_atom = McGregorChecks.searchCorrespondingAtom(newMapingSize, i_bond_setB.get(a * 3 + 0), 2, new_Mapping);
                                changeCharBonds(cor_atom, SignArray[counter], newNeighborNumA, new_i_neighborsA, new_c_neighborsA);
                                counter++;
                            } else {
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                cBondNeighborsB.add(c_tab2_copy.get(a * 4 + 2));
                                cBondNeighborsB.add("X");
                            }

                            normal_bond = false;
                            neighborBondNumB++;

                        }


                    }

                    if (normal_bond) {

                        new_i_bond_setB.add(i_bond_setB.get(a * 3 + 0));
                        new_i_bond_setB.add(i_bond_setB.get(a * 3 + 1));
                        new_i_bond_setB.add(i_bond_setB.get(a * 3 + 2));
                        new_c_bond_setB.add(c_tab2_copy.get(a * 4 + 0));
                        new_c_bond_setB.add(c_tab2_copy.get(a * 4 + 1));
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
    }

    private int searchCorrespondingAtom(int mapped_atoms_size, int atom_from_other_molecule, int molecule, List<Integer> mapped_atoms_org) {


        List<Integer> mapped_atoms = new ArrayList<Integer>(mapped_atoms_org);

        int corresponding_atom = 0;
        for (int a = 0; a < mapped_atoms_size; a++) {
            if ((molecule == 1) &&
                    (mapped_atoms.get(a * 2 + 0).intValue() == atom_from_other_molecule)) {
                corresponding_atom = mapped_atoms.get(a * 2 + 1);
            }
            if ((molecule == 2) && (mapped_atoms.get(a * 2 + 1).intValue() == atom_from_other_molecule)) {
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
