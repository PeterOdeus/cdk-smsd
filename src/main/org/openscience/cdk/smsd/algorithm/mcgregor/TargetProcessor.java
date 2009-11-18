/* Copyright (C) 2005-2006  Markus Leber
 *               2006-2009  Syed Asad Rahman {asad@ebi.ac.uk}
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
    private String[] SignROW;
    private int neighbor_bondnum_B = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
    private int set_bondnum_B = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
    private IAtomContainer target;
    private int neighbor_bondnum_A;
    private List<Integer> i_bond_neighborsA;
    private List<String> c_bond_neighborsA;

    /**
     * @param target
     * @param c_tab1_copy
     * @param c_tab2_copy
     * @param SignROW
     * @param neighbor_bondnum_B
     * @param set_bondnum_B
     * @param neighbor_bondnum_A
     * @param i_bond_neighborsA
     * @param c_bond_neighborsA 
     */
    protected TargetProcessor(IAtomContainer target, List<String> c_tab1_copy, List<String> c_tab2_copy, String[] SignROW, int neighbor_bondnum_B, int set_bondnum_B, int neighbor_bondnum_A, List<Integer> i_bond_neighborsA, List<String> c_bond_neighborsA) {

        this.target = target;

        this.c_tab1_copy = c_tab1_copy;
        this.c_tab2_copy = c_tab2_copy;
        this.SignROW = SignROW;
        this.neighbor_bondnum_B = neighbor_bondnum_B;
        this.set_bondnum_B = set_bondnum_B;
        this.c_bond_neighborsA = c_bond_neighborsA;
        this.i_bond_neighborsA = i_bond_neighborsA;
        this.neighbor_bondnum_A = neighbor_bondnum_A;
    }

    protected void process(
            List<Integer> unmapped_atoms_molB,
            int mappingSize,
            List<Integer> i_bond_neighborsB,
            List<Integer> i_bond_setB,
            List<String> c_bond_neighborsB,
            List<String> c_bond_setB,
            List<Integer> mapped_atoms,
            int SR_count) {

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
                            i_bond_neighborsB.add(indexI);
                            i_bond_neighborsB.add(indexJ);
                            i_bond_neighborsB.add(order);
                            if (c_tab2_copy.get(a * 4 + 3).compareToIgnoreCase("X") == 0) {
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add(SignROW[SR_count]);
                                c_bond_neighborsB.add("X");
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                change_char_bonds(indexJ, SignROW[SR_count], target.getBondCount(), target, c_tab2_copy);
                                int cor_atom = search_corresponding_atom(mappingSize, indexJ, 2, mapped_atoms);
                                //Commented by Asad
                                change_char_bonds(cor_atom, SignROW[SR_count], neighbor_bondnum_A, i_bond_neighborsA, c_bond_neighborsA);
//                                change_char_bonds(cor_atom, SignROW[SR_count], query.getBondCount(), query, c_tab1_copy);
                                SR_count++;
                            } else {
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                c_bond_neighborsB.add("X");
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 3));
                            }
                            normal_bond = false;
                            neighbor_bondnum_B++;
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
                        set_bondnum_B++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                if (unmapped_atoms_molB.get(b) == indexJ) {
                    for (int c = 0; c < mappingSize; c++) {
                        if (mapped_atoms.get(c * 2 + 1).equals(indexI)) {
                            i_bond_neighborsB.add(indexI);
                            i_bond_neighborsB.add(indexJ);
                            i_bond_neighborsB.add(order);
                            if (c_tab2_copy.get(a * 4 + 2).compareToIgnoreCase("X") == 0) {
                                c_bond_neighborsB.add(SignROW[SR_count]);
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add("X");
                                change_char_bonds(indexI, SignROW[SR_count], target.getBondCount(), target, c_tab2_copy);
                                int cor_atom = search_corresponding_atom(mappingSize, indexI, 2, mapped_atoms);
                                change_char_bonds(cor_atom, SignROW[SR_count], neighbor_bondnum_A, i_bond_neighborsA, c_bond_neighborsA);
//                                change_char_bonds(cor_atom, SignROW[SR_count], query.getBondCount(), query, c_tab1_copy);
                                SR_count++;
                            } else {
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 0));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 1));
                                c_bond_neighborsB.add(c_tab2_copy.get(a * 4 + 2));
                                c_bond_neighborsB.add("X");
                            }
                            normal_bond = false;
                            neighbor_bondnum_B++;
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

    }

    private int search_corresponding_atom(int mapped_atoms_size, int atom_from_other_molecule, int molecule, List<Integer> mapped_atoms_org) {


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

    private int change_char_bonds(int corresponding_atom, String new_symbol, int neighbor_bondnum, IAtomContainer ac, List<String> c_bond_neighbors) {
        //private int change_char_bonds(int corresponding_atom, String new_symbol, int neighbor_bondnum, Vector<Integer> i_bond_neighbors, Vector<String> c_bond_neighbors) {

        for (int a = 0; a < neighbor_bondnum; a++) {
            IBond bond = ac.getBond(a);
            if ((ac.getAtomNumber(bond.getAtom(0)) == corresponding_atom) && (c_bond_neighbors.get(a * 4 + 2).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 2, c_bond_neighbors.get(a * 4 + 0));
                c_bond_neighbors.set(a * 4 + 0, new_symbol);
            }

            if ((ac.getAtomNumber(bond.getAtom(1)) == corresponding_atom) && (c_bond_neighbors.get(a * 4 + 3).compareToIgnoreCase("X") == 0)) {
                c_bond_neighbors.set(a * 4 + 3, c_bond_neighbors.get(a * 4 + 1));
                c_bond_neighbors.set(a * 4 + 1, new_symbol);
            }

        }

        return 0;
    }

    private int change_char_bonds(int corresponding_atom, String new_symbol, int neighbor_bondnum, List<Integer> i_bond_neighbors, List<String> c_bond_neighbors) {

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
    public List<String> getCtab1() {
        return this.c_tab1_copy;
    }

    /**
     *
     * @return
     */
    public List<String> getCtab2() {
        return this.c_tab2_copy;
    }

    /**
     *
     * @return number of remaining molecule A bonds after the clique search,
     * which are neighbors of the MCS
     *
     */
    public int getNeighbor_bondnum_B() {
        return this.neighbor_bondnum_B;
    }

    /**
     *
     * @return number of remaining molecule A bonds after the clique search,
     * which aren't neighbors
     */
    public int geBondnum_B() {
        return this.set_bondnum_B;
    }

    public String[] getSigns() {
        return this.SignROW;
    }
}
