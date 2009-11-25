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

    private int neighborBondNumB = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
    private int setBondNumB = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
    private int neighborBondNumA = 0;
    private List<Integer> i_bond_neighborsA = null;
    private List<String> c_bond_neighborsA = null;
    private List<String> cBondNeighborsB = null;
    private List<Integer> iBondSetB = null;
    private List<String> cBondSetB = null;

    /**
     *
     * @param c_bond_neighborsB
     * @param iBondSetB
     * @param c_bond_setB
     * @param neighbor_bondnum_B
     * @param set_bondnum_B
     * @param neighbor_bondnum_A
     * @param i_bond_neighborsA
     * @param c_bond_neighborsA
     */
    protected TargetProcessor(
            List<String> c_bond_neighborsB,
            List<Integer> iBondSetB,
            List<String> c_bond_setB,
            int neighbor_bondnum_B,
            int set_bondnum_B,
            int neighbor_bondnum_A,
            List<Integer> i_bond_neighborsA,
            List<String> c_bond_neighborsA) {


        this.neighborBondNumB = neighbor_bondnum_B;
        this.setBondNumB = set_bondnum_B;
        this.c_bond_neighborsA = c_bond_neighborsA;
        this.i_bond_neighborsA = i_bond_neighborsA;
        this.neighborBondNumA = neighbor_bondnum_A;
        this.cBondNeighborsB = c_bond_neighborsB;
        this.iBondSetB = iBondSetB;
        this.cBondSetB = c_bond_setB;

    }

    protected void process(IAtomContainer target,
            List<Integer> unmappedAtomMolB,
            int mappingSize,
            List<Integer> iBondNeighborsB,
            List<Integer> mappedAtoms,
            int counter,
            List<String> cTab1,
            List<String> cTab2,
            String[] SignArray) {

        int unmapped_numB = unmappedAtomMolB.size();
        boolean bond_considered = false;
        boolean normal_bond = true;

        for (int a = 0; a < target.getBondCount(); a++) {

            Integer indexI = target.getAtomNumber(target.getBond(a).getAtom(0));
            Integer indexJ = target.getAtomNumber(target.getBond(a).getAtom(1));
            String AtomI = target.getBond(a).getAtom(0).getSymbol();
            String AtomJ = target.getBond(a).getAtom(1).getSymbol();
            Integer order = target.getBond(a).getOrder().ordinal() + 1;

            for (int b = 0; b < unmapped_numB; b++) {
                if (unmappedAtomMolB.get(b).equals(indexI)) {
                    for (int c = 0; c < mappingSize; c++) {
                        if (mappedAtoms.get(c * 2 + 1).equals(indexJ)) {
                            iBondNeighborsB.add(indexI);
                            iBondNeighborsB.add(indexJ);
                            iBondNeighborsB.add(order);
                            if (cTab2.get(a * 4 + 3).compareToIgnoreCase("X") == 0) {
                                cBondNeighborsB.add(cTab2.get(a * 4 + 0));
                                cBondNeighborsB.add(SignArray[counter]);
                                cBondNeighborsB.add("X");
                                cBondNeighborsB.add(cTab2.get(a * 4 + 1));
                                changeCharBonds(indexJ, SignArray[counter], target.getBondCount(), target, cTab2);
                                int cor_atom = searchCorrespondingAtom(mappingSize, indexJ, 2, mappedAtoms);
                                //Commented by Asad
                                changeCharBonds(cor_atom, SignArray[counter], neighborBondNumA, i_bond_neighborsA, c_bond_neighborsA);
//                                changeCharBonds(cor_atom, SignArray[counter], query.getBondCount(), query, cTab1);
                                counter++;
                            } else {
                                cBondNeighborsB.add(cTab2.get(a * 4 + 0));
                                cBondNeighborsB.add(cTab2.get(a * 4 + 1));
                                cBondNeighborsB.add("X");
                                cBondNeighborsB.add(cTab2.get(a * 4 + 3));
                            }
                            normal_bond = false;
                            neighborBondNumB++;
                        }
                    }
                    if (normal_bond) {
                        iBondSetB.add(indexI);
                        iBondSetB.add(indexJ);
                        iBondSetB.add(order);
                        cBondSetB.add(cTab2.get(a * 4 + 0));
                        cBondSetB.add(cTab2.get(a * 4 + 1));
                        cBondSetB.add("X");
                        cBondSetB.add("X");
                        setBondNumB++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                if (unmappedAtomMolB.get(b) == indexJ) {
                    for (int c = 0; c < mappingSize; c++) {
                        if (mappedAtoms.get(c * 2 + 1).equals(indexI)) {
                            iBondNeighborsB.add(indexI);
                            iBondNeighborsB.add(indexJ);
                            iBondNeighborsB.add(order);
                            if (cTab2.get(a * 4 + 2).compareToIgnoreCase("X") == 0) {
                                cBondNeighborsB.add(SignArray[counter]);
                                cBondNeighborsB.add(cTab2.get(a * 4 + 1));
                                cBondNeighborsB.add(cTab2.get(a * 4 + 0));
                                cBondNeighborsB.add("X");
                                changeCharBonds(indexI, SignArray[counter], target.getBondCount(), target, cTab2);
                                int cor_atom = searchCorrespondingAtom(mappingSize, indexI, 2, mappedAtoms);
                                changeCharBonds(cor_atom, SignArray[counter], neighborBondNumA, i_bond_neighborsA, c_bond_neighborsA);
//                                changeCharBonds(cor_atom, SignArray[counter], query.getBondCount(), query, cTab1);
                                counter++;
                            } else {
                                cBondNeighborsB.add(cTab2.get(a * 4 + 0));
                                cBondNeighborsB.add(cTab2.get(a * 4 + 1));
                                cBondNeighborsB.add(cTab2.get(a * 4 + 2));
                                cBondNeighborsB.add("X");
                            }
                            normal_bond = false;
                            neighborBondNumB++;
                        }
                    }
                    if (normal_bond) {

//                        System.out.println("\n\n normal_bond:" + normal_bond);

                        iBondSetB.add(indexI);
                        iBondSetB.add(indexJ);
                        iBondSetB.add(order);
                        cBondSetB.add(AtomI);
                        cBondSetB.add(AtomJ);
                        cBondSetB.add("X");
                        cBondSetB.add("X");
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
     * @param set_num_B
     * @param unmappedAtomMolB
     * @param mappingSize
     * @param iBondNeighborsB
     * @param mappedAtoms
     * @param counter
     * @param cTab1
     * @param cTab2
     * @param SignArray
     */
    protected void process(int set_num_B,
            List<Integer> unmappedAtomMolB,
            int mappingSize,
            List<Integer> iBondNeighborsB,
            List<Integer> mappedAtoms,
            int counter,
            List<String> cTab1,
            List<String> cTab2,
            String[] SignArray) {

        int unmapped_numB = unmappedAtomMolB.size();
        boolean bond_considered = false;
        boolean normal_bond = true;

        for (int a = 0; a < set_num_B; a++) {
            for (int b = 0; b < unmapped_numB; b++) {
                if (unmappedAtomMolB.get(b).equals(iBondSetB.get(a * 3 + 0))) {
                    for (int c = 0; c < mappingSize; c++) {
                        if (mappedAtoms.get(c * 2 + 1).equals(iBondSetB.get(a * 3 + 1))) {
                            iBondNeighborsB.add(iBondSetB.get(a * 3 + 0));
                            iBondNeighborsB.add(iBondSetB.get(a * 3 + 1));
                            iBondNeighborsB.add(iBondSetB.get(a * 3 + 2));
                            if (cTab2.get(a * 4 + 3).compareToIgnoreCase("X") == 0) {
                                cBondNeighborsB.add(cTab2.get(a * 4 + 0));
                                cBondNeighborsB.add(SignArray[counter]);
                                cBondNeighborsB.add("X");
                                cBondNeighborsB.add(cTab2.get(a * 4 + 1));
                                changeCharBonds(iBondSetB.get(a * 3 + 1), SignArray[counter], set_num_B, iBondSetB, cTab2);
                                int cor_atom = searchCorrespondingAtom(mappingSize, iBondSetB.get(a * 3 + 1), 2, mappedAtoms);
                                changeCharBonds(cor_atom, SignArray[counter], neighborBondNumA, i_bond_neighborsA, c_bond_neighborsA);
                                counter++;

                            } else {
                                cBondNeighborsB.add(cTab2.get(a * 4 + 0));
                                cBondNeighborsB.add(cTab2.get(a * 4 + 1));
                                cBondNeighborsB.add("X");
                                cBondNeighborsB.add(cTab2.get(a * 4 + 3));
                            }

                            normal_bond = false;
                            neighborBondNumB++;
                        }
                    }
                    if (normal_bond) {
                        iBondSetB.add(iBondSetB.get(a * 3 + 0));
                        iBondSetB.add(iBondSetB.get(a * 3 + 1));
                        iBondSetB.add(iBondSetB.get(a * 3 + 2));
                        cBondSetB.add(cTab2.get(a * 4 + 0));
                        cBondSetB.add(cTab2.get(a * 4 + 1));
                        cBondSetB.add("X");
                        cBondSetB.add("X");
                        setBondNumB++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                if (unmappedAtomMolB.get(b).equals(iBondSetB.get(a * 3 + 1))) {
                    for (int c = 0; c < mappingSize; c++) {

                        if (mappedAtoms.get(c * 2 + 1).equals(iBondSetB.get(a * 3 + 0))) {

                            iBondNeighborsB.add(iBondSetB.get(a * 3 + 0));
                            iBondNeighborsB.add(iBondSetB.get(a * 3 + 1));
                            iBondNeighborsB.add(iBondSetB.get(a * 3 + 2));

                            if (cTab2.get(a * 4 + 2).compareToIgnoreCase("X") == 0) {
                                cBondNeighborsB.add(SignArray[counter]);
                                cBondNeighborsB.add(cTab2.get(a * 4 + 1));
                                cBondNeighborsB.add(cTab2.get(a * 4 + 0));
                                cBondNeighborsB.add("X");
                                changeCharBonds(iBondSetB.get(a * 3 + 0), SignArray[counter], set_num_B, iBondSetB, cTab2);
                                int cor_atom = searchCorrespondingAtom(mappingSize, iBondSetB.get(a * 3 + 0), 2, mappedAtoms);
                                changeCharBonds(cor_atom, SignArray[counter], neighborBondNumA, i_bond_neighborsA, c_bond_neighborsA);
                                counter++;
                            } else {
                                cBondNeighborsB.add(cTab2.get(a * 4 + 0));
                                cBondNeighborsB.add(cTab2.get(a * 4 + 1));
                                cBondNeighborsB.add(cTab2.get(a * 4 + 2));
                                cBondNeighborsB.add("X");
                            }

                            normal_bond = false;
                            neighborBondNumB++;

                        }


                    }

                    if (normal_bond) {
                        iBondSetB.add(iBondSetB.get(a * 3 + 0));
                        iBondSetB.add(iBondSetB.get(a * 3 + 1));
                        iBondSetB.add(iBondSetB.get(a * 3 + 2));
                        cBondSetB.add(cTab2.get(a * 4 + 0));
                        cBondSetB.add(cTab2.get(a * 4 + 1));
                        cBondSetB.add("X");
                        cBondSetB.add("X");
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

    /**
     *
     * @return
     */
    protected List<String> getCBondNeighborsB() {
        return this.cBondNeighborsB;
    }

    protected List<Integer> getIBondSetB() {
        return iBondSetB;
    }

    protected List<String> getCBondSetB() {
        return cBondSetB;
    }
}
