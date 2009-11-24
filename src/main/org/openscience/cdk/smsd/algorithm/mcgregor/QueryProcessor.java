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
    private String[] SignROW;
    private int neighborBondNumA = 0; //number of remaining molecule A bonds after the clique search, which are neighbors of the MCS_1
    private int setBondNumA = 0; //number of remaining molecule A bonds after the clique search, which aren't neighbors
    private IAtomContainer query;
    private IAtomContainer target;
    private List<Integer> i_bond_neighborsA;
    private List<String> c_bond_neighborsA;

    /**
     * 
     * @param query
     * @param target
     * @param c_tab1_copy
     * @param c_tab2_copy
     * @param SignROW
     * @param neighbor_bondnum_A
     * @param set_bondnum_A
     */
    protected QueryProcessor(IAtomContainer query, IAtomContainer target, List<String> c_tab1_copy, List<String> c_tab2_copy, String[] SignROW, int neighbor_bondnum_A, int set_bondnum_A) {

        this.query = query;
        this.target = target;

        this.c_tab1_copy = c_tab1_copy;
        this.c_tab2_copy = c_tab2_copy;
        this.SignROW = SignROW;
        this.neighborBondNumA = neighbor_bondnum_A;
        this.setBondNumA = set_bondnum_A;
    }

    /**
     * 
     * @param unmapped_atoms_molA
     * @param mappingSize
     * @param i_bond_neighborsA
     * @param i_bond_setA
     * @param c_bond_neighborsA
     * @param c_bond_setA
     * @param mapped_atoms
     * @param SR_count
     */
    protected void process(
            List<Integer> unmapped_atoms_molA,
            int mappingSize,
            List<Integer> i_bond_neighborsA,
            List<Integer> i_bond_setA,
            List<String> c_bond_neighborsA,
            List<String> c_bond_setA,
            List<Integer> mapped_atoms,
            int SR_count) {

        this.i_bond_neighborsA = i_bond_neighborsA;
        this.c_bond_neighborsA = c_bond_neighborsA;

//        int mappingSize = clique_vector.size();

        boolean bond_considered = false;
        boolean normal_bond = true;

//        System.out.println("\n" + c_tab1_copy + "\n");


        for (int a = 0; a < query.getBondCount(); a++) {


            Integer indexI = query.getAtomNumber(query.getBond(a).getAtom(0));
            Integer indexJ = query.getAtomNumber(query.getBond(a).getAtom(1));
            String AtomI = query.getBond(a).getAtom(0).getSymbol();
            String AtomJ = query.getBond(a).getAtom(1).getSymbol();
            Integer order = query.getBond(a).getOrder().ordinal() + 1;

//            System.out.println(AtomI + "= , =" + AtomJ );
            for (Integer unMappedAtomIndex : unmapped_atoms_molA) {



                if (unMappedAtomIndex.equals(indexI)) {
//                    System.out.println("Unmapped Atom is equal to Reaction Bond table");
//                    System.out.println("unMappedAtomIndex " + unMappedAtomIndex + " " + "indexI: " + indexI);
//                    System.out.println("unMappedAtomIndex=IndexI " + query.getAtom(unMappedAtomIndex).getSymbol());
                    for (int c = 0; c < mappingSize; c++) {

                        if (mapped_atoms.get(c * 2).equals(indexJ)) {

                            i_bond_neighborsA.add(indexI);
                            i_bond_neighborsA.add(indexJ);
                            i_bond_neighborsA.add(order);

                            if (c_tab1_copy.get(a * 4 + 3).compareToIgnoreCase("X") == 0) {


                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add(SignROW[SR_count]);
                                c_bond_neighborsA.add("X");
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));

                                changeCharBonds(indexJ, SignROW[SR_count], query.getBondCount(), query, c_tab1_copy);

                                int cor_atom = searchCorrespondingAtom(mappingSize, indexJ, 1, mapped_atoms);
                                changeCharBonds(cor_atom, SignROW[SR_count], target.getBondCount(), target, c_tab2_copy);
                                SR_count++;
                            } else {
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 2));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 3));
                            }
                            normal_bond = false;
                            neighborBondNumA++;
                        }
                    }
                    if (normal_bond) {
                        i_bond_setA.add(indexI);
                        i_bond_setA.add(indexJ);
                        i_bond_setA.add(order);
                        c_bond_setA.add(c_tab1_copy.get(a * 4 + 0));
                        c_bond_setA.add(c_tab1_copy.get(a * 4 + 1));
                        c_bond_setA.add("X");
                        c_bond_setA.add("X");
                        setBondNumA++;
                    }
                    normal_bond = true;
                    bond_considered = true;
                }
                //Does a ungemaptes atom at second position in the connection occur?
                if (unMappedAtomIndex.equals(indexJ)) {
                    for (int c = 0; c < mappingSize; c++) {


                        if (mapped_atoms.get(c * 2 + 0).equals(indexI)) {



                            i_bond_neighborsA.add(indexI);
                            i_bond_neighborsA.add(indexJ);
                            i_bond_neighborsA.add(order);

                            if (c_tab1_copy.get(a * 4 + 2).compareToIgnoreCase("X") == 0) {


                                c_bond_neighborsA.add(SignROW[SR_count]);
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add("X");
                                changeCharBonds(indexI, SignROW[SR_count],
                                        query.getBondCount(), query, c_tab1_copy);

                                int cor_atom = searchCorrespondingAtom(mappingSize, indexI, 1, mapped_atoms);
                                changeCharBonds(cor_atom, SignROW[SR_count], target.getBondCount(), target, c_tab2_copy);
                                SR_count++;
                            } else {
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 0));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 1));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 2));
                                c_bond_neighborsA.add(c_tab1_copy.get(a * 4 + 3));
                            }
                            normal_bond = false;
                            neighborBondNumA++;
                            //System.out.println("Neighbor");
                            //System.out.println(neighborBondNumA);
                        }
                    }
                    if (normal_bond) {
                        i_bond_setA.add(indexI);
                        i_bond_setA.add(indexJ);
                        i_bond_setA.add(order);
                        c_bond_setA.add(AtomI);
                        c_bond_setA.add(AtomJ);
                        c_bond_setA.add("X");
                        c_bond_setA.add("X");
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

    }

//    protected void process(List<Integer> unmapped_atoms_molA,
//            int mappingSize,
//            List<Integer> i_bond_neighborsA,
//            List<Integer> i_bond_setA,
//            List<String> c_bond_neighborsA,
//            List<String> c_bond_setA,
//            List<Integer> mapped_atoms,
//            int SR_count){
//
//          boolean bond_considered = false;
//            boolean normal_bond = true;
//            for (int a = 0; a < set_num_A; a++) {
//
//                int _elementAt_a = i_bond_setA.get(a * 3 + 0).intValue();
//                for (int b = 0; b < unmapped_numA; b++) {
//                    Integer unMappedAtomIndex = unmapped_atoms_molA.get(b);
//                    if (unMappedAtomIndex == _elementAt_a) {
//                        for (int c = 0; c < new_MAPPING_size; c++) {
//
//                            if (new_MAPPING.get(c * 2 + 0).equals(i_bond_setA.get(a * 3 + 1))) {
//
//                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 0));
//                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 1));
//                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 2));
//                                new_c_neighborsA.add(cTab1Copy.get(a * 4 + 0));
//                                if (cTab1Copy.get(a * 4 + 3).compareToIgnoreCase("X") == 0) {
//
//                                    new_c_neighborsA.add(SignArray[counter]);
//                                    new_c_neighborsA.add("X");
//                                    new_c_neighborsA.add(cTab1Copy.get(a * 4 + 1));
//                                    changeCharBonds(i_bond_setA.get(a * 3 + 1), SignArray[counter], set_num_A, i_bond_setA, cTab1Copy);
//                                    int cor_atom = McGregorChecks.searchCorrespondingAtom(new_MAPPING_size, i_bond_setA.get(a * 3 + 1), 1, new_MAPPING);
//                                    changeCharBonds(cor_atom, SignArray[counter], set_num_B, i_bond_setB, cTab2Copy);
//                                    counter++;
//
//                                } else {
//
//                                    new_c_neighborsA.add(cTab1Copy.get(a * 4 + 1));
//                                    new_c_neighborsA.add("X");
//                                    new_c_neighborsA.add(cTab1Copy.get(a * 4 + 3));
//
//                                }
//
//                                normal_bond = false;
//                                new_neighbor_numA++;
//
//                            }
//                        }
//
//                        if (normal_bond) {
//
//                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 0));
//                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 1));
//                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 2));
//                            new_c_bond_setA.add(cTab1Copy.get(a * 4 + 0));
//                            new_c_bond_setA.add(cTab1Copy.get(a * 4 + 1));
//                            new_c_bond_setA.add("X");
//                            new_c_bond_setA.add("X");
//                            setBondNumA++;
//
//                        }
//                        normal_bond = true;
//                        bond_considered = true;
//                    }
//                    if (unMappedAtomIndex.equals(i_bond_setA.get(a * 3 + 1))) {
//                        for (int c = 0; c < new_MAPPING_size; c++) {
//
//                            if (new_MAPPING.get(c * 2 + 0).equals(i_bond_setA.get(a * 3 + 0))) {
//
//                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 0));
//                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 1));
//                                new_i_neighborsA.add(i_bond_setA.get(a * 3 + 2));
//                                if (cTab1Copy.get(a * 4 + 2).compareToIgnoreCase("X") == 0) {
//
//                                    new_c_neighborsA.add(SignArray[counter]);
//                                    new_c_neighborsA.add(cTab1Copy.get(a * 4 + 1));
//                                    new_c_neighborsA.add(cTab1Copy.get(a * 4 + 0));
//                                    new_c_neighborsA.add("X");
//                                    changeCharBonds(i_bond_setA.get(a * 3 + 0), SignArray[counter], set_num_A, i_bond_setA, cTab1Copy);
//                                    int cor_atom = McGregorChecks.searchCorrespondingAtom(new_MAPPING_size, i_bond_setA.get(a * 3 + 0), 1, new_MAPPING);
//                                    changeCharBonds(cor_atom, SignArray[counter], set_num_B, i_bond_setB, cTab2Copy);
//                                    counter++;
//
//                                } else {
//                                    new_c_neighborsA.add(cTab1Copy.get(a * 4 + 0));
//                                    new_c_neighborsA.add(cTab1Copy.get(a * 4 + 1));
//                                    new_c_neighborsA.add(cTab1Copy.get(a * 4 + 2));
//                                    new_c_neighborsA.add("X");
//                                }
//
//                                normal_bond = false;
//                                new_neighbor_numA++;
//
//                            }
//
//
//                        }
//                        if (normal_bond) {
//                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 0));
//                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 1));
//                            new_i_bond_setA.add(i_bond_setA.get(a * 3 + 2));
//                            new_c_bond_setA.add(cTab1Copy.get(a * 4 + 0));
//                            new_c_bond_setA.add(cTab1Copy.get(a * 4 + 1));
//                            new_c_bond_setA.add("X");
//                            new_c_bond_setA.add("X");
//                            setBondNumA++;
//                        }
//
//                        normal_bond = true;
//                        bond_considered = true;
//                    }
//
//                    if (bond_considered) {
//                        break;
//                    }
//                }
//                bond_considered = false;
//            }
//    }

    private int searchCorrespondingAtom(int mapped_atoms_size, int atom_from_other_molecule, int molecule, List<Integer> mapped_atoms_org) {


        List<Integer> mapped_atoms = new ArrayList<Integer>(mapped_atoms_org);

        int corresponding_atom = 0;
        for (int a = 0; a < mapped_atoms_size; a++) {
            if ((molecule == 1) &&
                    (mapped_atoms.get(a * 2 + 0).intValue() == atom_from_other_molecule)) {

                corresponding_atom = mapped_atoms.get(a * 2 + 1);


            }
            if ((molecule == 2) &&
                    (mapped_atoms.get(a * 2 + 1).intValue() == atom_from_other_molecule)) {
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
        return this.i_bond_neighborsA;
    }

    List<String> getCBondNeighborsA() {
        return this.c_bond_neighborsA;
    }
}
