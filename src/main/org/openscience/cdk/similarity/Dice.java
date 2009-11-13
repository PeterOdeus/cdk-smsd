/*  $RCSfile$
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright (C) 2009 Rajarshi Guha <rajarshi.guha@gmail.com>
 *
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *  All we ask is that proper credit is given for our work, which includes
 *  - but is not limited to - adding the above copyright notice to the beginning
 *  of your source code files, and to any copyright notice that you may distribute
 *  with programs based on this work.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */
package org.openscience.cdk.similarity;


import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;

import java.util.BitSet;

/**
 * Calculates the Dice coefficient for a given pair of two fingerprint bitsets or real valued feature vectors.
 * <p/>
 * The Dice coefficient is one way to quantitatively measure the "distance" or similarity of two chemical
 * structures.
 * <p/>
 * <p>You can use the FingerPrinter class to retrieve two fingerprint bitsets. We assume that you have two structures
 * stored in cdk.Molecule objects. A tanimoto coefficient can then be calculated like:
 * <pre>
 *   BitSet fingerprint1 = Fingerprinter.getFingerprint(molecule1);
 *   BitSet fingerprint2 = Fingerprinter.getFingerprint(molecule2);
 *   float dice_coefficient = Dice.calculate(fingerprint1, fingerprint2);
 *  </pre>
 * <p/>
 * <p>The FingerPrinter assumes that hydrogens are explicitely given, if this is desired!
 *
 * @author Rajarshi Guha
 * @cdk.githash
 * @cdk.created 2009-11-13
 * @cdk.keyword dice
 * @cdk.keyword similarity
 * @cdk.keyword sorensen
 * @cdk.module fingerprint
 */
@TestClass("org.openscience.cdk.similarity.CosineTest")
public class Dice {

    /**
     * Evaluates Dice coefficient for two bit sets.
     *
     * @param bitset1 A bitset (such as a fingerprint) for the first molecule
     * @param bitset2 A bitset (such as a fingerprint) for the second molecule
     * @return The Dice coefficient
     * @throws org.openscience.cdk.exception.CDKException
     *          if bitsets are not of the same length
     */
    @TestMethod("testDice1,testDice2")
    public static float calculate(BitSet bitset1, BitSet bitset2) throws CDKException {
        float _bitset1_cardinality = bitset1.cardinality();
        float _bitset2_cardinality = bitset2.cardinality();
        if (bitset1.size() != bitset2.size()) {
            throw new CDKException("Bisets must have the same bit length");
        }

        BitSet one_and_two = (BitSet) bitset1.clone();
        one_and_two.and(bitset2);
        float _common_bit_count = one_and_two.cardinality();
        return 2 * _common_bit_count / (_bitset1_cardinality + _bitset2_cardinality);                
    }
}