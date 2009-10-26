/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.core.tools;

import java.util.BitSet;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 * 
 * 
 * <PRE>
 *
 *   The Measures program takes as input any fixed length bit strings,
 *   these can be from the Mesa Fingerprint programs or user supplied fingerprints.
 *   User supplied fingerprints must take the form of ASCII 1's and 0's, (e.g. 011100001111000....),
 *   ASCII CDK fingerprints inside the FP<> data type are also valid input to  Measures .
 *   The Measures program  produces a similarity or dissimilarity matrix (user's choice)
 *   using one of the following user selected measures:
 *   <B> Tversky, Tanimoto, Euclidean, Hamman, or Ochia (1-Cosine).
 *
 * In similarity form:
 *
 *                            Tanimoto(A,B)  = c / [a + b - c]  (symmetric)
 *
 *                            Euclidean(A,B) = 1 - {[(a + b)] / n}(1/2)   (symmetric)
 *
 *                            Hamman(A,B)  = [c + d] /n  (symmetric)
 *
 *                            Ochia(A,B) = 1 - Cosine(A,B) = c / [(c + a) * (c + b)](1/2)  (symmetric)
 *
 *                            Tversky(A,B) = c / [(alpha) * a + (beta) * b + c]  (asymmetric)
 *
 *                            a : Unique bits turned on in molecule "A"
 *                            b:  Unique bits turned on in molecule "B"
 *                            c:  Common bits turned on in both molecule "A" and molecule "B"
 *                            d:  Common bits turned off in both molecule "A" and molecule "B"
 *                            n:  The total number of bits in the fingerprint
 * </B>
 *
 *  <B> Note:</B>The Tanimoto, Euclidean, Hamman, and Ochai are all symmetric measures.
 *
 *  <U> This means that the comparison of A to B yields the same number as the comparison of compound B to compound A.</U>
 *  <B> Note:</B> The dissimilarity is just 1 - similarity.
 *
 *
 *
 * </PRE>
 *
 * @ref <B>Willett et.al., Chemical Similarity Searching,</B> <I>J.Chem. Inf. Comput. Sci.</I>, Vol. 38, No. 6, 1998
 *
 *
 */
public class Similarity {

    private static BitSet A;
    private static BitSet B;

    /**
     * 
     * @param Molecule1 BitSet
     * @param Molecule2 BitSet
     * @return <B>Similarity <U>Tanimoto, Jaccard</U> </B>
     * <B>c/(a+b-c)></B>
     * @throws java.lang.Exception
     */
    public static float getTanimotoSimilarity(BitSet Molecule1, BitSet Molecule2) throws Exception {

        A = (BitSet) Molecule1.clone();
        B = (BitSet) Molecule2.clone();

        float _bitset1_cardinality = A.cardinality();
        float _bitset2_cardinality = B.cardinality();

//        System.out.println("A: "+ A.size() + " " + " B" + B.size());
        if (A.size() != B.size()) {
            throw new Exception("BitSets must have the same bit length");
        }
        BitSet one_and_two = (BitSet) A.clone();
        one_and_two.and(B);
        float _common_bit_count = one_and_two.cardinality();
        return _common_bit_count / (_bitset1_cardinality + _bitset2_cardinality - _common_bit_count);

    }

    /**
     *
     * @param Molecule1 
     * @param Molecule2 
     * @return <B>Similarity <U>Cosine,Ochiai,Carbo</U></B>
     * <B>c/sqrt(a*b)</B>
     * @throws Exception 
     */
    public static double getCosineSimilarity(BitSet Molecule1, BitSet Molecule2) throws Exception {

        A = (BitSet) Molecule1.clone();
        B = (BitSet) Molecule2.clone();

        float _bitset1_cardinality = A.cardinality();
        float _bitset2_cardinality = B.cardinality();

        if (A.size() != B.size()) {
            throw new Exception("Bisets must have the same bit length");
        }
        BitSet one_and_two = (BitSet) A.clone();
        one_and_two.and(B);
        float _common_bit_count = one_and_two.cardinality();

        return _common_bit_count / (Math.sqrt(_bitset1_cardinality * _bitset2_cardinality));

    }

    /**
     *
     * @param Molecule1 
     * @param Molecule2 
     * @return <B>Similarity <U>Dice, Sorensen, Czekanowski, Hodgkin-Richards</U></B>
     * <B>2c/(a+b)</B>
     * @throws Exception 
     *
     */
    public static double getDiceSimilarity(BitSet Molecule1, BitSet Molecule2) throws Exception {

        A = (BitSet) Molecule1.clone();
        B = (BitSet) Molecule2.clone();

        float _bitset1_cardinality = A.cardinality();
        float _bitset2_cardinality = B.cardinality();

        if (A.size() != B.size()) {
            throw new Exception("Bisets must have the same bit length");
        }
        BitSet one_and_two = (BitSet) A.clone();
        one_and_two.and(B);
        float _common_bit_count = one_and_two.cardinality();

        return 2 * _common_bit_count / (_bitset1_cardinality + _bitset2_cardinality);

    }
}
