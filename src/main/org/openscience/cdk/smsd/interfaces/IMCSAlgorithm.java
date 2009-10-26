/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.interfaces;

import java.io.IOException;
import org.openscience.cdk.smsd.core.tools.EBIException;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @contact e-mail: asad@ebi.ac.uk
 */
public interface IMCSAlgorithm extends IMCSBase {

    Double getEnergyScore(int Key);

    Integer getFragmentSize(int Key);

    IAtomContainer getProductMolecule();

    IAtomContainer getReactantMolecule();

    Integer getStereoScore(int Key);

    boolean isStereoMisMatch();

    boolean isSubgraph();

    double getTanimotoSimilarity() throws IOException;

    double getEuclideanDistance() throws IOException;

//    /**
//     *
//     * @param Reactant
//     * @param Product
//     * @param ReactantFingerprint
//     * @param ProductFingerprint
//     * @throws EBIException
//     */
//    void init(IMolecule Reactant, IMolecule Product, BitSet ReactantFingerprint, BitSet ProductFingerprint)throws EBIException;
//
//    /**
//     *
//     * @param Reactant
//     * @param Product
//     * @param ReactantFingerprint
//     * @param ProductFingerprint
//     * @throws EBIException
//     */
//    void init(IAtomContainer Reactant, IAtomContainer Product, BitSet ReactantFingerprint, BitSet ProductFingerprint)throws EBIException;
//
//    /**
//     *
//     * @param Reactant
//     * @param Product
//     * @param ReactantFingerprint
//     * @param ProductFingerprint
//     * @throws EBIException
//     */
//    void init(MolHandler Reactant, MolHandler Product, BitSet ReactantFingerprint, BitSet ProductFingerprint)throws EBIException;
    /**
     *
     * @param Reactant
     * @param Product
     * @throws EBIException 
     *
     */
    void init(MolHandler Reactant, MolHandler Product) throws EBIException;

    /**
     *
     * @param Reactant
     * @param Product
     * @throws EBIException 
     */
    void init(IMolecule Reactant, IMolecule Product) throws EBIException;

    /**
     *
     * @param Reactant
     * @param Product
     * @throws EBIException 
     */
    void init(IAtomContainer Reactant, IAtomContainer Product) throws EBIException;
}
