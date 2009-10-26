/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.openscience.cdk.smsd.interfaces;


import java.io.IOException;
import org.openscience.cdk.smsd.core.tools.EBIException;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 */
public interface IMCS extends IMCSBase {

  

    /**
     *
     * @param removeHydrogen 
     * @return
     * @throws java.io.IOException
     * @throws EBIException
     */
    int search_MCS(boolean removeHydrogen) throws IOException, EBIException;

    /**
     * Creates a new instance of SearchCliques
     * @param Reactant
     * @param Product
     * @throws java.io.IOException
     *
     *
     */
    void set(MolHandler Reactant, MolHandler Product) throws IOException;

    /**
     * Creates a new instance of SearchCliques
     * @param ReactantMolFileName
     * @param ProductMolFileName
     * @throws java.io.IOException
     */
    void set(String ReactantMolFileName, String ProductMolFileName) throws IOException;

    /**
     * Creates a new instance of SearchCliques
     * @param ReactantMol
     * @param ProductMol
     * @throws java.io.IOException
     */
    void set(IAtomContainer ReactantMol, IAtomContainer ProductMol) throws IOException;


}
