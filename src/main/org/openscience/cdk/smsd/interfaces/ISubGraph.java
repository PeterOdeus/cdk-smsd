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
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public interface ISubGraph extends IMCSBase {

    /**
     *
     * @return true if Query/Reactant is a subgraph of Target/Product
     * else false
     * @throws java.io.IOException
     * @throws EBIException 
     */
    boolean isSubgraph() throws IOException, EBIException;

    /**
     * Creates a new instance of SearchCliques
     * @param Query
     * @param Target
     * @param removeHydrogen 
     * @throws java.io.IOException
     *
     *
     */
    void set(MolHandler Query, MolHandler Target, boolean removeHydrogen) throws IOException;

    /**
     * Creates a new instance of SearchCliques
     * @param ReactantMolFileName
     * @param ProductMolFileName
     * @param removeHydrogen 
     * @throws java.io.IOException
     */
    void set(String ReactantMolFileName, String ProductMolFileName, boolean removeHydrogen) throws IOException;

    /**
     * Creates a new instance of SearchCliques
     * @param ReactantMol
     * @param ProductMol
     * @param removeHydrogen 
     * @throws java.io.IOException
     */
    void set(IAtomContainer ReactantMol, IAtomContainer ProductMol, boolean removeHydrogen) throws IOException;
//    void setAllAtomMapping();
//
//    void setAllMapping();
//
//    void setFirstAtomMapping();
//
//    void setFirstMapping();
}
