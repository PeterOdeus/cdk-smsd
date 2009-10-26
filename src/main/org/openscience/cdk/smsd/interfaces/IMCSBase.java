/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.interfaces;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.interfaces.IAtom;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @contact e-mail: asad@ebi.ac.uk
 */
public interface IMCSBase {

    /**
     * 
     * @return All possible MCS atom Mappings
     */
    List<Map<IAtom, IAtom>> getAllAtomMapping();

    /**
     * 
     * @return All possible MCS Mapping Index
     */
    List<TreeMap<Integer, Integer>> getAllMapping();

    /**
     *
     * @return Best Atom Mapping
     */
    Map<IAtom, IAtom> getFirstAtomMapping();

    /**
     *
     * @return Best Mapping Index
     */
    TreeMap<Integer, Integer> getFirstMapping();
}
