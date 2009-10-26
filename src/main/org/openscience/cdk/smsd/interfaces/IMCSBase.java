/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.interfaces;

import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.interfaces.IAtom;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @contact e-mail: asad@ebi.ac.uk
 */
public interface IMCSBase {

    Vector<Map<IAtom, IAtom>> getAllAtomMapping();

    Vector<TreeMap<Integer, Integer>> getAllMapping();

    Map<IAtom, IAtom> getFirstAtomMapping();

    TreeMap<Integer, Integer> getFirstMapping();
}
