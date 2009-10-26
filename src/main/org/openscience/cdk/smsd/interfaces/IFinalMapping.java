/*
 * IFinalMapping.java
 *
 * Created on Sep 12, 2007, 8:58:45 PM
 *
 */
package org.openscience.cdk.smsd.interfaces;

import java.util.Iterator;
import java.util.TreeMap;
import java.util.Vector;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public interface IFinalMapping {

    /**
     * 
     * @param v MCS mapping between a given
     * reactant and product 
     * @throws java.lang.Exception
     */
    public void addElement(TreeMap<Integer, Integer> v) throws Exception;

    /**
     * 
     * @param v List of all MCS mapping between a given
     * reactant and product 
     * @throws java.lang.Exception
     */
    public void add(TreeMap<Integer, Integer> v) throws Exception;

    /**
     * 
     * @param v List of all MCS mapping between a given
     * reactant and product 
     */
    public void set(Vector<TreeMap<Integer, Integer>> v);

    /**
     * 
     * @return Iterator of mappings
     */
    public Iterator getIterator();

    /**
     * clear the maping
     */
    public void Clear();

    /**
     * 
     * @return get of MCS mapping List
     */
    public Vector<TreeMap<Integer, Integer>> getFinalMapping();

    /**
     * 
     * @return size of the mapping
     */
    public int getSize();
}
