/*
 * FinalMapping.java
 *
 * Created on 09 February 2007, 11:38
 *
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
package org.openscience.cdk.smsd.helper;

import org.openscience.cdk.smsd.interfaces.IFinalMapping;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import java.util.Vector;

public class FinalMappings implements IFinalMapping {

    private List<TreeMap<Integer, Integer>> _mappings;
    private static FinalMappings INSTANCE = null;

    protected FinalMappings() {
        _mappings = new Vector<TreeMap<Integer, Integer>>();
    }

    synchronized public static FinalMappings getInstance() {
        if (INSTANCE == null) {

            INSTANCE = new FinalMappings();

        }
        return INSTANCE;
    }

    @Override
    synchronized public void add(TreeMap<Integer, Integer> v) {
        //System.err.println("Vector added");
        _mappings.add(v);
    }

    /**
     * 
     * @param v
     */
    @Override
    synchronized public final void set(List<TreeMap<Integer, Integer>> v) {

        _mappings.clear();
        for (TreeMap<Integer, Integer> M : v) {
            _mappings.add(M);
        }
    }

    @Override
    synchronized public Iterator<TreeMap<Integer, Integer>> getIterator() {
        Iterator<TreeMap<Integer, Integer>> it = _mappings.iterator();
        return it;
    }

    @Override
    synchronized public void Clear() {
        _mappings.clear();
    }

    @Override
    synchronized public List<TreeMap<Integer, Integer>> getFinalMapping() {
        return _mappings;
    }

    @Override
    synchronized public int getSize() {
        return _mappings.size();
    }
}
