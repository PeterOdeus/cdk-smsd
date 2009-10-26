/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.openscience.cdk.smsd.interfaces;

import java.io.IOException;
import java.util.BitSet;
import java.util.Map;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public interface IInChIContainer {

    /**
     *
     * @throws java.io.IOException
     */
    void Clear() throws IOException;

    /**
     *
     * @param Key
     * @throws java.io.IOException
     */
    void Erase(String Key) throws IOException;

    Object clone() throws CloneNotSupportedException;

    /**
     *
     * @param Key
     * @throws java.io.IOException
     * @return
     */
    String getInChI(String Key) throws IOException;

    /**
     *
     * @throws java.io.IOException
     * @return
     */
    Map<String, String> getInChIMap() throws IOException;

    /**
     *
     * @param Value
     * @return
     * @throws java.io.IOException
     */
    String getMoleculeID(BitSet Value) throws IOException;

    /**
     *
     * @param Key
     * @throws java.io.IOException
     * @return
     */
    boolean isKeyPresent(String Key) throws IOException;

    /**
     *
     * @param Value
     * @throws java.io.IOException
     * @return
     */
    boolean isValuePresent(String Value) throws IOException;

    /**
     *
     * @param Key
     * @param Value
     * @throws java.io.IOException
     */
    void put(String Key, String Value) throws IOException;

    /**
     *
     * @param Key
     * @param Value
     * @throws java.io.IOException
     */
    void setValue(String Key, String Value) throws IOException;

    /**
     *
     * @throws java.io.IOException
     */
    void write() throws IOException;

}
