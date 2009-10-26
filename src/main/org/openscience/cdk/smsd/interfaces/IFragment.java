/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.openscience.cdk.smsd.interfaces;

/**
 *
 * @author sar
 */
public interface IFragment extends IMCSBase{

//    Vector<List<IAtom>> getAllAtomMapping();
//
//    Vector<List<Integer>> getAllMapping();
//
//    List<IAtom> getFirstAtomMapping();
//
//    List<Integer> getFirstMapping();

    boolean getFlag();

    boolean[][] getFlagMatrix();

//    double getSimilarity() throws IOException;

}
