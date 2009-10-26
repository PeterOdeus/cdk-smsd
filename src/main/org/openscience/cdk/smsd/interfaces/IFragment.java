/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.interfaces;

/**
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 */
public interface IFragment extends IMCSBase {

    boolean getFlag();

    boolean[][] getFlagMatrix();
}
