/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.global;

/**
 *
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
public class BondType {

    private static BondType INSTANCE = null;
    private static boolean bondtype = false;

    /**
     *
     * @return
     */
    public static synchronized BondType getInstance() {
        if (INSTANCE == null) {

            // it's ok, we can call this constructor
            INSTANCE = new BondType();
        }

        return INSTANCE;
    }

    protected BondType() {
    }

    /**
     *
     * @param isBondSensitive (set true if bondsensetive else false)
     */
    public void setBondSensitiveFlag(boolean isBondSensitive) {

        BondType.bondtype = isBondSensitive;
    }

    /**
     *
     * @return true if bondsensitive else false
     */
    public boolean getBondSensitiveFlag() {

        return bondtype;
    }

    public void reset() {
        bondtype = false;
    }
}
