/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.global;
/**
 *
 * Copyright (C) 2006-2009  Syed Asad Rahman {asad@ebi.ac.uk}
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/**
 * @cdk.module smsd
 */
public class BondType {

    private static BondType _instance = null;
    private static boolean bondSensitive = false;

    /**
     *
     * @return
     */
    public static synchronized BondType getInstance() {
        if (_instance == null) {
            // it's ok, we can call this constructor
            _instance = new BondType();
        }
        return _instance;
    }

    protected BondType() {
    }

    /**
     *
     * @param isBondSensitive (set true if bondsensetive else false)
     */
    public void setBondSensitiveFlag(boolean isBondSensitive) {
        BondType.bondSensitive = isBondSensitive;
    }

    /**
     *
     * @return true if bondsensitive else false
     */
    public boolean getBondSensitiveFlag() {
        return bondSensitive;
    }

    /**
     * reset bond sensitive flag to default
     */
    public void reset() {
        bondSensitive = false;
    }
}
