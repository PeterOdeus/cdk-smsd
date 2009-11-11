/* Copyright (C) 2006-2009  Syed Asad Rahman {asad@ebi.ac.uk}
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
package org.openscience.cdk.smsd.interfaces;

import java.io.IOException;
import org.openscience.cdk.smsd.core.tools.EBIException;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @cdk.module smsd
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
