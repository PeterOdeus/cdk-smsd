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

package org.openscience.cdk.smsd.interfaces;


import java.io.IOException;
import org.openscience.cdk.smsd.core.tools.EBIException;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 */
public interface IMCS extends IMCSBase {

  

    /**
     *
     * @param removeHydrogen 
     * @return
     * @throws java.io.IOException
     * @throws EBIException
     */
    int search_MCS(boolean removeHydrogen) throws IOException, EBIException;

    /**
     * Creates a new instance of SearchCliques
     * @param Reactant
     * @param Product
     * @throws java.io.IOException
     *
     *
     */
    void set(MolHandler Reactant, MolHandler Product) throws IOException;

    /**
     * Creates a new instance of SearchCliques
     * @param ReactantMolFileName
     * @param ProductMolFileName
     * @throws java.io.IOException
     */
    void set(String ReactantMolFileName, String ProductMolFileName) throws IOException;

    /**
     * Creates a new instance of SearchCliques
     * @param ReactantMol
     * @param ProductMol
     * @throws java.io.IOException
     */
    void set(IAtomContainer ReactantMol, IAtomContainer ProductMol) throws IOException;


}
