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
import org.openscience.cdk.interfaces.IMolecule;

/**
 * @cdk.module smsd
 */
public interface IMCSAlgorithm extends IMCSBase {

    Double getEnergyScore(int Key);

    Integer getFragmentSize(int Key);

    IAtomContainer getProductMolecule();

    IAtomContainer getReactantMolecule();

    Integer getStereoScore(int Key);

    boolean isStereoMisMatch();

    boolean isSubgraph();

    double getTanimotoSimilarity() throws IOException;

    double getEuclideanDistance() throws IOException;

//    /**
//     *
//     * @param Reactant
//     * @param Product
//     * @param ReactantFingerprint
//     * @param ProductFingerprint
//     * @throws EBIException
//     */
//    void init(IMolecule Reactant, IMolecule Product, BitSet ReactantFingerprint, BitSet ProductFingerprint)throws EBIException;
//
//    /**
//     *
//     * @param Reactant
//     * @param Product
//     * @param ReactantFingerprint
//     * @param ProductFingerprint
//     * @throws EBIException
//     */
//    void init(IAtomContainer Reactant, IAtomContainer Product, BitSet ReactantFingerprint, BitSet ProductFingerprint)throws EBIException;
//
//    /**
//     *
//     * @param Reactant
//     * @param Product
//     * @param ReactantFingerprint
//     * @param ProductFingerprint
//     * @throws EBIException
//     */
//    void init(MolHandler Reactant, MolHandler Product, BitSet ReactantFingerprint, BitSet ProductFingerprint)throws EBIException;
    /**
     *
     * @param Reactant
     * @param Product
     * @throws EBIException 
     *
     */
    void init(MolHandler Reactant, MolHandler Product) throws EBIException;

    /**
     *
     * @param Reactant
     * @param Product
     * @throws EBIException 
     */
    void init(IMolecule Reactant, IMolecule Product) throws EBIException;

    /**
     *
     * @param Reactant
     * @param Product
     * @throws EBIException 
     */
    void init(IAtomContainer Reactant, IAtomContainer Product) throws EBIException;
}
