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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smsd.helper.MolHandler;

/**
 * @cdk.module smsd
 */
public interface IMCSAlgorithm extends IMCSBase {

    /**
     * 0: default, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS
     */
    public enum Algorithm {

        DEFAULT(0, "Default SMSD algorithm"),
        MCSPlus(1, "MCS Plus algorithm"),
        VFLibMCS(2, "VF Lib based MCS algorithm"),
        CDKMCS(3, "CDK UIT MCS"),
        SubStructure(4, "Turbo Mode based Substructure search");
        private final int type;
        private final String description;

        Algorithm(int aStatus, String desc) {
            this.type = aStatus;
            this.description = desc;
        }

        public int type() {
            return this.type;
        }

        /**
         *
         * @return
         */
        public String description() {
            return this.description;
        }

        /**
         * 
         * @param <status>
         * @param obj
         * @return
         */
        public <status> int compareTo(Algorithm obj) {
            return 0;
        }
    }

    /**
     *
     * @param source
     * @param target
     * @param removeHydrogen true if remove H before mapping
     * @throws CDKException
     */
    void init(IMolecule source, IMolecule target, boolean removeHydrogen) throws CDKException;

    /**
     *
     * @param source
     * @param target
     * @param removeHydrogen true if remove H before mapping
     * @throws CDKException
     */
    void init(IAtomContainer source, IAtomContainer target, boolean removeHydrogen) throws CDKException;

    /**
     *
     * @param stereoFilter
     * @param fragmentFilter
     * @param energyFilter
     */
    void setChemFilters(boolean stereoFilter, boolean fragmentFilter, boolean energyFilter);

    /**
     * 
     * @param Key Index of the mapping solution
     * @return Total bond breaking energy required to remove the mapped part
     */
    Double getEnergyScore(int Key);

    /**
     *
     * @param Key Index of the mapping solution
     * @return Fragment count(s) generated after removing the mapped parts
     */
    Integer getFragmentSize(int Key);

    /**
     *
     * @return return modified Product Molecule
     */
    IAtomContainer getProductMolecule();

    /**
     *
     * @return return modified Reactant Molecule
     */
    IAtomContainer getReactantMolecule();

    /**
     *
     * @param Key Index of the mapping solution
     * @return true if no stereo mismatch occures
     * else false if stereo mismatch occures
     */
    Integer getStereoScore(int Key);

    /**
     *
     * @return true if two molecules have same stereo match
     */
    boolean isStereoMisMatch();

    /**
     *
     * @return true if query molecule is a subgraph of the target molecule
     */
    boolean isSubgraph();

    /**
     *
     * @return Tanimoto Similarity between 0 and 1
     * @throws IOException
     */
    double getTanimotoSimilarity() throws IOException;

    /**
     *
     * @return Euclidean Distance (lower the score, better the match)
     * @throws IOException
     */
    double getEuclideanDistance() throws IOException;
}
