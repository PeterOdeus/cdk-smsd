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
 * 
 * MX Cheminformatics Tools for Java
 *
 * Copyright (c) 2007-2009 Metamolecular, LLC
 *
 * http://metamolecular.com
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.openscience.cdk.smsd.algorithm.vflib.query;


import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.smsd.algorithm.vflib.builder.VFQueryBuilder;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IQuery;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IQueryCompiler;
import org.openscience.cdk.smsd.algorithm.vflib.validator.VFAtomMatcher;
import org.openscience.cdk.smsd.algorithm.vflib.validator.VFBondMatcher;

/**
 * @author Richard L. Apodaca <rapodaca at metamolecular.com>
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk> (modified the orignal code)
 * 
 * @cdk.module smsd
 */
public class TemplateCompiler implements IQueryCompiler {

//    private Reducer reducer;
//    private Map<IAtom, Integer> reductions;
    private IAtomContainer molecule;

    public TemplateCompiler() {
    }

    public static IQuery compile(IAtomContainer molecule) {
        TemplateCompiler compiler = new TemplateCompiler();

        compiler.setMolecule(molecule);

        return compiler.compile();
    }

    public void setMolecule(IAtomContainer molecule) {
        this.molecule = molecule;
    }

    public IAtomContainer getMolecule() {
        return molecule;
    }

    @Override
    public IQuery compile() {
//        try {
//            IAtomContainer copy = (IAtomContainer) queryMolecule.clone();
//        reductions.clear();
        //Remove Hydrogen if Necesarry Asad
//            reducer.reduce(copy, reductions);
        return build(molecule);
//        } catch (CloneNotSupportedException ex) {
//            Logger.getLogger(TemplateCompiler.class.getName()).log(Level.SEVERE, null, ex);
//        }

//        return null;
    }

    private IQuery build(IAtomContainer queryMolecule) {
        VFQueryBuilder result = new VFQueryBuilder();

        for (int i = 0; i < queryMolecule.getAtomCount(); i++) {
            IAtom atom = queryMolecule.getAtom(i);
            IQueryAtom matcher = createMatcher(atom);

            if (matcher != null) {
                result.addNode(matcher, atom);
            }
        }

        for (int i = 0; i < queryMolecule.getBondCount(); i++) {
            IBond bond = queryMolecule.getBond(i);

            IAtom s = bond.getAtom(0);
            IAtom t = bond.getAtom(1);

//            int IndexI = queryMolecule.getAtomNumber(s);
//            int IndexJ = queryMolecule.getAtomNumber(t);
//
//            System.out.println("\nIndexI " + IndexI);
//            System.out.println("IndexJ " + IndexJ);

//            result.connect(result.getNode(IndexI), result.getNode(IndexJ), createBondMatcher(bond));

            result.connect(result.getNode(s), result.getNode(t), createBondMatcher(bond));
        }

        return result;
    }

    private IQueryAtom createMatcher(IAtom atom) {
        return new VFAtomMatcher(atom);
    }

    private IQueryBond createBondMatcher(IBond bond) {
        return new VFBondMatcher(bond);
    }
}
