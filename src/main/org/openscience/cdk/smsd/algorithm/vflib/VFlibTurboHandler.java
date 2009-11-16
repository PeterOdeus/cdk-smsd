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
package org.openscience.cdk.smsd.algorithm.vflib;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IMapper;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.INode;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IQuery;
import org.openscience.cdk.smsd.algorithm.vflib.map.VFMapper;
import org.openscience.cdk.smsd.algorithm.vflib.query.TemplateCompiler;
import org.openscience.cdk.smsd.core.tools.EBIAtomContainerManipulator;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.ISubGraph;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @cdk.module smsd
 */
public class VFlibTurboHandler implements ISubGraph {

    private IAtomContainer Reactant;
    private IAtomContainer Product;
    private static Vector<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static Vector<TreeMap<Integer, Integer>> allMCS = null;

    public VFlibTurboHandler() {


        allAtomMCS = new Vector<Map<IAtom, IAtom>>();
        atomsMCS = new HashMap<IAtom, IAtom>();
        firstMCS = new TreeMap<Integer, Integer>();
        allMCS = new Vector<TreeMap<Integer, Integer>>();

    }

    /**
     *
     * @return true if Query/Reactant is a subgraph of Target/Product
     * else false
     */
    @Override
    public boolean isSubgraph(){

        IQuery query = TemplateCompiler.compile(Reactant);

        IMapper mapper = new VFMapper(query);

        Map<INode, IAtom> vfLibSolutions = mapper.getFirstMap(Product);

        Map<IAtom, IAtom> atomatomMapping = new HashMap<IAtom, IAtom>();
        TreeMap<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();

        for (Map.Entry<INode, IAtom> mapping : vfLibSolutions.entrySet()) {

            IAtom qAtom = query.getAtom(mapping.getKey());
            IAtom tAtom = mapping.getValue();

            Integer qIndex = Reactant.getAtomNumber(qAtom);
            Integer tIndex = Product.getAtomNumber(tAtom);

            atomatomMapping.put(qAtom, tAtom);
            indexindexMapping.put(qIndex, tIndex);


        }
        if (atomatomMapping.size() > 0) {
            allAtomMCS.add(atomatomMapping);
            allMCS.add(indexindexMapping);
        }



        if (allAtomMCS.size() > 0) {
            atomsMCS.putAll(allAtomMCS.firstElement());
            firstMCS.putAll(allMCS.firstElement());
        }

        if (firstMCS.size() == Reactant.getAtomCount()) {
            return true;
        } else {
            return false;
        }

    }

    private int checkForH(IAtomContainer mol) {
        int hCount = 0;


        for (int i = 0; i <
                mol.getAtomCount(); i++) {


            if (mol.getAtom(i).getSymbol().equals("H")) {
                hCount++;
            }

        }


        return hCount;
    }

    /**
     * Set the JMCS software
     *
     * @param reactant
     * @param product
     * @param removeHydrogen
     */
    @Override
    public void set(IAtomContainer reactant, IAtomContainer product, boolean removeHydrogen) {

        this.Reactant = reactant;
        this.Product = product;

        /*Remove Hydrogen by Asad*/
        if (checkForH(Reactant) > 0 && removeHydrogen) {
            Reactant = EBIAtomContainerManipulator.removeHydrogens(reactant);
        }
        if (checkForH(Product) > 0 && removeHydrogen) {
            Product = EBIAtomContainerManipulator.removeHydrogens(product);
        }

    }

    /**
     * Set the VFLib software
     *
     * @param reactant
     * @param product
     * @param removeHydrogen
     */
    @Override
    public void set(MolHandler reactant, MolHandler product, boolean removeHydrogen) {


        this.Reactant = reactant.getMolecule();
        this.Product = product.getMolecule();

        /*Remove Hydrogen by Asad*/
        if (checkForH(Reactant) > 0 && removeHydrogen) {
            Reactant = EBIAtomContainerManipulator.removeHydrogens(reactant.getMolecule());
        }
        if (checkForH(Product) > 0 && removeHydrogen) {
            Product = EBIAtomContainerManipulator.removeHydrogens(product.getMolecule());
        }


    }

    /**
     * Creates a new instance of SearchCliques
     * @param ReactantMolFileName
     * @param ProductMolFileName
     * @param removeHydrogen
     */
    @Override
    public void set(String ReactantMolFileName, String ProductMolFileName, boolean removeHydrogen) {


        String mol1 = ReactantMolFileName;
        String mol2 = ProductMolFileName;

        this.Reactant = new MolHandler(mol1, false, removeHydrogen).getMolecule();
        this.Product = new MolHandler(mol2, false, removeHydrogen).getMolecule();

        if (checkForH(Reactant) > 0 && removeHydrogen) {
            Reactant = EBIAtomContainerManipulator.removeHydrogens(Reactant);
        }
        if (checkForH(Product) > 0 && removeHydrogen) {
            Product = EBIAtomContainerManipulator.removeHydrogens(Product);
        }

    }

    /**
     *
     * @return
     */
    @Override
    public Vector<Map<IAtom, IAtom>> getAllAtomMapping() {
        return allAtomMCS;
    }

    @Override
    public Vector<TreeMap<Integer, Integer>> getAllMapping() {
        return allMCS;
    }

    @Override
    public Map<IAtom, IAtom> getFirstAtomMapping() {
        return atomsMCS;
    }

    @Override
    public TreeMap<Integer, Integer> getFirstMapping() {
        return firstMCS;
    }
}
