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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.smsd.algorithm.mcgregor.McGregor;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IMapper;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.INode;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IQuery;
import org.openscience.cdk.smsd.algorithm.vflib.map.VFMCSMapper;
import org.openscience.cdk.smsd.algorithm.vflib.query.TemplateCompiler;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.IMCS;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @cdk.module smsd
 */
public class VFlibMCSHandler implements IMCS {

    private static List<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static List<Map<IAtom, IAtom>> allAtomMCS_copy = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static List<TreeMap<Integer, Integer>> allMCS = null;
    private static List<TreeMap<Integer, Integer>> allMCS_copy = null;
    private IAtomContainer ac1 = null;
    private IAtomContainer ac2 = null;
    private List<Map<INode, IAtom>> vfLibSolutions = null;
    private int VFMCSSize = 0;

    public VFlibMCSHandler() {


        allAtomMCS = new Vector<Map<IAtom, IAtom>>();
        allAtomMCS_copy = new Vector<Map<IAtom, IAtom>>();
        atomsMCS = new HashMap<IAtom, IAtom>();
        firstMCS = new TreeMap<Integer, Integer>();
        allMCS = new Vector<TreeMap<Integer, Integer>>();

        allMCS_copy = new Vector<TreeMap<Integer, Integer>>();


    }

    /**
     *
     * @return true if Query/Reactant is a subgraph of Target/Product
     * else false
     * @throws java.io.IOException
     * @throws CDKException 
     */
    @Override
    public int searchMCS(boolean removeHydrogen) throws IOException, CDKException {


//        System.out.println("VF Solution Count: " + vfLibSolutions.size());


//
//        System.out.println("count of the VFLib Mapping Stored: " + allAtomMCS_copy.size());
//        System.out.println(
//                "Mapping Size of the VFLib Mapping Stored: " + allAtomMCS_copy.firstElement().size());

        boolean flag = mcgregorFlag();



        if (flag) {
            List<List<Integer>> _mappings = new Vector<List<Integer>>();

            for (TreeMap<Integer, Integer> firstPassMappings : allMCS_copy) {
//
//                System.out.println("\n\n---------------\n");
//                System.out.println("Calling McGregor");
//                System.out.println("for Start size " + firstPassMappings.size());
                McGregor mgit = new McGregor(ac1, ac2, _mappings);
                mgit.startMcGregorIteration(mgit.getMCSSize(), firstPassMappings); //Start McGregor search

                _mappings = mgit.getMappings();
//                System.out.println("");
                mgit = null;
//                System.out.println("After Calling McGregor, solution Size" + _mappings.firstElement().size() / 2);
//

            }
//            System.out.println("Solution Count After McGregor " + _mappings.size());
//            System.out.println("Solution Size After Calling McGregor" + _mappings.firstElement().size() / 2);

            int counter = 0;
            for (List<Integer> mapping : _mappings) {

                Map<IAtom, IAtom> atomatomMapping = new HashMap<IAtom, IAtom>();
                TreeMap<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();

                for (int index = 0; index < mapping.size(); index += 2) {
                    IAtom qAtom = null;
                    IAtom tAtom = null;

                    qAtom = ac1.getAtom(mapping.get(index));
                    tAtom = ac2.getAtom(mapping.get(index + 1));


                    Integer qIndex = mapping.get(index);
                    Integer tIndex = mapping.get(index + 1);


                    if (qIndex != null && tIndex != null) {
                        atomatomMapping.put(qAtom, tAtom);
                        indexindexMapping.put(qIndex, tIndex);
                    } else {

                        throw new CDKException("Atom index pointing to NULL");
                    }
                }

                if (!atomatomMapping.isEmpty()) {
                    if (!hasMap(indexindexMapping, allMCS)) {
                        allAtomMCS.add(counter, atomatomMapping);
                        allMCS.add(counter, indexindexMapping);
                        counter++;
                    }
                }

            }


        } else if (!allAtomMCS_copy.isEmpty()) {


            allAtomMCS.addAll(allAtomMCS_copy);
            allMCS.addAll(allMCS_copy);
        }



        if (!allAtomMCS.isEmpty()) {
            atomsMCS.putAll(allAtomMCS.get(0));
            firstMCS.putAll(allMCS.get(0));


        }

//        System.out.println("Size of the Atom Mapping Stored: " + allAtomMCS.firstElement().size());
//        System.out.println("Size of the Atom Index Mapping Stored: " + allMCS.firstElement().size());




        return 0;

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

    private boolean mcgregorFlag() {
        int commonAtomCount = checkCommonAtomCount(ac1, ac2);
//        System.out.println("ra: " + ratomCount);
//        System.out.println("pa: " + patomCount);
//        System.out.println("ca: " + commonAtomCount);

        if (commonAtomCount > VFMCSSize && commonAtomCount > VFMCSSize) {
            return true;

        } else {
            return false;
        }
//       

    }

    /**
     * Set the JMCS software
     *
     * @param reactant
     * @param product
     */
    @Override
    public void set(IAtomContainer reactant, IAtomContainer product) {

        IAtomContainer mol1 = reactant;
        IAtomContainer mol2 = product;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);

        //System.out.println(atom_number1+" " + atom_num_H_1);

        //best_clique_size = 0;

        this.set(Reactant, Product);
        // System.exit(1);

    }

    /**
     * Set the JMCS software
     *
     */
    @Override
    public void set(MolHandler Reactant, MolHandler Product) {

        IQuery query = null;
        IMapper mapper = null;
        boolean RONP = false;
        ac1 = Reactant.getMolecule();
        ac2 = Product.getMolecule();



//        System.out.println("R Atom count " + ac1.getAtomCount());
//        System.out.println("P Atom count " + ac2.getAtomCount());

        if (ac1.getAtomCount() <= ac2.getAtomCount()) {

//            System.out.println("Query is Reactant");

            query = TemplateCompiler.compile(ac1);
            mapper = new VFMCSMapper(query);
            vfLibSolutions = new ArrayList<Map<INode, IAtom>>(mapper.getMaps(ac2));
            RONP = true;

        } else {
//            System.out.println("Query is Prouct");
            query = TemplateCompiler.compile(ac2);
            mapper = new VFMCSMapper(query);
            vfLibSolutions = new ArrayList<Map<INode, IAtom>>(mapper.getMaps(ac1));
            RONP = false;
        }

        for (Map<INode, IAtom> solution : vfLibSolutions) {

            Map<IAtom, IAtom> atomatomMapping = new HashMap<IAtom, IAtom>();
            TreeMap<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();

            for (Map.Entry<INode, IAtom> mapping : solution.entrySet()) {
                IAtom qAtom = null;
                IAtom tAtom = null;
                if (RONP) {
                    qAtom = query.getAtom(mapping.getKey());
                    tAtom = mapping.getValue();

                } else {
                    tAtom = query.getAtom(mapping.getKey());
                    qAtom = mapping.getValue();
                }

                Integer qIndex = Integer.valueOf(ac1.getAtomNumber(qAtom));
                Integer tIndex = Integer.valueOf(ac2.getAtomNumber(tAtom));

//
//                System.out.println("i:" + qIndex + " j:" + tIndex);
//                System.out.println("i:" + qAtom.getSymbol() + " j:" + tAtom.getSymbol());

                if (qIndex != null && tIndex != null) {
                    atomatomMapping.put(qAtom, tAtom);
                    indexindexMapping.put(qIndex, tIndex);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to NULL");
                    } catch (CDKException ex) {
                        Logger.getLogger(VFlibMCSHandler.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

            }
            if (!atomatomMapping.isEmpty()) {
                allAtomMCS_copy.add(atomatomMapping);
                allMCS_copy.add(indexindexMapping);
                this.VFMCSSize = atomatomMapping.size();
            }


        }


    }

    /**
     * Creates a new instance of SearchCliques
     * @param ReactantMolFileName
     * @param ProductMolFileName
     */
    @Override
    public void set(String ReactantMolFileName, String ProductMolFileName) {

        String mol1 = ReactantMolFileName;
        String mol2 = ProductMolFileName;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);
        this.set(Reactant, Product);


    }

    private boolean hasMap(Map<Integer, Integer> map, List<TreeMap<Integer, Integer>> mapGlobal) {

        for (Map<Integer, Integer> test : mapGlobal) {

            if (test.equals(map)) {
                return true;
            }

        }


        return false;
    }

    /**
     *
     * @return
     */
    @Override
    public List<Map<IAtom, IAtom>> getAllAtomMapping() {
        return allAtomMCS;
    }

    @Override
    public List<TreeMap<Integer, Integer>> getAllMapping() {
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

    private int checkCommonAtomCount(IAtomContainer reactantMolecule, IAtomContainer productMolecule) {
        ArrayList<String> a = new ArrayList<String>();
        for (int i = 0; i <
                reactantMolecule.getAtomCount(); i++) {
            a.add(reactantMolecule.getAtom(i).getSymbol());

        }


        int common = 0;
        for (int i = 0; i <
                productMolecule.getAtomCount(); i++) {

            if (a.contains(productMolecule.getAtom(i).getSymbol())) {
                a.remove(productMolecule.getAtom(i).getSymbol());
                common++;

            }


        }
        return common;
    }
}
