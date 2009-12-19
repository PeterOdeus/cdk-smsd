/* Copyright (C) 2009  Syed Asad Rahman {asad@ebi.ac.uk}
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
 * You should have received commonAtomList copy of the GNU Lesser General Public License
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
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smsd.algorithm.mcgregor.McGregor;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IMapper;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.INode;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IQuery;
import org.openscience.cdk.smsd.algorithm.vflib.map.VFMCSMapper;
import org.openscience.cdk.smsd.algorithm.vflib.query.TemplateCompiler;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.IMCSAlgorithm;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @cdk.module smsd
 */
public class VFlibMCSHandler implements IMCSAlgorithm {

    private static List<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static List<Map<IAtom, IAtom>> allAtomMCSCopy = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static List<TreeMap<Integer, Integer>> allMCS = null;
    private static List<TreeMap<Integer, Integer>> allMCSCopy = null;
    private IAtomContainer ac1 = null;
    private IAtomContainer ac2 = null;
    private List<Map<INode, IAtom>> vfLibSolutions = null;
    private int vfMCSSize = 0;

    public VFlibMCSHandler() {
        allAtomMCS = new ArrayList<Map<IAtom, IAtom>>();
        allAtomMCSCopy = new ArrayList<Map<IAtom, IAtom>>();
        atomsMCS = new HashMap<IAtom, IAtom>();
        firstMCS = new TreeMap<Integer, Integer>();
        allMCS = new ArrayList<TreeMap<Integer, Integer>>();
        allMCSCopy = new ArrayList<TreeMap<Integer, Integer>>();
    }

    /**
     *
     * @return true if Query/Reactant is atoms subgraph of Target/Product
     * else false
     * @throws java.io.IOException
     */
    @Override
    public int searchMCS() throws IOException, CDKException {
        searchVFMCSMappings();
        boolean flag = mcgregorFlag();
        if (flag) {
            searchMcGregorMapping();
        } else if (!allAtomMCSCopy.isEmpty()) {
            allAtomMCS.addAll(allAtomMCSCopy);
            allMCS.addAll(allMCSCopy);
        }
        setFirstMappings();
        return 0;

    }

    private void setFirstMappings() {
        if (!allAtomMCS.isEmpty()) {
            atomsMCS.putAll(allAtomMCS.get(0));
            firstMCS.putAll(allMCS.get(0));
        }

    }

    private boolean mcgregorFlag() {
        int commonAtomCount = checkCommonAtomCount(ac1, ac2);
        if (commonAtomCount > vfMCSSize && commonAtomCount > vfMCSSize) {
            return true;

        } else {
            return false;
        }

    }

    /**
     * Set the VFLibMCS software
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
        this.set(Reactant, Product);

    }

    /**
     * Set the VFLib MCS software
     *
     * @param Reactant
     * @param Product 
     */
    @Override
    public void set(MolHandler Reactant, MolHandler Product) {
        ac1 = Reactant.getMolecule();
        ac2 = Product.getMolecule();
    }

    /**
     * Creates atoms new instance of SearchCliques
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

    /**
     * Set the VFLib MCS software
     *
     * @param source
     * @param target 
     */
    public void set(IMolecule source, IMolecule target) throws CDKException {
        MolHandler Reactant = new MolHandler(source, false);
        MolHandler Product = new MolHandler(target, false);
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
        ArrayList<String> atoms = new ArrayList<String>();
        for (int i = 0; i <
                reactantMolecule.getAtomCount(); i++) {
            atoms.add(reactantMolecule.getAtom(i).getSymbol());

        }
        int common = 0;
        for (int i = 0; i < productMolecule.getAtomCount(); i++) {
            String symbol = productMolecule.getAtom(i).getSymbol();
            if (atoms.contains(symbol)) {
                atoms.remove(symbol);
                common++;

            }
        }
        return common;
    }

    private void searchVFMCSMappings() {
        IQuery query = null;
        IMapper mapper = null;
        boolean RONP = false;
        if (ac1.getAtomCount() <= ac2.getAtomCount()) {
            query = TemplateCompiler.compile(ac1);
            mapper = new VFMCSMapper(query);
            vfLibSolutions = new ArrayList<Map<INode, IAtom>>(mapper.getMaps(ac2));
            RONP = true;

        } else {
            query = TemplateCompiler.compile(ac2);
            mapper = new VFMCSMapper(query);
            vfLibSolutions = new ArrayList<Map<INode, IAtom>>(mapper.getMaps(ac1));
            RONP = false;
        }

        setVFMCSMappings(RONP, query);
        this.vfMCSSize = allMCSCopy.isEmpty() ? 0 : allMCSCopy.get(0).size();
    }

    private void searchMcGregorMapping() throws CDKException, IOException {
        List<List<Integer>> _mappings = new ArrayList<List<Integer>>();

        for (TreeMap<Integer, Integer> firstPassMappings : allMCSCopy) {
            McGregor mgit = new McGregor(ac1, ac2, _mappings);
            mgit.startMcGregorIteration(mgit.getMCSSize(), firstPassMappings); //Start McGregor search
            _mappings = mgit.getMappings();
            mgit = null;
        }
        setMcGregorMappings(_mappings);
    }

    private void setVFMCSMappings(boolean RONP, IQuery query) {
        int counter = 0;

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
                allAtomMCSCopy.add(counter, atomatomMapping);
                allMCSCopy.add(counter, indexindexMapping);
                counter++;
            }
        }
    }

    private void setMcGregorMappings(List<List<Integer>> _mappings) throws CDKException {
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

            if (!atomatomMapping.isEmpty() && !hasMap(indexindexMapping, allMCS)) {
                allAtomMCS.add(counter, atomatomMapping);
                allMCS.add(counter, indexindexMapping);
                counter++;
            }
        }
    }
}
