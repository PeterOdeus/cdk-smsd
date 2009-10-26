/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.core.tools;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK. e-mail: asad@ebi.ac.uk
 */
public class SubGraphMapper {

    /**
     *
     * @param mainMolecule
     * @param matchedCore
     * @return
     */
    public static IAtomContainer getNeedle(IAtomContainer mainMolecule, List<Integer> matchedCore) {
        IAtomContainer needle = DefaultChemObjectBuilder.getInstance().newAtomContainer();

        assert matchedCore != null;
//        System.out.println("Number of matched atoms = " + matchedCore.size());
        for (int i = 0; i < matchedCore.size(); i++) {
            for (int j = 0; j < matchedCore.size(); j++) {

                int iIndex = matchedCore.get(i);
                int jIndex = matchedCore.get(j);

                if (i != j && iIndex != jIndex) {


                    IAtom atom0 = mainMolecule.getAtom(iIndex);
                    IAtom atom1 = mainMolecule.getAtom(jIndex);

                    atom0.setID(String.valueOf(iIndex));
                    atom1.setID(String.valueOf(jIndex));

                    IBond bond = mainMolecule.getBond(atom0, atom1);
                    if (bond != null) {
                        needle.addBond(bond);

                    }
                }
            }
        }

        return needle;
    }

    public static Vector<List<Integer>> getNeedleList(List<Integer> matchedCore) {
        Vector<List<Integer>> matchedList = new Vector<List<Integer>>();
        List<Integer> rList = new LinkedList<Integer>();
        List<Integer> pList = new LinkedList<Integer>();

        for (int i = 0; i < matchedCore.size(); i += 2) {

            rList.add(matchedCore.get(i));
            pList.add(matchedCore.get(i + 1));

        }

        matchedList.add(0, rList);
        matchedList.add(1, pList);
        return matchedList;
    }

    public static Vector<List<Integer>> getNeedleList(Map<Integer, Integer> matchedCore) {
        Vector<List<Integer>> matchedList = new Vector<List<Integer>>();
        List<Integer> rList = new LinkedList<Integer>();
        List<Integer> pList = new LinkedList<Integer>();

        for (Map.Entry<Integer, Integer> M : matchedCore.entrySet()) {

            rList.add(M.getKey());
            pList.add(M.getValue());

        }

        matchedList.add(0, rList);
        matchedList.add(1, pList);
        return matchedList;
    }

//    /**
//     *
//     * @param ac1 Query Container
//     * @param ac2 Target Container
//     * @param mappings Atom mapping
//     * @return Bond mapping
//     */
//    public static Map<IBond, IBond> makeBondMapsOfAtomMaps(IAtomContainer ac1, IAtomContainer ac2, TreeMap<Integer, Integer> mappings) {
//
//        HashMap<IBond, IBond> maps = new HashMap<IBond, IBond>();
//
//        for (IAtom atoms : ac1.atoms()) {
//
//            int ac1AtomNumber = ac1.getAtomNumber(atoms);
//
//            if (mappings.containsKey(ac1AtomNumber)) {
//
//                int ac2AtomNumber = mappings.get(ac1AtomNumber);
//
//                List<IAtom> connectedAtoms = ac1.getConnectedAtomsList(atoms);
//
//                for (IAtom cAtoms : connectedAtoms) {
//                    int ac1ConnectedAtomNumber = ac1.getAtomNumber(cAtoms);
//
//                    if (mappings.containsKey(ac1ConnectedAtomNumber)) {
//                        {
//
//                            int ac2ConnectedAtomNumber = mappings.get(ac1ConnectedAtomNumber);
//
//                            IBond ac1Bond = ac1.getBond(atoms, cAtoms);
//
//                            IBond ac2Bond = ac2.getBond(ac2.getAtom(ac2AtomNumber), ac2.getAtom(ac2ConnectedAtomNumber));
//
//                            if (ac2Bond == null) {
//                                ac2Bond = ac2.getBond(ac2.getAtom(ac2ConnectedAtomNumber), ac2.getAtom(ac2AtomNumber));
//
//                            }
//
//                            if (ac1Bond != null && ac2Bond != null) {
//                                maps.put(ac1Bond, ac2Bond);
//                            }
//
//                        }
//
//
//                    }
//                }
//            }
//        }
//
//        return maps;
//
//    }
    public static Map<IBond, IBond> makeBondMapsOfAtomMaps(IAtomContainer ac1, IAtomContainer ac2, TreeMap<Integer, Integer> mappings) {

        HashMap<IBond, IBond> maps = new HashMap<IBond, IBond>();

        for (Map.Entry<Integer, Integer> mapS : mappings.entrySet()) {
            int indexI = mapS.getKey();
            int indexJ = mapS.getValue();

            for (Map.Entry<Integer, Integer> mapD : mappings.entrySet()) {

                int indexIPlus = mapD.getKey();
                int indexJPlus = mapD.getValue();

                if (indexI != indexIPlus && indexJ != indexJPlus) {

                    IAtom atomI0 = ac1.getAtom(indexI);
                    IAtom atomI1 = ac1.getAtom(indexIPlus);

                    atomI0.setID(String.valueOf(indexI));
                    atomI1.setID(String.valueOf(indexIPlus));

                    IBond ac1Bond = ac1.getBond(atomI0, atomI1);

                    if (ac1Bond != null) {

                        IAtom atomJ0 = ac1.getAtom(indexJ);
                        IAtom atomJ1 = ac1.getAtom(indexJPlus);

                        atomJ0.setID(String.valueOf(indexJ));
                        atomJ1.setID(String.valueOf(indexJPlus));

                        IBond ac2Bond = ac1.getBond(atomJ0, atomJ1);

                        if (ac2Bond != null) {

                            maps.put(ac1Bond, ac2Bond);
                        }

                    }

                }

            }
        }

//        System.out.println("bond Map size:" + maps.size());

        return maps;

    }
}
