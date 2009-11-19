/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.tools;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import java.util.Map;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.interfaces.IAtomParity;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.ILonePair;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.mcss.RMap;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

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
 * <p>This is an extension of CDK AtomContainer.
 * Some part of this code was taken from CDK source code and modified.</p>
 * @cdk.module smsd
 */
public class ExtAtomContainerManipulator extends AtomContainerManipulator {

    /**
     * 
     * @param container
     * @return
     */
    public static IMolecule makeDeepCopy(IAtomContainer container) {

        IMolecule newAtomContainer = DefaultChemObjectBuilder.getInstance().newMolecule();

        int atomCount = container.getAtomCount();
        int bondCount = container.getBondCount();
        int lonePairCount = container.getLonePairCount();
        int singleElectronCount = container.getSingleElectronCount();
        IAtom[] atoms = new IAtom[atomCount];
        IBond[] bonds = new IBond[bondCount];
        ILonePair[] lonePairs = new ILonePair[lonePairCount];
        ISingleElectron[] singleElectrons = new ISingleElectron[singleElectronCount];

//        Deep copy of the Atoms

        for (int f = 0; f < container.getAtomCount(); f++) {

            if (container.getAtom(f) instanceof PseudoAtom) {

                atoms[f] = new PseudoAtom(container.getAtom(f));
            } else {
                atoms[f] = new Atom(container.getAtom(f));
            }


            if ((container.getAtom(f)).getPoint2d() != null) {
                atoms[f].setPoint2d(new Point2d(container.getAtom(f).getPoint2d()));
            }
            if ((container.getAtom(f)).getPoint3d() != null) {
                atoms[f].setPoint3d(new Point3d(container.getAtom(f).getPoint3d()));
            }
            if ((container.getAtom(f)).getFractionalPoint3d() != null) {
                atoms[f].setFractionalPoint3d(new Point3d(container.getAtom(f).getFractionalPoint3d()));
            }

            if (container.getAtom(f).getID() != null) {
                atoms[f].setID(new String(container.getAtom(f).getID()));
            }
            if (container.getAtom(f).getHydrogenCount() != null) {
                atoms[f].setHydrogenCount(new Integer(container.getAtom(f).getHydrogenCount()));
            }
            if (container.getAtom(f).getCharge() != null) {
                atoms[f].setCharge(new Double(container.getAtom(f).getCharge()));
            }
            if (container.getAtom(f).getStereoParity() != null) {
                atoms[f].setStereoParity(new Integer(container.getAtom(f).getStereoParity()));
            }

            newAtomContainer.addAtom(atoms[f]);

            if (container.getAtomParity(container.getAtom(f)) != null) {
                IAtomParity ap = container.getAtomParity(container.getAtom(f));
                newAtomContainer.addAtomParity(ap);
            }

        }

//  Deep copy of the bonds
        for (int f = 0; f < container.getBondCount(); f++) {


            bonds[f] = new Bond();


            int IndexI = 999;
            for (int i = 0; i < container.getAtomCount(); i++) {
                if (container.getBond(f).getAtom(0) == container.getAtom(i)) {
                    IndexI = i;
                    break;
                }
            }
            int IndexJ = 999;
            for (int j = 0; j < container.getAtomCount(); j++) {
                if (container.getBond(f).getAtom(1) == container.getAtom(j)) {
                    IndexJ = j;
                    break;
                }
            }

            IAtom atom1 = atoms[IndexI];
            IAtom atom2 = atoms[IndexJ];

            Order order = container.getBond(f).getOrder();
            int stereo = container.getBond(f).getStereo();
            bonds[f] = new Bond(atom1, atom2, order, stereo);
            if (container.getBond(f).getID() != null) {
                bonds[f].setID(new String(container.getBond(f).getID()));
            }
            newAtomContainer.addBond(bonds[f]);

        }

// Deep copy of the LonePairs

        for (int f = 0; f < container.getLonePairCount(); f++) {

            if (container.getAtom(f).getSymbol().equalsIgnoreCase("R")) {
                lonePairs[f] = DefaultChemObjectBuilder.getInstance().newLonePair(container.getAtom(f));
            }
            newAtomContainer.addLonePair(lonePairs[f]);

        }

        for (int f = 0; f < container.getSingleElectronCount(); f++) {
            singleElectrons[f] = DefaultChemObjectBuilder.getInstance().newSingleElectron(container.getAtom(f));
            newAtomContainer.addSingleElectron(singleElectrons[f]);

        }
        newAtomContainer.setProperties(container.getProperties());
        newAtomContainer.setFlags(container.getFlags());

        newAtomContainer.setID(container.getID());

        newAtomContainer.notifyChanged();
        return newAtomContainer;

    }

   
    /**
     *  The method takes an xml files like the following:<br>
     *  &lt;replace-set&gt;<br>
     *  &lt;replace&gt;O=N=O&lt;/replace&gt;<br>
     *  &lt;replacement&gt;[O-][N+]=O&lt;/replacement&gt;<br>
     *  &lt;/replace-set&gt;<br>
     *  All parts in ac which are the same as replace will be changed according to replacement.
     *  Currently the following changes are done: BondOrder, FormalCharge.
     *  For detection of fragments like replace, we rely on UniversalIsomorphismTester.
     *  doc may contain several replace-sets and a replace-set may contain several replace fragments, which will all be normalized according to replacement.
     *
     * @param  ac                          The atomcontainer to normalize.
     * @param  doc                         The configuration file.
     * @return                             Did a replacement take place?
     * @exception  InvalidSmilesException  doc contains an invalid smiles.
     * @throws CDKException 
     */
    public static boolean normalize(IAtomContainer ac, Document doc) throws InvalidSmilesException, CDKException {
        NodeList nl = doc.getElementsByTagName("replace-set");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        boolean change = false;
        for (int i = 0; i < nl.getLength(); i++) {
            Element child = (Element) nl.item(i);
            NodeList replaces = child.getElementsByTagName("replace");
            NodeList replacement = child.getElementsByTagName("replacement");
            String replacementstring = replacement.item(0).getFirstChild().getNodeValue();
            if (replacementstring.indexOf("\n") > -1 || replacementstring.length() < 1) {
                replacementstring = replacement.item(0).getFirstChild().getNextSibling().getNodeValue();
            }
            IAtomContainer replacementStructure = sp.parseSmiles(replacementstring);
            for (int k = 0; k < replaces.getLength(); k++) {
                Element replace = (Element) replaces.item(k);
                String replacestring = replace.getFirstChild().getNodeValue();
                if (replacestring.indexOf("\n") > -1 || replacestring.length() < 1) {
                    replacestring = replace.getFirstChild().getNextSibling().getNodeValue();
                }
                IAtomContainer replaceStructure = sp.parseSmiles(replacestring);
                List l = null;
                while ((l = UniversalIsomorphismTester.getSubgraphMap(ac, replaceStructure)) != null) {
                    @SuppressWarnings("unchecked")
                    List<RMap> l2 = UniversalIsomorphismTester.makeAtomsMapOfBondsMap(l, ac, replaceStructure);
                    Iterator bondit = l.iterator();
                    while (bondit.hasNext()) {
                        RMap rmap = (RMap) bondit.next();
                        org.openscience.cdk.interfaces.IBond acbond = ac.getBond(rmap.getId1());
                        org.openscience.cdk.interfaces.IBond replacebond = replacementStructure.getBond(rmap.getId2());
                        acbond.setOrder(replacebond.getOrder());
                        change = true;
                    }
                    Iterator atomit = l2.iterator();
                    while (atomit.hasNext()) {
                        RMap rmap = (RMap) atomit.next();
                        org.openscience.cdk.interfaces.IAtom acatom = ac.getAtom(rmap.getId1());
                        org.openscience.cdk.interfaces.IAtom replaceatom = replacementStructure.getAtom(rmap.getId2());
                        acatom.setFormalCharge(replaceatom.getFormalCharge());
                        change = true;
                    }
                }
            }
        }
        return (change);
    }

    public static String fixSmiles(String Smiles) {
        Smiles = Smiles.replaceAll("CL", "Cl");
        Smiles = Smiles.replaceAll("(H)", "([H])");
//		Smiles=Smiles.replace("N=N#N","N=[N+]=[N-]");
//		Smiles=Smiles.replace("#N=O","#[N+][O-]");
        Smiles = Smiles.trim();

        return Smiles;

    }

    private static boolean fixNitroGroups(IMolecule m) {
        // changes nitros given by N(=O)(=O) to [N+](=O)[O-]
        boolean changed = false;
        try {
            for (int i = 0; i <= m.getAtomCount() - 1; i++) {
                IAtom a = m.getAtom(i);
//                boolean nitro = false;

                if (a.getSymbol().equals("N")) {
                    List ca = m.getConnectedAtomsList(a);

                    if (ca.size() == 3) {

                        IAtom[] cao = new IAtom[2];

                        int count = 0;

                        for (int j = 0; j <= 2; j++) {
                            if (((IAtom) ca.get(j)).getSymbol().equals("O")) {
                                count++;
                            }
                        }

                        if (count > 1) {

                            count = 0;
                            for (int j = 0; j <= 2; j++) {
                                IAtom caj = (IAtom) ca.get(j);
                                if ((caj.getSymbol().equals("O")) && (m.getConnectedAtomsCount(caj) == 1)) {// account for possibility of ONO2
                                    cao[count] = caj;
                                    count++;
                                }
                            }


                            IBond.Order order1 = m.getBond(a, cao[0]).getOrder();
                            IBond.Order order2 = m.getBond(a, cao[1]).getOrder();


                            //if (totalobonds==4) { // need to fix (FIXME)
                            if (order1 == IBond.Order.SINGLE &&
                                    order2 == IBond.Order.DOUBLE) {
                                a.setFormalCharge(1);
                                cao[0].setFormalCharge(-1); // pick first O arbitrarily
                                m.getBond(a, cao[0]).setOrder(IBond.Order.SINGLE);
                                changed = true;
                            }
                        } //else if (count==1) {// end if count>1

                    }// end ca==3 if

                } // end symbol == N


            }

            return changed;


        } catch (Exception e) {
            return changed;
        }

    }

    public static boolean fixNitroGroups2(IMolecule m) {
        // changes nitros given by [N+](=O)[O-] to N(=O)(=O) 
        boolean changed = false;
        try {
            for (int i = 0; i <= m.getAtomCount() - 1; i++) {
                IAtom a = m.getAtom(i);
//                boolean nitro = false;

                if (a.getSymbol().equals("N")) {
                    List ca = m.getConnectedAtomsList(a);

                    if (ca.size() == 3) {

                        IAtom[] cao = new IAtom[2];

                        int count = 0;

                        for (int j = 0; j <= 2; j++) {
                            IAtom caj = (IAtom) ca.get(j);
                            if (caj.getSymbol().equals("O")) {
                                count++;
                            }
                        }

                        if (count > 1) {

                            count = 0;
                            for (int j = 0; j <= 2; j++) {
                                IAtom caj = (IAtom) ca.get(j);
                                if ((caj.getSymbol().equals("O")) && (m.getConnectedAtomsCount(caj) == 1)) {// account for possibility of ONO2
                                    cao[count] = caj;
                                    count++;
                                }
                            }
                            IBond.Order order1 = m.getBond(a, cao[0]).getOrder();
                            IBond.Order order2 = m.getBond(a, cao[1]).getOrder();

                            //int totalobonds=0;						
                            //totalobonds+=m.getBond(a,cao[0]).getOrder();
//						totalobonds+=m.getBond(a,cao[1]).getOrder();

                            //if (totalobonds==4) { // need to fix
                            if ((order1 == IBond.Order.SINGLE && order2 == IBond.Order.DOUBLE) ||
                                    (order1 == IBond.Order.DOUBLE && order2 == IBond.Order.SINGLE)) {
                                a.setFormalCharge(0);
                                cao[0].setFormalCharge(0); // pick first O arbitrarily
                                cao[1].setFormalCharge(0); // pick first O arbitrarily
                                m.getBond(a, cao[0]).setOrder(IBond.Order.DOUBLE);
                                m.getBond(a, cao[1]).setOrder(IBond.Order.DOUBLE);
                                changed = true;
                            }
                        } // end if count>1


                    }// end ca==3 if

                } // end symbol == N


            }

            return changed;
        } catch (Exception e) {
            return changed;
        }
    }

    /**
     * 
     * @param m Molecule to be aromatized
     */
    public static void aromatizeMolecule(IAtomContainer m) {

        // need to find rings and aromaticity again since added H's

        IRingSet rs = null;
        try {
            AllRingsFinder arf = new AllRingsFinder();
            rs = arf.findAllRings(m);

            // SSSRFinder s = new SSSRFinder(m);
            // srs = s.findEssentialRings();

        } catch (Exception e) {
            e.printStackTrace();
        }

        try {
            // figure out which atoms are in aromatic rings:
//            printAtoms(m);
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(m);
//            printAtoms(m);
            CDKHueckelAromaticityDetector.detectAromaticity(m);
//            printAtoms(mol);
            // figure out which rings are aromatic:
            RingSetManipulator.markAromaticRings(rs);
//            printAtoms(mol);
            // figure out which simple (non cycles) rings are aromatic:
            // HueckelAromaticityDetector.detectAromaticity(m, srs);
        } catch (Exception e) {
            e.printStackTrace();
        }


        // only atoms in 6 membered rings are aromatic
        // determine largest ring that each atom is a part of

        for (int i = 0; i <= m.getAtomCount() - 1; i++) {

            m.getAtom(i).setFlag(CDKConstants.ISAROMATIC, false);

            jloop:
            for (int j = 0; j <= rs.getAtomContainerCount() - 1; j++) {
                //logger.debug(i+"\t"+j);
                IRing r = (IRing) rs.getAtomContainer(j);
                if (!r.getFlag(CDKConstants.ISAROMATIC)) {
                    continue jloop;
                }

                boolean haveatom = r.contains(m.getAtom(i));

                //logger.debug("haveatom="+haveatom);

                if (haveatom && r.getAtomCount() == 6) {
                    m.getAtom(i).setFlag(CDKConstants.ISAROMATIC, true);
                }

            }

        }


    }

    public static void fixSulphurH(IMolecule m) {
        // removes extra H's attached to sulphurs
        //logger.debug("EnterFixSulphur");

        for (int i = 0; i <= m.getAtomCount() - 1; i++) {
            IAtom a = m.getAtom(i);

            if (a.getSymbol().equals("S")) {
                List connectedAtoms = m.getConnectedAtomsList(a);


                int bondOrderSum = 0;

//                double oldBondOrderSum = m.getBondOrderSum(a); // includes H's

                for (int j = 0; j < connectedAtoms.size(); j++) {
                    IAtom conAtom = (IAtom) connectedAtoms.get(j);
                    if (!conAtom.getSymbol().equals("H")) {
                        IBond bond = m.getBond(a, conAtom);
                        if (bond.getOrder() == IBond.Order.SINGLE) {
                            bondOrderSum += 1;
                        } else if (bond.getOrder() == IBond.Order.DOUBLE) {
                            bondOrderSum += 2;
                        } else if (bond.getOrder() == IBond.Order.TRIPLE) {
                            bondOrderSum += 3;
                        } else if (bond.getOrder() == IBond.Order.QUADRUPLE) {
                            bondOrderSum += 4;
                        }
                    }
                }

                if (bondOrderSum > 1) {
                    for (int j = 0; j < connectedAtoms.size(); j++) {
                        IAtom conAtom = (IAtom) connectedAtoms.get(j);
                        if (conAtom.getSymbol().equals("H")) {
                            m.removeAtomAndConnectedElectronContainers(conAtom);
                        }
                    }
                }

            }

        }
    }

    /**
     * @param atomContainer
     * @param atom
     * @return The number of explicit hydrogens on the given IAtom.
     */
    public static int getExplicitHydrogenCount(IAtomContainer atomContainer, IAtom atom) {
        int hCount = 0;
        for (IAtom iAtom : atomContainer.getConnectedAtomsList(atom)) {
            IAtom connectedAtom = iAtom;
            if (connectedAtom.getSymbol().equals("H")) {
                hCount++;
            }
        }
        return hCount;
    }

    /**
     * @param atomContainer
     * @param atom
     * @return The summed implicit + explicit hydrogens of the given IAtom.
     */
    public static int getTotalHydrogenCount(IAtomContainer atomContainer, IAtom atom) {
        int hCount = atom.getHydrogenCount() == CDKConstants.UNSET ? 0 : atom.getHydrogenCount();
        hCount += getExplicitHydrogenCount(atomContainer, atom);
        return hCount;
    }


    /**
     *
     * @param atomContainer
     * @note function added by Asad
     * @return IAtomContainer without Hydrogen. If an AtomContainer has a single atom which
     * is a Hydrogen then its not removed.
     * @throws Exception
     */
    public static IAtomContainer removeHydrogensAndPreserveAtomID(IAtomContainer atomContainer) throws Exception {
        Map<IAtom, IAtom> map = new HashMap<IAtom, IAtom>();        // maps original atoms to clones.
        List<IAtom> remove = new ArrayList<IAtom>();  // lists removed Hs.
        IMolecule mol = null;
        if (atomContainer.getBondCount() > 0) {
            // Clone atoms except those to be removed.
            mol = atomContainer.getBuilder().newMolecule();
            int count = atomContainer.getAtomCount();
            for (int i = 0;
                    i < count;
                    i++) {
                // Clone/remove this atom?
                IAtom atom = atomContainer.getAtom(i);
                if (!atom.getSymbol().equals("H")) {
                    IAtom clonedAtom = null;
                    try {
                        clonedAtom = (IAtom) atom.clone();
                    } catch (CloneNotSupportedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                    //added by Asad to preserve the Atom ID for atom mapping without Hydrogen
                    clonedAtom.setID(atom.getID());
                    clonedAtom.setFlags(atom.getFlags());
                    int countH = 0;
                    if (atom.getHydrogenCount() != null) {
                        countH = atom.getHydrogenCount();
                    }
                    clonedAtom.setHydrogenCount(countH);
                    mol.addAtom(clonedAtom);
                    map.put(atom, clonedAtom);

                } else {
                    remove.add(atom);   // maintain list of removed H.
                }
            }

            // Clone bonds except those involving removed atoms.
            count = atomContainer.getBondCount();
            for (int i = 0;
                    i < count;
                    i++) {
                // Check bond.
                final IBond bond = atomContainer.getBond(i);
                boolean removedBond = false;
                final int length = bond.getAtomCount();
                for (int k = 0;
                        k < length;
                        k++) {
                    if (remove.contains(bond.getAtom(k))) {
                        removedBond = true;
                        break;
                    }
                }

                // Clone/remove this bond?
                if (!removedBond) // if (!remove.contains(atoms[0]) && !remove.contains(atoms[1]))
                {
                    IBond clone = null;
                    try {
                        clone = (IBond) atomContainer.getBond(i).clone();
                    } catch (CloneNotSupportedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                    assert clone != null;
                    clone.setAtoms(new IAtom[]{map.get(bond.getAtom(0)), map.get(bond.getAtom(1))});
                    mol.addBond(clone);
                }
            }



            // Recompute hydrogen counts of neighbours of removed Hydrogens.
            for (IAtom aRemove : remove) {
                // Process neighbours.
                for (IAtom iAtom : atomContainer.getConnectedAtomsList(aRemove)) {
                    final IAtom neighb = map.get(iAtom);
                    if (neighb == null) {
                        continue; // since for the case of H2, neight H has a heavy atom neighbor
                    }
                    //Added by Asad
                    if (!(neighb instanceof PseudoAtom)) {
                        neighb.setHydrogenCount(
                                (neighb.getHydrogenCount() == null ? 0 : neighb.getHydrogenCount()) + 1);
                    } else {
                        neighb.setHydrogenCount(0);
                    }
                }
            }
            mol.setProperties(atomContainer.getProperties());
            mol.setFlags(atomContainer.getFlags());
            if (atomContainer.getID() != null) {
                mol.setID(atomContainer.getID());
            }

        } else {
            mol = atomContainer.getBuilder().newMolecule(atomContainer);
            mol.setProperties(atomContainer.getProperties());
            mol.setFlags(atomContainer.getFlags());
            if (atomContainer.getID() != null) {
                mol.setID(atomContainer.getID());
            }
            if (atomContainer.getAtom(0).getSymbol().equalsIgnoreCase("H")) {
                System.err.println("WARNING: single hydrogen atom removal not supported!");
            }

        }

        return mol;
    }

    /**
     *
     * @param atomContainer
     * @note function added by Asad
     * @return IAtomContainer without Hydrogen. If an AtomContainer has a single atom which
     * is a Hydrogen then its not removed.
     * @throws Exception
     */
    public static IAtomContainer convertExplicitToImplicitHydrogens(IAtomContainer atomContainer) throws Exception {
        Map<IAtom, IAtom> map = new HashMap<IAtom, IAtom>();        // maps original atoms to clones.
        List<IAtom> remove = new ArrayList<IAtom>();  // lists removed Hs.
        IMolecule mol = null;
        if (atomContainer.getBondCount() > 0) {
            // Clone atoms except those to be removed.
            mol = atomContainer.getBuilder().newMolecule();
            int count = atomContainer.getAtomCount();
            for (int i = 0;
                    i < count;
                    i++) {
                // Clone/remove this atom?
                IAtom atom = atomContainer.getAtom(i);
                if (!atom.getSymbol().equals("H")) {
                    IAtom clonedAtom = null;
                    try {
                        clonedAtom = (IAtom) atom.clone();
                    } catch (CloneNotSupportedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                    //added by Asad to preserve the Atom ID for atom mapping without Hydrogen
                    clonedAtom.setID(atom.getID());
                    clonedAtom.setFlags(atom.getFlags());


                    mol.addAtom(clonedAtom);
                    map.put(atom, clonedAtom);
                } else {
                    remove.add(atom);   // maintain list of removed H.
                }
            }

            // Clone bonds except those involving removed atoms.
            count = atomContainer.getBondCount();
            for (int i = 0;
                    i < count;
                    i++) {
                // Check bond.
                final IBond bond = atomContainer.getBond(i);
                boolean removedBond = false;
                final int length = bond.getAtomCount();
                for (int k = 0;
                        k < length;
                        k++) {
                    if (remove.contains(bond.getAtom(k))) {
                        removedBond = true;
                        break;
                    }
                }

                // Clone/remove this bond?
                if (!removedBond) // if (!remove.contains(atoms[0]) && !remove.contains(atoms[1]))
                {
                    IBond clone = null;
                    try {
                        clone = (IBond) atomContainer.getBond(i).clone();
                    } catch (CloneNotSupportedException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                    assert clone != null;
                    clone.setAtoms(new IAtom[]{map.get(bond.getAtom(0)), map.get(bond.getAtom(1))});
                    mol.addBond(clone);
                }
            }



            // Recompute hydrogen counts of neighbours of removed Hydrogens.
            for (IAtom aRemove : remove) {
                // Process neighbours.
                for (IAtom iAtom : atomContainer.getConnectedAtomsList(aRemove)) {
                    final IAtom neighb = map.get(iAtom);
                    if (neighb == null) {
                        continue; // since for the case of H2, neight H has a heavy atom neighbor
                    }
                    //Added by Asad
                    if (!(neighb instanceof PseudoAtom)) {
                        neighb.setHydrogenCount(
                                (neighb.getHydrogenCount() == null ? 0 : neighb.getHydrogenCount()) + 1);
                    } else {
                        neighb.setHydrogenCount(0);
                    }
                }
            }
            mol.setProperties(atomContainer.getProperties());
            mol.setFlags(atomContainer.getFlags());
            if (atomContainer.getID() != null) {
                mol.setID(atomContainer.getID());
            }

        } else {
            mol = atomContainer.getBuilder().newMolecule(atomContainer);
            mol.setProperties(atomContainer.getProperties());
            mol.setFlags(atomContainer.getFlags());
            if (atomContainer.getID() != null) {
                mol.setID(atomContainer.getID());
            }
            if (atomContainer.getAtom(0).getSymbol().equalsIgnoreCase("H")) {
                System.err.println("WARNING: single hydrogen atom removal not supported!");
            }

        }

        return mol;
    }

    /**
     * Convenience method to perceive atom types for all <code>IAtom</code>s in the
     * <code>IAtomContainer</code>, using the <code>CDKAtomTypeMatcher</code>. If the
     * matcher finds a matching atom type, the <code>IAtom</code> will be configured
     * to have the same properties as the <code>IAtomType</code>. If no matching atom
     * type is found, no configuration is performed.
     * @see function added by Asad to fix the PseudoAtom configration
     * @param container
     * @throws CDKException
     */
    public static void percieveAtomTypesAndConfigureAtoms(IAtomContainer container) throws CDKException {
        CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(container.getBuilder());
        for (IAtom atom : container.atoms()) {
            if (!(atom instanceof PseudoAtom)) {

                IAtomType matched = matcher.findMatchingAtomType(container, atom);
                if (matched != null) {
                    AtomTypeManipulator.configure(atom, matched);
                }

            }
        }
    }
}



