/**
 *
 * Copyright (C) 2006-2009  Syed Asad Rahman {asad@ebi.atomContainer.uk}
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
 * You should have received atom copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.smsd.tools;

import java.util.ArrayList;
import java.util.HashMap;
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
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;

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



        int lonePairCount = container.getLonePairCount();
        int singleElectronCount = container.getSingleElectronCount();


        ILonePair[] lonePairs = new ILonePair[lonePairCount];
        ISingleElectron[] singleElectrons = new ISingleElectron[singleElectronCount];

//      Deep copy of the Atoms
        IAtom[] atoms = copyAtoms(container, newAtomContainer);

//      Deep copy of the bonds
        copyBonds(atoms, container, newAtomContainer);

//      Deep copy of the LonePairs
        for (int index = 0; index < container.getLonePairCount(); index++) {

            if (container.getAtom(index).getSymbol().equalsIgnoreCase("R")) {
                lonePairs[index] = DefaultChemObjectBuilder.getInstance().newLonePair(container.getAtom(index));
            }
            newAtomContainer.addLonePair(lonePairs[index]);
        }

        for (int index = 0; index < container.getSingleElectronCount(); index++) {
            singleElectrons[index] = DefaultChemObjectBuilder.getInstance().newSingleElectron(container.getAtom(index));
            newAtomContainer.addSingleElectron(singleElectrons[index]);

        }
        newAtomContainer.setProperties(container.getProperties());
        newAtomContainer.setFlags(container.getFlags());

        newAtomContainer.setID(container.getID());

        newAtomContainer.notifyChanged();
        return newAtomContainer;

    }

    /**
     *
     * @param mol Molecule to be aromatized
     */
    public static void aromatizeMolecule(IAtomContainer mol) {

        // need to find rings and aromaticity again since added H's

        IRingSet ringSet = null;
        try {
            AllRingsFinder arf = new AllRingsFinder();
            ringSet = arf.findAllRings(mol);

            // SSSRFinder s = new SSSRFinder(mol);
            // srs = s.findEssentialRings();

        } catch (Exception e) {
            e.printStackTrace();
        }

        try {
            // figure out which atoms are in aromatic rings:
//            printAtoms(mol);
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
//            printAtoms(mol);
            CDKHueckelAromaticityDetector.detectAromaticity(mol);
//            printAtoms(mol);
            // figure out which rings are aromatic:
            RingSetManipulator.markAromaticRings(ringSet);
//            printAtoms(mol);
            // figure out which simple (non cycles) rings are aromatic:
            // HueckelAromaticityDetector.detectAromaticity(mol, srs);
        } catch (Exception e) {
            e.printStackTrace();
        }


        // only atoms in 6 membered rings are aromatic
        // determine largest ring that each atom is atom part of

        for (int i = 0; i <= mol.getAtomCount() - 1; i++) {

            mol.getAtom(i).setFlag(CDKConstants.ISAROMATIC, false);

            jloop:
            for (int j = 0; j <= ringSet.getAtomContainerCount() - 1; j++) {
                //logger.debug(i+"\t"+j);
                IRing ring = (IRing) ringSet.getAtomContainer(j);
                if (!ring.getFlag(CDKConstants.ISAROMATIC)) {
                    continue jloop;
                }

                boolean haveatom = ring.contains(mol.getAtom(i));

                //logger.debug("haveatom="+haveatom);

                if (haveatom && ring.getAtomCount() == 6) {
                    mol.getAtom(i).setFlag(CDKConstants.ISAROMATIC, true);
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
     * @return IAtomContainer without Hydrogen. If an AtomContainer has atom single atom which
     * is atom Hydrogen then its not removed.
     */
    public static IAtomContainer removeHydrogensAndPreserveAtomID(IAtomContainer atomContainer) {
        Map<IAtom, IAtom> map = new HashMap<IAtom, IAtom>();        // maps original atoms to clones.
        List<IAtom> remove = new ArrayList<IAtom>();  // lists removed Hs.
        IMolecule mol = null;
        if (atomContainer.getBondCount() > 0) {
            // Clone atoms except those to be removed.
            mol = atomContainer.getBuilder().newMolecule();
            int count = atomContainer.getAtomCount();
            for (int i = 0; i < count; i++) {
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
//            Clone bonds except those involving removed atoms.
            mol = cloneAndMarkNonHBonds(mol, atomContainer, remove, map);
//            Recompute hydrogen counts of neighbours of removed Hydrogens.
            mol = reComputeHydrogens(mol, atomContainer, remove, map);

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
     * @return IAtomContainer without Hydrogen. If an AtomContainer has atom single atom which
     * is atom Hydrogen then its not removed.
     */
    public static IAtomContainer convertExplicitToImplicitHydrogens(IAtomContainer atomContainer) {
        Map<IAtom, IAtom> map = new HashMap<IAtom, IAtom>();        // maps original atoms to clones.
        List<IAtom> remove = new ArrayList<IAtom>();  // lists removed Hs.
        IMolecule mol = null;
        if (atomContainer.getBondCount() > 0) {
            // Clone atoms except those to be removed.
            mol = atomContainer.getBuilder().newMolecule();
            int count = atomContainer.getAtomCount();
            for (int i = 0; i < count; i++) {
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
            mol = cloneAndMarkNonHBonds(mol, atomContainer, remove, map);
            // Recompute hydrogen counts of neighbours of removed Hydrogens.
            mol = reComputeHydrogens(mol, atomContainer, remove, map);
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
     * matcher finds atom matching atom type, the <code>IAtom</code> will be configured
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

    private static IAtom[] copyAtoms(IAtomContainer container, IMolecule newAtomContainer) {
        int atomCount = container.getAtomCount();
        IAtom[] atoms = new IAtom[atomCount];
        for (int index = 0; index < container.getAtomCount(); index++) {

            if (container.getAtom(index) instanceof PseudoAtom) {
                atoms[index] = new PseudoAtom(container.getAtom(index));
            } else {
                atoms[index] = new Atom(container.getAtom(index));
            }

            set2D(container, index, atoms);
            set3D(container, index, atoms);
            setFractionalPoint3d(container, index, atoms);
            setID(container, index, atoms);
            setHydrogenCount(container, index, atoms);
            setCharge(container, index, atoms);
            setStereoParity(container, index, atoms);
            newAtomContainer.addAtom(atoms[index]);
            setAtomParity(container, index, newAtomContainer);

        }

        return atoms;
    }

    private static void copyBonds(IAtom[] atoms, IAtomContainer container, IMolecule newAtomContainer) {
        int bondCount = container.getBondCount();
        IBond[] bonds = new IBond[bondCount];
        for (int index = 0; index < container.getBondCount(); index++) {
            bonds[index] = new Bond();
            int IndexI = 999;
            for (int i = 0; i < container.getAtomCount(); i++) {
                if (container.getBond(index).getAtom(0) == container.getAtom(i)) {
                    IndexI = i;
                    break;
                }
            }
            int IndexJ = 999;
            for (int j = 0; j < container.getAtomCount(); j++) {
                if (container.getBond(index).getAtom(1) == container.getAtom(j)) {
                    IndexJ = j;
                    break;
                }
            }

            IAtom atom1 = atoms[IndexI];
            IAtom atom2 = atoms[IndexJ];

            Order order = container.getBond(index).getOrder();
            IBond.Stereo stereo = container.getBond(index).getStereo();
            bonds[index] = new Bond(atom1, atom2, order, stereo);
            if (container.getBond(index).getID() != null) {
                bonds[index].setID(new String(container.getBond(index).getID()));
            }
            newAtomContainer.addBond(bonds[index]);

        }
    }

    private static IMolecule reComputeHydrogens(
            IMolecule mol,
            IAtomContainer atomContainer,
            List<IAtom> remove,
            Map<IAtom, IAtom> map) {

        // Recompute hydrogen counts of neighbours of removed Hydrogens.
        for (IAtom aRemove : remove) {
            // Process neighbours.
            for (IAtom iAtom : atomContainer.getConnectedAtomsList(aRemove)) {
                final IAtom neighb = map.get(iAtom);
                if (neighb == null) {
                    continue; // since for the case of H2, neight H has atom heavy atom neighbor
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
        return mol;
    }

    private static IMolecule cloneAndMarkNonHBonds(
            IMolecule mol,
            IAtomContainer atomContainer,
            List<IAtom> remove,
            Map<IAtom, IAtom> map) {
        // Clone bonds except those involving removed atoms.
        int count = atomContainer.getBondCount();
        for (int i = 0; i < count; i++) {
            // Check bond.
            final IBond bond = atomContainer.getBond(i);
            boolean removedBond = false;
            final int length = bond.getAtomCount();
            for (int k = 0; k < length; k++) {
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

        return mol;
    }

    private static void set2D(IAtomContainer container, int index, IAtom[] atoms) {
        if ((container.getAtom(index)).getPoint2d() != null) {
            atoms[index].setPoint2d(new Point2d(container.getAtom(index).getPoint2d()));
        }
    }

    private static void set3D(IAtomContainer container, int index, IAtom[] atoms) {
        if ((container.getAtom(index)).getPoint3d() != null) {
            atoms[index].setPoint3d(new Point3d(container.getAtom(index).getPoint3d()));
        }
    }

    private static void setFractionalPoint3d(IAtomContainer container, int index, IAtom[] atoms) {
        if ((container.getAtom(index)).getFractionalPoint3d() != null) {
            atoms[index].setFractionalPoint3d(new Point3d(container.getAtom(index).getFractionalPoint3d()));
        }
    }

    private static void setID(IAtomContainer container, int index, IAtom[] atoms) {

        if (container.getAtom(index).getID() != null) {
            atoms[index].setID(new String(container.getAtom(index).getID()));
        }
    }

    private static void setHydrogenCount(IAtomContainer container, int index, IAtom[] atoms) {
        if (container.getAtom(index).getHydrogenCount() != null) {
            atoms[index].setHydrogenCount(Integer.valueOf(container.getAtom(index).getHydrogenCount()));
        }
    }

    private static void setCharge(IAtomContainer container, int index, IAtom[] atoms) {
        if (container.getAtom(index).getCharge() != null) {
            atoms[index].setCharge(new Double(container.getAtom(index).getCharge()));
        }
    }

    private static void setStereoParity(IAtomContainer container, int index, IAtom[] atoms) {
        if (container.getAtom(index).getStereoParity() != null) {
            atoms[index].setStereoParity(Integer.valueOf(container.getAtom(index).getStereoParity()));
        }
    }

    private static void setAtomParity(IAtomContainer container, int index, IAtomContainer newAtomContainer) {
        if (container.getAtomParity(container.getAtom(index)) != null) {
            IAtomParity atomParity = container.getAtomParity(container.getAtom(index));
            newAtomContainer.addAtomParity(atomParity);
        }
    }
}



