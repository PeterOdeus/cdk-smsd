/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.core.canonical;

import java.util.TreeMap;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.PseudoAtom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomParity;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.ILonePair;
import org.openscience.cdk.interfaces.ISingleElectron;

/**
 *
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class ReOrganizeAtomContainer {

    public final static String INVARIANCE_PAIR = "InvariancePair";
    public final static String CANONICAL_LABEL = "CanonicalLable";

    /**
     *
     * @param container
     * @return
     */
    public synchronized static IAtomContainer makeDeepCopy(IAtomContainer container) {

        IAtomContainer newAtomContainer = DefaultChemObjectBuilder.getInstance().newAtomContainer();

        int lonePairCount = container.getLonePairCount();
        int singleElectronCount = container.getSingleElectronCount();
        ILonePair[] lonePairs = new ILonePair[lonePairCount];
        ISingleElectron[] singleElectrons = new ISingleElectron[singleElectronCount];

//        Deep copy of the Atoms

        TreeMap<Integer, Integer> noncanonicalSortedAtoms = new TreeMap<Integer, Integer>();
        TreeMap<Integer, Integer> canonicalSortedAtoms = new TreeMap<Integer, Integer>();

        for (Integer f = 0; f < container.getAtomCount(); f++) {
            IAtom a = container.getAtom(f);
            if (a.getProperty(CANONICAL_LABEL) != null) {
                Integer val = Integer.parseInt(a.getProperty(CANONICAL_LABEL).toString()) - 1;
                noncanonicalSortedAtoms.put(f, val);
                canonicalSortedAtoms.put(val, f);
            } else {
                return new Molecule(container);
            }


        }

        for (int f = 0; f < canonicalSortedAtoms.size(); f++) {

            IAtom oAtom = container.getAtom(canonicalSortedAtoms.get(f));

            IAtom nAtom = null;
            if (oAtom instanceof PseudoAtom) {

                nAtom = new PseudoAtom(oAtom);
            } else {
                nAtom = new Atom(oAtom);
            }

            if (oAtom.getID() != null) {

                nAtom.setID(oAtom.getID());
            }


            if ((oAtom).getPoint2d() != null) {
                nAtom.setPoint2d(new Point2d(oAtom.getPoint2d()));
            }
            if ((oAtom).getPoint3d() != null) {
                nAtom.setPoint3d(new Point3d(oAtom.getPoint3d()));
            }
            if ((oAtom).getFractionalPoint3d() != null) {
                nAtom.setFractionalPoint3d(new Point3d(oAtom.getFractionalPoint3d()));
            }

            if (oAtom.getID() != null) {
                nAtom.setID(new String(oAtom.getID()));
            }
            if (oAtom.getHydrogenCount() != null) {
                nAtom.setHydrogenCount(new Integer(oAtom.getHydrogenCount()));
            }
            if (oAtom.getCharge() != null) {
                nAtom.setCharge(new Double(oAtom.getCharge()));
            }
            if (oAtom.getStereoParity() != null) {
                nAtom.setStereoParity(new Integer(oAtom.getStereoParity()));
            }

            nAtom.setProperties(oAtom.getProperties());
            newAtomContainer.addAtom(nAtom);

            if (container.getAtomParity(oAtom) != null) {
                IAtomParity ap = container.getAtomParity(oAtom);
                newAtomContainer.addAtomParity(ap);
            }



        }

//  Deep copy of the bonds
        for (int iIndex = 0; iIndex < container.getAtomCount(); iIndex++) {

            IAtom a = container.getAtom(iIndex);

            for (int jIndex = 0; jIndex < container.getAtomCount(); jIndex++) {


                IAtom b = container.getAtom(jIndex);

                IBond oBond = null;
                if (container.getBond(a, b) != null) {

                    oBond = container.getBond(a, b);


                    IBond nBond = new Bond();



                    IAtom atom1 = newAtomContainer.getAtom(noncanonicalSortedAtoms.get(iIndex));
                    IAtom atom2 = newAtomContainer.getAtom(noncanonicalSortedAtoms.get(jIndex));

                    Order order = oBond.getOrder();
                    int stereo = oBond.getStereo();
                    nBond = new Bond(atom1, atom2, order, stereo);
                    if (oBond.getID() != null) {
                        nBond.setID(new String(oBond.getID()));
                    }

                    newAtomContainer.addBond(nBond);
                }

            }
        }

// Deep copy of the LonePairs

        for (int f = 0; f < container.getLonePairCount(); f++) {

            if (container.getAtom(f).getSymbol().equalsIgnoreCase("R")) {
                lonePairs[noncanonicalSortedAtoms.get(f).intValue()] = DefaultChemObjectBuilder.getInstance().newLonePair(container.getAtom(f));
            }
            newAtomContainer.addLonePair(lonePairs[noncanonicalSortedAtoms.get(f).intValue()]);

        }

        for (int f = 0; f < container.getSingleElectronCount(); f++) {
            singleElectrons[noncanonicalSortedAtoms.get(f).intValue()] = DefaultChemObjectBuilder.getInstance().newSingleElectron(container.getAtom(f));
            newAtomContainer.addSingleElectron(singleElectrons[noncanonicalSortedAtoms.get(f).intValue()]);

        }
        newAtomContainer.setProperties(container.getProperties());
        newAtomContainer.setFlags(container.getFlags());
        newAtomContainer.setID(container.getID());
        newAtomContainer.notifyChanged();
        return newAtomContainer;

    }
}
