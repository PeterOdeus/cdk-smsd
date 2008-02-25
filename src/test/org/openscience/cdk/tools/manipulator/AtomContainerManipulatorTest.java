/* $RCSfile$    
 * $Author$    
 * $Date$    
 * $Revision$
 * 
 * Copyright (C) 1997-2007  The Chemistry Development Kit (CDK) project
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
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
package org.openscience.cdk.tools.manipulator;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.*;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.templates.MoleculeFactory;
import org.openscience.cdk.NewCDKTestCase;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.IOException;
import java.util.List;

/**
 * @cdk.module test-standard
 */
public class AtomContainerManipulatorTest extends NewCDKTestCase {
    IAtomContainer ac;

    public AtomContainerManipulatorTest()
    {
        super();
    }

    @Before
    public void setUp()
    {
        ac = MoleculeFactory.makeAlphaPinene();
    }


    @Test
    public void testGetTotalHydrogenCount_IAtomContainer() throws IOException, ClassNotFoundException, CDKException {
        Molecule mol = new Molecule(); // ethene
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addBond(0, 1, IBond.Order.DOUBLE);
        mol.addBond(0, 2, IBond.Order.SINGLE);
        mol.addBond(0, 3, IBond.Order.SINGLE);
        mol.addBond(1, 4, IBond.Order.SINGLE);
        mol.addBond(1, 5, IBond.Order.SINGLE);
        Assert.assertEquals(6, mol.getAtomCount());
        Assert.assertEquals(5, mol.getBondCount());
        Assert.assertEquals(0, AtomContainerManipulator.getTotalHydrogenCount(mol));
    }
        
    @Test public void testConvertImplicitToExplicitHydrogens_IAtomContainer() throws Exception {
        Molecule mol = new Molecule(); // ethene
        mol.addAtom(new Atom("C"));
        mol.getAtom(0).setHydrogenCount(2);
        mol.addAtom(new Atom("C"));
        mol.getAtom(1).setHydrogenCount(2);
        mol.addBond(0,1, CDKConstants.BONDORDER_DOUBLE);
        Assert.assertEquals(2, mol.getAtomCount());
        Assert.assertEquals(1, mol.getBondCount());
        
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
        Assert.assertEquals(6, mol.getAtomCount());
        Assert.assertEquals(5, mol.getBondCount());
    }
        
    @Test public void testGetTotalHydrogenCount_IAtomContainer_zeroImplicit() throws IOException, ClassNotFoundException, CDKException {
        Molecule mol = new Molecule(); // ethene
        mol.addAtom(new Atom("C"));
        mol.getAtom(0).setHydrogenCount(0);
        mol.addAtom(new Atom("C"));
        mol.getAtom(1).setHydrogenCount(0);
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addBond(0, 1, IBond.Order.DOUBLE);
        mol.addBond(0, 2, IBond.Order.SINGLE);
        mol.addBond(0, 3, IBond.Order.SINGLE);
        mol.addBond(1, 4, IBond.Order.SINGLE);
        mol.addBond(1, 5, IBond.Order.SINGLE);
        Assert.assertEquals(6, mol.getAtomCount());
        Assert.assertEquals(5, mol.getBondCount());

        AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
        Assert.assertEquals(6, mol.getAtomCount());
        Assert.assertEquals(5, mol.getBondCount());
    }

    @Test public void testGetTotalHydrogenCount_IAtomContainer_nullImplicit() throws IOException, ClassNotFoundException, CDKException {
        Molecule mol = new Molecule(); // ethene
        mol.addAtom(new Atom("C"));
        mol.getAtom(0).setHydrogenCount(null);
        mol.addAtom(new Atom("C"));
        mol.getAtom(1).setHydrogenCount(null);
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addBond(0, 1, IBond.Order.DOUBLE);
        mol.addBond(0, 2, IBond.Order.SINGLE);
        mol.addBond(0, 3, IBond.Order.SINGLE);
        mol.addBond(1, 4, IBond.Order.SINGLE);
        mol.addBond(1, 5, IBond.Order.SINGLE);
        Assert.assertEquals(6, mol.getAtomCount());
        Assert.assertEquals(5, mol.getBondCount());

        AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
        Assert.assertEquals(6, mol.getAtomCount());
        Assert.assertEquals(5, mol.getBondCount());
    }

    @Test public void testGetTotalHydrogenCount_ImplicitHydrogens() throws Exception {
        Molecule mol = new Molecule();
        Atom carbon = new Atom("C");
        carbon.setHydrogenCount(4);
        mol.addAtom(carbon);
        Assert.assertEquals(4, AtomContainerManipulator.getTotalHydrogenCount(mol));
    }
        
    @Test public void testRemoveHydrogens_IAtomContainer() throws IOException, ClassNotFoundException, CDKException{
        Molecule mol = new Molecule(); // ethene
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addBond(0, 1, IBond.Order.DOUBLE);
        mol.addBond(0, 2, IBond.Order.SINGLE);
        mol.addBond(0, 3, IBond.Order.SINGLE);
        mol.addBond(1, 4, IBond.Order.DOUBLE);
        mol.addBond(1, 5, IBond.Order.DOUBLE);
        mol.setFlag(5,true);
        
        Assert.assertEquals(6, mol.getAtomCount());
        IAtomContainer ac = AtomContainerManipulator.removeHydrogens(mol);
        Assert.assertEquals(2, ac.getAtomCount());
        Assert.assertTrue(ac.getFlag(5));
    }

    @Test public void testGetAllIDs_IAtomContainer() {
        Molecule mol = new Molecule(); // ethene
        mol.addAtom(new Atom("C")); mol.getAtom(0).setID("a1");
        mol.addAtom(new Atom("C")); mol.getAtom(1).setID("a2");
        mol.addAtom(new Atom("H")); mol.getAtom(2).setID("a3");
        mol.addAtom(new Atom("H")); mol.getAtom(3).setID("a4");
        mol.addAtom(new Atom("H")); mol.getAtom(4).setID("a5");
        mol.addAtom(new Atom("H")); mol.getAtom(5).setID("a6");
        
        List ids = AtomContainerManipulator.getAllIDs(mol);
        Assert.assertEquals(6, ids.size());
        Assert.assertTrue(ids.contains("a1"));
        Assert.assertTrue(ids.contains("a2"));
        Assert.assertTrue(ids.contains("a3"));
        Assert.assertTrue(ids.contains("a4"));
        Assert.assertTrue(ids.contains("a5"));
        Assert.assertTrue(ids.contains("a6"));
    }

    @Test public void testGetAtomArray_IAtomContainer() {
        Molecule mol = new Molecule(); // ethene
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        
        IAtom[] atoms = AtomContainerManipulator.getAtomArray(mol);
        Assert.assertEquals(6, atoms.length);
        Assert.assertEquals(mol.getAtom(0), atoms[0]);
        Assert.assertEquals(mol.getAtom(1), atoms[1]);
        Assert.assertEquals(mol.getAtom(2), atoms[2]);
        Assert.assertEquals(mol.getAtom(3), atoms[3]);
        Assert.assertEquals(mol.getAtom(4), atoms[4]);
        Assert.assertEquals(mol.getAtom(5), atoms[5]);
    }

    @Test public void testGetAtomArray_List() {
        Molecule mol = new Molecule(); // ethene
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addBond(0, 1, IBond.Order.DOUBLE);
        mol.addBond(0, 2, IBond.Order.SINGLE);
        mol.addBond(0, 3, IBond.Order.SINGLE);
        mol.addBond(1, 4, IBond.Order.DOUBLE);
        mol.addBond(1, 5, IBond.Order.DOUBLE);
        mol.setFlag(5,true);
        
        IAtom[] atoms = AtomContainerManipulator.getAtomArray(
        	mol.getConnectedAtomsList(mol.getAtom(0))
        );
        Assert.assertEquals(3, atoms.length);
        Assert.assertEquals(mol.getAtom(1), atoms[0]);
        Assert.assertEquals(mol.getAtom(2), atoms[1]);
        Assert.assertEquals(mol.getAtom(3), atoms[2]);
    }

    @Test public void testGetBondArray_List() {
        Molecule mol = new Molecule(); // ethene
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addBond(0, 1, IBond.Order.DOUBLE);
        mol.addBond(0, 2, IBond.Order.SINGLE);
        mol.addBond(0, 3, IBond.Order.SINGLE);
        mol.addBond(1, 4, IBond.Order.DOUBLE);
        mol.addBond(1, 5, IBond.Order.DOUBLE);
        mol.setFlag(5,true);
        
        IBond[] bonds = AtomContainerManipulator.getBondArray(
        	mol.getConnectedBondsList(mol.getAtom(0))
        );
        Assert.assertEquals(3, bonds.length);
        Assert.assertEquals(mol.getBond(0), bonds[0]);
        Assert.assertEquals(mol.getBond(1), bonds[1]);
        Assert.assertEquals(mol.getBond(2), bonds[2]);
    }

    @Test public void testGetBondArray_IAtomContainer() {
        Molecule mol = new Molecule(); // ethene
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("C"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addAtom(new Atom("H"));
        mol.addBond(0, 1, IBond.Order.DOUBLE);
        mol.addBond(0, 2, IBond.Order.SINGLE);
        mol.addBond(0, 3, IBond.Order.SINGLE);
        mol.addBond(1, 4, IBond.Order.DOUBLE);
        mol.addBond(1, 5, IBond.Order.DOUBLE);
        mol.setFlag(5,true);
        
        IBond[] bonds = AtomContainerManipulator.getBondArray(mol);
        Assert.assertEquals(5, bonds.length);
        Assert.assertEquals(mol.getBond(0), bonds[0]);
        Assert.assertEquals(mol.getBond(1), bonds[1]);
        Assert.assertEquals(mol.getBond(2), bonds[2]);
        Assert.assertEquals(mol.getBond(3), bonds[3]);
        Assert.assertEquals(mol.getBond(4), bonds[4]);
    }

    @Test public void testGetAtomById_IAtomContainer_String() throws CDKException {
        Molecule mol = new Molecule(); // ethene
        mol.addAtom(new Atom("C")); mol.getAtom(0).setID("a1");
        mol.addAtom(new Atom("C")); mol.getAtom(1).setID("a2");
        mol.addAtom(new Atom("H")); mol.getAtom(2).setID("a3");
        mol.addAtom(new Atom("H")); mol.getAtom(3).setID("a4");
        mol.addAtom(new Atom("H")); mol.getAtom(4).setID("a5");
        mol.addAtom(new Atom("H")); mol.getAtom(5).setID("a6");
        mol.addBond(0, 1, IBond.Order.DOUBLE);
        mol.addBond(0, 2, IBond.Order.SINGLE);
        mol.addBond(0, 3, IBond.Order.SINGLE);
        mol.addBond(1, 4, IBond.Order.DOUBLE);
        mol.addBond(1, 5, IBond.Order.DOUBLE);
        mol.setFlag(5,true);
        
        Assert.assertEquals(mol.getAtom(0), AtomContainerManipulator.getAtomById(mol, "a1"));
        Assert.assertEquals(mol.getAtom(1), AtomContainerManipulator.getAtomById(mol, "a2"));
        Assert.assertEquals(mol.getAtom(2), AtomContainerManipulator.getAtomById(mol, "a3"));
        Assert.assertEquals(mol.getAtom(3), AtomContainerManipulator.getAtomById(mol, "a4"));
        Assert.assertEquals(mol.getAtom(4), AtomContainerManipulator.getAtomById(mol, "a5"));
        Assert.assertEquals(mol.getAtom(5), AtomContainerManipulator.getAtomById(mol, "a6"));
    }

    /**
     * Test removeHydrogens for B2H6, which contains two multiply bonded H.
     *
     * @throws IOException
     * @throws ClassNotFoundException
     * @throws CDKException
     */
    @Test public void testRemoveHydrogensBorane() throws IOException, ClassNotFoundException, CDKException
    {
    	IAtomContainer borane = new Molecule();
    	borane.addAtom(borane.getBuilder().newAtom("H"));
    	borane.addAtom(borane.getBuilder().newAtom("H"));
    	borane.addAtom(borane.getBuilder().newAtom("B"));
    	borane.addAtom(borane.getBuilder().newAtom("H"));
    	borane.addAtom(borane.getBuilder().newAtom("H"));
    	borane.addAtom(borane.getBuilder().newAtom("B"));
    	borane.addAtom(borane.getBuilder().newAtom("H"));
    	borane.addAtom(borane.getBuilder().newAtom("H"));
    	borane.addBond(0,2,CDKConstants.BONDORDER_SINGLE);
    	borane.addBond(1,2,CDKConstants.BONDORDER_SINGLE);
    	borane.addBond(2,3,CDKConstants.BONDORDER_SINGLE); // REALLY 3-CENTER-2-ELECTRON
    	borane.addBond(2,4,CDKConstants.BONDORDER_SINGLE); // REALLY 3-CENTER-2-ELECTRON
    	borane.addBond(3,5,CDKConstants.BONDORDER_SINGLE); // REALLY 3-CENTER-2-ELECTRON
    	borane.addBond(4,5,CDKConstants.BONDORDER_SINGLE); // REALLY 3-CENTER-2-ELECTRON
    	borane.addBond(5,6,CDKConstants.BONDORDER_SINGLE);
    	borane.addBond(5,7,CDKConstants.BONDORDER_SINGLE);
        IAtomContainer ac = AtomContainerManipulator.removeHydrogens(borane);

        // Should be two disconnected Bs with H-count == 4
        Assert.assertEquals("incorrect atom count", 2, ac.getAtomCount());
        Assert.assertEquals("incorrect bond count", 0, ac.getBondCount());
        Assert.assertEquals("incorrect hydrogen count", 4, ac.getAtom(0).getHydrogenCount());
        Assert.assertEquals("incorrect hydrogen count", 4, ac.getAtom(1).getHydrogenCount());
    }
    /**
     * Test total formal charge.
     *
     * @throws IOException
     * @throws ClassNotFoundException
     * @throws CDKException
     */
    @Test public void testGetTotalFormalCharge_IAtomContainer() throws IOException, ClassNotFoundException, CDKException
    {
        SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IMolecule mol = parser.parseSmiles("[O-]C([N+])C([N+])C");
        int totalCharge = AtomContainerManipulator.getTotalFormalCharge(mol);

        Assert.assertEquals(1,totalCharge);
    }
    
    /**
     * Test total Exact Mass.
     *
     * @throws IOException
     * @throws ClassNotFoundException
     * @throws CDKException
     */
    @Test public void testGetTotalExactMass_IAtomContainer() throws IOException, ClassNotFoundException, CDKException{
    	
        SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IMolecule mol = parser.parseSmiles("CCl");
        mol.getAtom(0).setExactMass(12.00);
        mol.getAtom(1).setExactMass(34.96885268);
        double totalExactMass = AtomContainerManipulator.getTotalExactMass(mol);

        Assert.assertEquals(46.96885268,totalExactMass,0.000001);
    }
    
    /**
     * Test total natural abundance.
     *
     * @throws IOException
     * @throws ClassNotFoundException
     * @throws CDKException
     */
    @Test public void testGetTotalNaturalAbundance_IAtomContainer() throws IOException, ClassNotFoundException, CDKException{
    	
        SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IMolecule mol = parser.parseSmiles("CCl");
        mol.getAtom(0).setNaturalAbundance(98.93);
        mol.getAtom(1).setNaturalAbundance(75.78);
        double totalAbudance = AtomContainerManipulator.getTotalNaturalAbundance(mol);

        Assert.assertEquals(0.74969154,totalAbudance,0.000001);
    }
    
    /**
     * Test total positive formal charge.
     *
     * @throws IOException
     * @throws ClassNotFoundException
     * @throws CDKException
     */
    @Test public void testGetTotalPositiveFormalCharge_IAtomContainer() throws IOException, ClassNotFoundException, CDKException
    {
        SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IMolecule mol = parser.parseSmiles("[O-]C([N+])C([N+])C");
        int totalCharge = AtomContainerManipulator.getTotalPositiveFormalCharge(mol);

        Assert.assertEquals(2,totalCharge);
    }
    
    /**
     * Test total negative formal charge.
     *
     * @throws IOException
     * @throws ClassNotFoundException
     * @throws CDKException
     */
    @Test public void testGetTotalNegativeFormalCharge_IAtomContainer() throws IOException, ClassNotFoundException, CDKException
    {
        SmilesParser parser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IMolecule mol = parser.parseSmiles("[O-]C([N+])C([N+])C");
        int totalCharge = AtomContainerManipulator.getTotalNegativeFormalCharge(mol);

        Assert.assertEquals(-1,totalCharge);
    }

    @Test public void testGetIntersection_IAtomContainer_IAtomContainer() {
    	IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IAtom c1 = builder.newAtom("C");
        IAtom o = builder.newAtom("O");
        IAtom c2 = builder.newAtom("C");
        IAtom c3 = builder.newAtom("C");
        
        IBond b1 = builder.newBond(c1, o);
        IBond b2 = builder.newBond(o, c2);
        IBond b3 = builder.newBond(c2, c3);
        
        IAtomContainer container1 = new org.openscience.cdk.AtomContainer();
        container1.addAtom(c1);
        container1.addAtom(o);
        container1.addAtom(c2);
        container1.addBond(b1);
        container1.addBond(b2);
        IAtomContainer container2 = new org.openscience.cdk.AtomContainer();
        container2.addAtom(o);
        container2.addAtom(c3);
        container2.addAtom(c2);
        container2.addBond(b3);
        container2.addBond(b2);

        IAtomContainer intersection = AtomContainerManipulator.getIntersection(container1, container2);
        Assert.assertEquals(2, intersection.getAtomCount());
        Assert.assertEquals(1, intersection.getBondCount());
        Assert.assertTrue(intersection.contains(b2));
        Assert.assertTrue(intersection.contains(o));
        Assert.assertTrue(intersection.contains(c2));
    }
 
    @Test
    public void testPerceiveAtomTypesAndConfigureAtoms() {
    	IAtomContainer container = new AtomContainer();
        container.addAtom(new Atom("R"));
        
        // the next should not throw an exception
        try {
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
		} catch (CDKException e) {
			Assert.fail("The percieveAtomTypesAndConfigureAtoms must not throw exceptions when no atom type is perceived.");
		}
    }

    @Test
    public void testGetMaxBondOrder() {
        Assert.assertEquals(CDKConstants.BONDORDER_DOUBLE, AtomContainerManipulator.getMaximumBondOrder(ac));
    }

    @Test
    public void testGetSBE() {
        Assert.assertEquals(12, AtomContainerManipulator.getSingleBondEquivalentSum(ac));
    }

    @Test
    public void testGetTotalCharge() throws IOException, ClassNotFoundException, CDKException {
        IAtomContainer container = DefaultChemObjectBuilder.getInstance().newAtomContainer();
        IAtom atom1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        atom1.setCharge(1.0);
        IAtom atom2 = DefaultChemObjectBuilder.getInstance().newAtom("N");

        container.addAtom(atom1);
        container.addAtom(atom2);
        container.addBond(new Bond(atom1, atom2, CDKConstants.BONDORDER_SINGLE));

        double totalCharge = AtomContainerManipulator.getTotalCharge(container);

        Assert.assertEquals(1.0, totalCharge);
    }

    @Test
    public void testCountExplicitH() {
        IAtomContainer container = DefaultChemObjectBuilder.getInstance().newAtomContainer();
        IAtom atom1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        atom1.setCharge(1.0);
        IAtom atom2 = DefaultChemObjectBuilder.getInstance().newAtom("N");

        container.addAtom(atom1);
        container.addAtom(atom2);
        container.addBond(new Bond(atom1, atom2, CDKConstants.BONDORDER_SINGLE));

        Assert.assertEquals(0, AtomContainerManipulator.countExplicitHydrogens(container, atom1));
        Assert.assertEquals(0, AtomContainerManipulator.countExplicitHydrogens(container, atom2));

        for (int i = 0; i < 3; i++) {
            IAtom h = DefaultChemObjectBuilder.getInstance().newAtom("H");
            container.addAtom(h);
            container.addBond(new Bond(atom1, h, CDKConstants.BONDORDER_SINGLE));
        }
        Assert.assertEquals(3, AtomContainerManipulator.countExplicitHydrogens(container, atom1));                
    }

    @Test
    public void testCountH() throws CDKException {
        IAtomContainer container = DefaultChemObjectBuilder.getInstance().newAtomContainer();
        IAtom atom1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        atom1.setCharge(1.0);
        IAtom atom2 = DefaultChemObjectBuilder.getInstance().newAtom("N");

        container.addAtom(atom1);
        container.addAtom(atom2);
        container.addBond(new Bond(atom1, atom2, CDKConstants.BONDORDER_SINGLE));

        // no atom type perception, so implicit count is 0
        Assert.assertEquals(0, AtomContainerManipulator.countHydrogens(container, atom1));
        Assert.assertEquals(0, AtomContainerManipulator.countHydrogens(container, atom2));

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
        CDKHydrogenAdder ha = CDKHydrogenAdder.getInstance(DefaultChemObjectBuilder.getInstance());
        ha.addImplicitHydrogens(container);

        Assert.assertEquals(3, AtomContainerManipulator.countHydrogens(container, atom1));
        Assert.assertEquals(2, AtomContainerManipulator.countHydrogens(container, atom2));

        
        for (int i = 0; i < 3; i++) {
            IAtom h = DefaultChemObjectBuilder.getInstance().newAtom("H");
            container.addAtom(h);
            container.addBond(new Bond(atom1, h, CDKConstants.BONDORDER_SINGLE));
        }
        Assert.assertEquals(6, AtomContainerManipulator.countHydrogens(container, atom1));
        
    }

    @Test
    public void testReplaceAtom() {
        IAtomContainer container = DefaultChemObjectBuilder.getInstance().newAtomContainer();
        IAtom atom1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        atom1.setCharge(1.0);
        IAtom atom2 = DefaultChemObjectBuilder.getInstance().newAtom("N");

        container.addAtom(atom1);
        container.addAtom(atom2);
        container.addBond(new Bond(atom1, atom2, CDKConstants.BONDORDER_SINGLE));

        IAtom atom3 = DefaultChemObjectBuilder.getInstance().newAtom("Br");

        AtomContainerManipulator.replaceAtomByAtom(container, atom2, atom3);
        Assert.assertEquals(atom3, container.getAtom(1));
    }
}

