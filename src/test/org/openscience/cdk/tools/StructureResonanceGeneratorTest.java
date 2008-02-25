/* $Revision: 6221 $ $Author: miguelrojasch $ $Date: 2006-05-11 14:25:07 +0200 (Do, 11 Mai 2006) $
 *
 * Copyright (C) 2004-2007  Miguel Rojas <miguel.rojas@uni-koeln.de>
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
package org.openscience.cdk.tools;

import java.util.List;

import junit.framework.Assert;
import junit.framework.Test;
import junit.framework.TestSuite;

import org.openscience.cdk.Atom;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.LonePair;
import org.openscience.cdk.SingleElectron;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.ILonePair;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.tools.LonePairElectronChecker;
import org.openscience.cdk.tools.StructureResonanceGenerator;

/**
* TestSuite that runs all QSAR tests.
*
* @cdk.module test-reaction
*/
public class StructureResonanceGeneratorTest  extends CDKTestCase
{

	private StructureResonanceGenerator gR;

	/**
	 *  Constructor for the StructureResonanceGeneratorTest object
	 *
	 *@param  name  Description of the Parameter
	 */
	public StructureResonanceGeneratorTest(String name)
	{
		super(name);
	}

    /**
    *  The JUnit setup method
    */
    public void setUp() throws Exception {
    	gR = new StructureResonanceGenerator();
    }

	/**
	 * A unit test suite for JUnit
	 *
	 * @return    The test suite
	 */
    public static Test suite() {
        TestSuite suite = new TestSuite(StructureResonanceGeneratorTest.class);
        return suite;
	}
    
    /**
	 * <p>A unit test suite for JUnit: Resonance - CC(=[O*+])C=O</p>
	 * <p>CC(=[O*+])C=O <=> C[C+]([O*])C=O <=> CC([O*])=CO <=> CC(=O)[C*][O+] <=> CC(=O)C=[O*+]</p>
	 *
	 * @return    The test suite
	 */
	public void testGetAllStructures() throws ClassNotFoundException, CDKException, java.lang.Exception {
		IMolecule molecule = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("CC(=O)C=O");
        addExplicitHydrogens(molecule);
        LonePairElectronChecker lpcheck = new LonePairElectronChecker();
        lpcheck.saturate(molecule);
		
        IAtom atom =  molecule.getAtom(2);
        molecule.addSingleElectron(new SingleElectron(atom));
        atom.setFormalCharge(1);
        List selectron = molecule.getConnectedLonePairsList(atom);
		molecule.removeLonePair((ILonePair)selectron.get(0));

		StructureResonanceGenerator gRI = new StructureResonanceGenerator(true,true,true,true,false,false,-1);
		IAtomContainerSet setOfMolecules = gRI.getAllStructures(molecule);

		Assert.assertEquals(8,setOfMolecules.getAtomContainerCount());
		
		/*1*/
        IMolecule molecule1 = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("C[C+](O)C=O");
        for(int i = 0; i < 4; i++)
			molecule1.addAtom(new Atom("H"));
		molecule1.addBond(0, 5, IBond.Order.SINGLE);
	    molecule1.addBond(0, 6, IBond.Order.SINGLE);
	    molecule1.addBond(0, 7, IBond.Order.SINGLE);
	    molecule1.addBond(3, 8, IBond.Order.SINGLE);
        lpcheck.saturate(molecule1);
        IAtom atom1 =  molecule1.getAtom(2);
        molecule1.addSingleElectron(new SingleElectron(atom1));
        QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule1);
		Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(1),qAC));
		
//		/*2*/
//		Molecule molecule2 = (new SmilesParser()).parseSmiles("CC(O)=CO");
//		for(int i = 0; i < 4; i++)
//			molecule2.addAtom(new Atom("H"));
//		molecule2.addBond(0, 5, IBond.Order.SINGLE);
//	    molecule2.addBond(0, 6, IBond.Order.SINGLE);
//	    molecule2.addBond(0, 7, IBond.Order.SINGLE);
//	    molecule2.addBond(3, 8, IBond.Order.SINGLE);
//        lpcheck.newSaturate(molecule2);
//		IAtom atom2a =  molecule2.getAtom(2);
//		molecule2.addElectronContainer(new SingleElectron(atom2a));
//
//		IAtom atom2b =  molecule2.getAtom(4);
//		atom2b.setHydrogenCount(0);
//		atom2b.setFormalCharge(1);
//		
//		qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule2);
//		Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(3),qAC));
	}
	/**
	 * A unit test suite for JUnit: Resonance CC(=[O*+])C=O <=> CC(=O)C=[O*+]
	 *
	 * @return    The test suite
	 */
	public void testGetStructures1() throws ClassNotFoundException, CDKException, java.lang.Exception {
		IMolecule molecule = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("CC(=O)C=O");
        addExplicitHydrogens(molecule);
        LonePairElectronChecker lpcheck = new LonePairElectronChecker();
        lpcheck.saturate(molecule);
		
        IAtom atom =  molecule.getAtom(2);
        molecule.addSingleElectron(new SingleElectron(atom));
        atom.setFormalCharge(1);
        List selectron = molecule.getConnectedLonePairsList(atom);
		molecule.removeLonePair((ILonePair)selectron.get(selectron.size()-1));

		StructureResonanceGenerator gRI = new StructureResonanceGenerator();
		IAtomContainerSet setOfMolecules = gRI.getStructures(molecule);

		Assert.assertEquals(2,setOfMolecules.getAtomContainerCount());
		
		IMolecule molecule1 = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("CC(=O)C=O");
		addExplicitHydrogens(molecule1);
		lpcheck.saturate(molecule1);
		IAtom atom1 =  molecule1.getAtom(4);
		molecule1.addSingleElectron(new SingleElectron(atom1));	
		selectron = molecule1.getConnectedLonePairsList(atom1);
		molecule1.removeLonePair((ILonePair)selectron.get(0));
		atom1.setFormalCharge(1);
		

		QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule1);
		Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(1),qAC));
	
	}	
	/**
	 * A unit test suite for JUnit: Resonance CCC(=[O*+])C(C)=O <=> CCC(=O)C(C)=[O*+]
	 *
	 * @return    The test suite
	 */
	public void testGetStructures2() throws ClassNotFoundException, CDKException, java.lang.Exception {
		IMolecule molecule = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("CCC(=O)C(C)=O");
        addExplicitHydrogens(molecule);
        LonePairElectronChecker lpcheck = new LonePairElectronChecker();
        lpcheck.saturate(molecule);
        
        IAtom atom =  molecule.getAtom(3);
        molecule.addSingleElectron(new SingleElectron(atom));
        atom.setFormalCharge(1);
        List selectron = molecule.getConnectedLonePairsList(atom);
		molecule.removeLonePair((ILonePair)selectron.get(0));

		StructureResonanceGenerator gRI = new StructureResonanceGenerator();
		IAtomContainerSet setOfMolecules = gRI.getStructures(molecule);

		Assert.assertEquals(2,setOfMolecules.getAtomContainerCount());
		
		IMolecule molecule1 = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("CCC(=O)C(C)=O");
		addExplicitHydrogens(molecule1);
		lpcheck.saturate(molecule1);

        IAtom atom1 =  molecule1.getAtom(6);
        molecule1.addSingleElectron(new SingleElectron(atom1));
        atom1.setFormalCharge(1);
        selectron = molecule.getConnectedLonePairsList(atom);
		molecule.removeLonePair((ILonePair)selectron.get(0));

		QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule1);
		Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(1),qAC));
	
	}
	/**
	 * A unit test suite for JUnit: Resonance C-C=C-[C+]-C-C=C-[C+] <=> C-[C+]-C=C-C-C=C-[C+]
	 *
	 * @return    The test suite
	 */
	public void testFlagActiveCenter1() throws ClassNotFoundException, CDKException, java.lang.Exception {
		IMolecule molecule = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("CC=C[C+]-C-C=C-C=C");
        
        molecule.getAtom(1).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getBond(1).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getAtom(2).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getBond(2).setFlag(CDKConstants.REACTIVE_CENTER,true);
		molecule.getAtom(3).setFlag(CDKConstants.REACTIVE_CENTER,true);
		
        StructureResonanceGenerator gRI = new StructureResonanceGenerator();
		IAtomContainerSet setOfMolecules = gRI.getStructures(molecule);
        
		Assert.assertEquals(2,setOfMolecules.getAtomContainerCount());

        IMolecule molecule2 = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("C-[C+]-C=C-C-C=C-C=C");
        
        QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule2);
        Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(1),qAC));
	}
	/**
	 * A unit test suite for JUnit: Resonance C-C=C-[C-] <=> C=C-[C-]-C
	 *
	 * @return    The test suite
	 */
	public void testFlagActiveCenter2() throws ClassNotFoundException, CDKException, java.lang.Exception {
		IMolecule molecule = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("C-C=C-[C-]");
		molecule.addLonePair(new LonePair(molecule.getAtom(3)));
		
        IAtomContainerSet setOfMolecules = gR.getStructures(molecule);
        
		Assert.assertEquals(2,setOfMolecules.getAtomContainerCount());

        IMolecule molecule2 = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("C=C-[C-]-C");
        molecule.addLonePair(new LonePair(molecule.getAtom(2)));
		
        QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule2);
        Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(1),qAC));
	}
	/**
	 * A unit test suite for JUnit: Resonance Formic acid  C(=O)O <=> [C-](-[O+])O <=> [C+](-[O-])O <=> C([O-])=[O+]
	 *
	 * @return    The test suite
	 */
	public void testResonanceFormicAcid() throws ClassNotFoundException, CDKException, java.lang.Exception {
		IMolecule molecule = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("C(=O)O");
        addImplicitHydrogens(molecule);
        LonePairElectronChecker lpcheck = new LonePairElectronChecker();
        lpcheck.saturate(molecule);
        
        IAtomContainerSet setOfMolecules = gR.getAllStructures(molecule);
        
		Assert.assertEquals(3,setOfMolecules.getAtomContainerCount());

//        Molecule molecule1 = (new SmilesParser()).parseSmiles("[C-](-[O+])O");
//        addImplicitHydrogens(molecule1);
//        lpcheck.newSaturate(molecule1);
//        QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule1);
//        Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(1),qAC));

		IMolecule molecule2 = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("[C+](-[O-])O");
        addImplicitHydrogens(molecule2);
        lpcheck.saturate(molecule2);
        QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule2);
        Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(1),qAC));
        
        IMolecule molecule3 = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("C([O-])=[O+]");
        addImplicitHydrogens(molecule3);
        lpcheck.saturate(molecule3);
        
        qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule3);
        Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(2),qAC));
	}

	/**
	 * A unit test suite for JUnit: Resonance Formic acid  F-C=C <=> [F+]=C-[C-]
	 *
	 * @return    The test suite
	 */
	public void testResonanceFluoroethylene() throws ClassNotFoundException, CDKException, java.lang.Exception {
		IMolecule molecule = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("F-C=C");
        addImplicitHydrogens(molecule);
        LonePairElectronChecker lpcheck = new LonePairElectronChecker();
        lpcheck.saturate(molecule);
        
        IAtomContainerSet setOfMolecules = gR.getAllStructures(molecule);
        
		Assert.assertEquals(2,setOfMolecules.getAtomContainerCount());

        IMolecule molecule1 = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("[F+]=C-[C-]");
        addImplicitHydrogens(molecule1);
        lpcheck.saturate(molecule1);
        QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule1);
        Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(1),qAC));
	}
	/**
	 * A unit test suite for JUnit: Resonance Fluorobenzene  Fc1ccccc1 <=> ...
	 *
	 * @return    The test suite
	 */
	public void testResonanceFluorobenzene() throws ClassNotFoundException, CDKException, java.lang.Exception {
		IMolecule molecule = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("Fc1ccccc1");
        addImplicitHydrogens(molecule);
        LonePairElectronChecker lpcheck = new LonePairElectronChecker();
        lpcheck.saturate(molecule);
        
        StructureResonanceGenerator gRI = new StructureResonanceGenerator();
		IAtomContainerSet setOfMolecules = gRI.getAllStructures(molecule);
        
		Assert.assertEquals(4,setOfMolecules.getAtomContainerCount());

        IMolecule molecule1 = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("[F+]=C1C=CC=C[C-]1");
        addImplicitHydrogens(molecule1);
        lpcheck.saturate(molecule1);
        QueryAtomContainer qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule1);
        Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(1),qAC));
        
        IMolecule molecule2 = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("[F+]=C1-C=C-[C-]-C=C1");
        addImplicitHydrogens(molecule2);
        lpcheck.saturate(molecule2);
        qAC = QueryAtomContainerCreator.createSymbolAndChargeQueryContainer(molecule2);
        Assert.assertTrue(UniversalIsomorphismTester.isIsomorph(setOfMolecules.getAtomContainer(2),qAC));
	}
	/**
	 * A unit test suite for JUnit: Resonance   n1ccccc1 <=> ...
	 *
	 * @return    The test suite
	 */
	public void test_n1ccccc1() throws ClassNotFoundException, CDKException, java.lang.Exception {
		IMolecule molecule = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("n1ccccc1");
        addImplicitHydrogens(molecule);
        LonePairElectronChecker lpcheck = new LonePairElectronChecker();
        lpcheck.saturate(molecule);
        
        StructureResonanceGenerator gRI = new StructureResonanceGenerator();
		 IAtomContainerSet setOfMolecules = gRI.getAllStructures(molecule);
        
		Assert.assertEquals(10,setOfMolecules.getAtomContainerCount());
	}
}