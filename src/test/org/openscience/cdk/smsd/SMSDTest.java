/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smsd.factory.SubGraphFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 * 
 * @cdk.module test-smsd
 */
public class SMSDTest extends CDKTestCase {

    private static IMolecule Napthalene;
    private static IMolecule Cyclohexane;
    private static IMolecule Benzene;

    @BeforeClass
    public static void setUp() throws CDKException {
        Napthalene = createNaphthalene();
        Cyclohexane = createCyclohexane();
        Benzene = createBenzene();
    }

    @Test
    public void testVFLib() throws Exception {
        SubGraphFactory sbf = new SubGraphFactory(false, true, true, true);
        sbf.init(Benzene, Benzene, true);
        Assert.assertEquals(true, sbf.isSubgraph());

    }

    @Test
    public void testCDKMCS() throws Exception {
        SMSD ebimcs = new SMSD(3, false, true, true, true);
        ebimcs.init(Cyclohexane, Benzene, true);
        Assert.assertEquals(true, ebimcs.isSubgraph());
    }

    @Test
    public void testMCSPlus() throws Exception {
        SMSD ebimcs = new SMSD(1, false, true, true, true);
        ebimcs.init(Cyclohexane, Benzene, true);
        Assert.assertEquals(true, ebimcs.isSubgraph());
    }

    @Test
    public void testSMSD() throws Exception {
        SMSD ebimcs = new SMSD(0, false, true, true, true);
        ebimcs.init(Cyclohexane, Benzene, true);
        Assert.assertEquals(6, ebimcs.getFirstMapping().size());

        SMSD ebimcs1 = new SMSD(0, true, true, true, true);
        ebimcs1.init(Benzene, Napthalene, true);
        Assert.assertEquals(6, ebimcs.getFirstAtomMapping().size());
    }

    @Test
    public void testSMSDBondSensitive() throws Exception {
        SMSD ebimcs = new SMSD(3, false, true, true, true);
        ebimcs.init(Cyclohexane, Benzene, true);
        Assert.assertEquals(true, ebimcs.isSubgraph());

        SMSD ebimcs1 = new SMSD(3, true, true, true, true);
        ebimcs1.init(Cyclohexane, Benzene, true);
        Assert.assertEquals(false, ebimcs1.isSubgraph());
    }

    @Test
    public void testSMSDBondInSensitive() throws Exception {
        SMSD ebimcs = new SMSD(3, false, true, true, true);
        ebimcs.init(Cyclohexane, Benzene, true);
        SMSD ebimcs1 = new SMSD(3, true, true, true, true);
        ebimcs1.init(Benzene, Napthalene, true);
        Assert.assertEquals(6, ebimcs.getFirstAtomMapping().size());
        Assert.assertEquals(6, ebimcs1.getFirstAtomMapping().size());
    }

    @Test
    public void testSMSDChemicalFilters() throws Exception {
        SMSD ebimcs = new SMSD(3, false, true, true, true);
        ebimcs.init(Cyclohexane, Benzene, true);
        Assert.assertEquals(6, ebimcs.getFirstAtomMapping().size());
        Assert.assertEquals(true, ebimcs.isSubgraph());
    }

    @Test
    public void testSMSDScores() throws Exception {
        SMSD ebimcs = new SMSD(3, false, true, true, true);
        ebimcs.init(Cyclohexane, Benzene, true);
        Assert.assertEquals(Cyclohexane.getAtomCount(), ebimcs.getFirstMapping().size());
    }

    private IMolecule create4Toluene() throws CDKException {
        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();
        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c1.setID("1");
        IAtom c2 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c2.setID("2");
        IAtom c3 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c3.setID("3");
        IAtom c4 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c4.setID("4");
        IAtom c5 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c5.setID("5");
        IAtom c6 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c6.setID("6");
        IAtom c7 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c7.setID("7");

        result.addAtom(c1);
        result.addAtom(c2);
        result.addAtom(c3);
        result.addAtom(c4);
        result.addAtom(c5);
        result.addAtom(c6);
        result.addAtom(c7);



        IBond bond1 = new Bond(c1, c2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(c2, c3, IBond.Order.DOUBLE);
        IBond bond3 = new Bond(c3, c4, IBond.Order.SINGLE);
        IBond bond4 = new Bond(c4, c5, IBond.Order.DOUBLE);
        IBond bond5 = new Bond(c5, c6, IBond.Order.SINGLE);
        IBond bond6 = new Bond(c6, c1, IBond.Order.DOUBLE);
        IBond bond7 = new Bond(c7, c4, IBond.Order.SINGLE);

        result.addBond(bond1);
        result.addBond(bond2);
        result.addBond(bond3);
        result.addBond(bond4);
        result.addBond(bond5);
        result.addBond(bond6);
        result.addBond(bond7);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(result);

        return result;
    }

    public IMolecule createMethane() {
        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();
        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        result.addAtom(c1);

        return result;
    }

    public IMolecule createPropane() {
        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();
        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        IAtom c2 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        IAtom c3 = DefaultChemObjectBuilder.getInstance().newAtom("C");


        result.addAtom(c1);
        result.addAtom(c2);
        result.addAtom(c3);



        IBond bond1 = new Bond(c1, c2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(c2, c3, IBond.Order.SINGLE);

        result.addBond(bond1);
        result.addBond(bond2);



        return result;
    }

    public IMolecule createHexane() throws CDKException {
        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();

        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c1.setID("1");
        IAtom c2 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c2.setID("2");
        IAtom c3 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c3.setID("3");
        IAtom c4 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c4.setID("4");
        IAtom c5 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c5.setID("5");
        IAtom c6 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c6.setID("6");

        result.addAtom(c1);
        result.addAtom(c2);
        result.addAtom(c3);
        result.addAtom(c4);
        result.addAtom(c5);
        result.addAtom(c6);

        IBond bond1 = new Bond(c1, c2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(c2, c3, IBond.Order.SINGLE);
        IBond bond3 = new Bond(c3, c4, IBond.Order.SINGLE);
        IBond bond4 = new Bond(c4, c5, IBond.Order.SINGLE);
        IBond bond5 = new Bond(c5, c6, IBond.Order.SINGLE);


        result.addBond(bond1);
        result.addBond(bond2);
        result.addBond(bond3);
        result.addBond(bond4);
        result.addBond(bond5);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(result);

        return result;
    }

    public static IMolecule createBenzene() throws CDKException {
        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();

        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c1.setID("1");
        IAtom c2 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c2.setID("2");
        IAtom c3 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c3.setID("3");
        IAtom c4 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c4.setID("4");
        IAtom c5 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c5.setID("5");
        IAtom c6 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c6.setID("6");

        result.addAtom(c1);
        result.addAtom(c2);
        result.addAtom(c3);
        result.addAtom(c4);
        result.addAtom(c5);
        result.addAtom(c6);

        IBond bond1 = new Bond(c1, c2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(c2, c3, IBond.Order.DOUBLE);
        IBond bond3 = new Bond(c3, c4, IBond.Order.SINGLE);
        IBond bond4 = new Bond(c4, c5, IBond.Order.DOUBLE);
        IBond bond5 = new Bond(c5, c6, IBond.Order.SINGLE);
        IBond bond6 = new Bond(c6, c1, IBond.Order.DOUBLE);


        result.addBond(bond1);
        result.addBond(bond2);
        result.addBond(bond3);
        result.addBond(bond4);
        result.addBond(bond5);
        result.addBond(bond6);


        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(result);



        return result;
    }
    //    public static Molecule createPyridine() {
//        Molecule result = new DefaultMolecule();
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//        Atom c6 = result.addAtom("N");
//
//        result.connect(c1, c2, 1);
//        result.connect(c2, c3, 2);
//        result.connect(c3, c4, 1);
//        result.connect(c4, c5, 2);
//        result.connect(c5, c6, 1);
//        result.connect(c6, c1, 2);
//
//        return result;
//    }
//
//    public static Molecule createToluene() {
//        Molecule result = new DefaultMolecule();
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//        Atom c6 = result.addAtom("C");
//        Atom c7 = result.addAtom("C");
//
//        result.connect(c1, c2, 1);
//        result.connect(c2, c3, 2);
//        result.connect(c3, c4, 1);
//        result.connect(c4, c5, 2);
//        result.connect(c5, c6, 1);
//        result.connect(c6, c1, 2);
//        result.connect(c7, c1, 1);
//
//        return result;
//    }
//
//    public static Molecule createPhenol() {
//        Molecule result = new DefaultMolecule();
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//        Atom c6 = result.addAtom("C");
//        Atom c7 = result.addAtom("O");
//
//        result.connect(c1, c2, 1);
//        result.connect(c2, c3, 2);
//        result.connect(c3, c4, 1);
//        result.connect(c4, c5, 2);
//        result.connect(c5, c6, 1);
//        result.connect(c6, c1, 2);
//        result.connect(c7, c1, 1);
//
//        return result;
//    }
//

    public static IMolecule createNaphthalene() throws CDKException {
        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();

        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c1.setID("1");
        IAtom c2 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c2.setID("2");
        IAtom c3 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c3.setID("3");
        IAtom c4 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c4.setID("4");
        IAtom c5 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c5.setID("5");
        IAtom c6 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c6.setID("6");
        IAtom c7 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c4.setID("7");
        IAtom c8 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c5.setID("8");
        IAtom c9 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c6.setID("9");
        IAtom c10 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c6.setID("10");

        result.addAtom(c1);
        result.addAtom(c2);
        result.addAtom(c3);
        result.addAtom(c4);
        result.addAtom(c5);
        result.addAtom(c6);
        result.addAtom(c7);
        result.addAtom(c8);
        result.addAtom(c9);
        result.addAtom(c10);



        IBond bond1 = new Bond(c1, c2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(c2, c3, IBond.Order.DOUBLE);
        IBond bond3 = new Bond(c3, c4, IBond.Order.SINGLE);
        IBond bond4 = new Bond(c4, c5, IBond.Order.DOUBLE);
        IBond bond5 = new Bond(c5, c6, IBond.Order.SINGLE);
        IBond bond6 = new Bond(c6, c1, IBond.Order.DOUBLE);
        IBond bond7 = new Bond(c5, c7, IBond.Order.SINGLE);
        IBond bond8 = new Bond(c7, c8, IBond.Order.DOUBLE);
        IBond bond9 = new Bond(c8, c9, IBond.Order.SINGLE);
        IBond bond10 = new Bond(c9, c10, IBond.Order.DOUBLE);
        IBond bond11 = new Bond(c10, c6, IBond.Order.SINGLE);


        result.addBond(bond1);
        result.addBond(bond2);
        result.addBond(bond3);
        result.addBond(bond4);
        result.addBond(bond5);
        result.addBond(bond6);
        result.addBond(bond7);
        result.addBond(bond8);
        result.addBond(bond9);
        result.addBond(bond10);
        result.addBond(bond11);


        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(result);



        return result;
    }
//
//    public static Molecule createAcetone() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom o3 = result.addAtom("O");
//
//        result.connect(c0, c1, 1);
//        result.connect(c1, c2, 1);
//        result.connect(c1, o3, 2);
//
//        return result;
//    }
//
//    public static Molecule createNeopentane() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//
//        result.connect(c0, c1, 1);
//        result.connect(c0, c2, 1);
//        result.connect(c0, c3, 1);
//        result.connect(c0, c4, 1);
//
//        return result;
//    }
//
//    public static Molecule createCubane() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//        Atom c6 = result.addAtom("C");
//        Atom c7 = result.addAtom("C");
//
//        result.connect(c0, c1, 1);
//        result.connect(c1, c2, 1);
//        result.connect(c2, c3, 1);
//        result.connect(c3, c0, 1);
//
//        result.connect(c4, c5, 1);
//        result.connect(c5, c6, 1);
//        result.connect(c6, c7, 1);
//        result.connect(c7, c4, 1);
//
//        result.connect(c0, c4, 1);
//        result.connect(c1, c5, 1);
//        result.connect(c2, c6, 1);
//        result.connect(c3, c7, 1);
//
//        return result;
//    }
//
//    public static Molecule createBicyclo220hexane() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//        Atom c3 = result.addAtom("C");
//        Atom c4 = result.addAtom("C");
//        Atom c5 = result.addAtom("C");
//
//        result.connect(c0, c1, 1);
//        result.connect(c1, c2, 1);
//        result.connect(c2, c3, 1);
//        result.connect(c3, c4, 1);
//        result.connect(c4, c5, 1);
//        result.connect(c5, c0, 1);
//        result.connect(c2, c5, 1);
//
//        return result;
//    }
//
//    public static Molecule createEthylbenzeneWithSuperatom() {
//        Molecule result = Molecules.createBenzene();
//        Atom carbon1 = result.addAtom("C");
//        Atom carbon2 = result.addAtom("C");
//        Bond crossingBond = result.connect(result.getAtom(0), carbon1, 1);
//        result.connect(carbon1, carbon2, 1);
//
//        Superatom substructure = result.addSuperatom();
//        substructure.addAtom(carbon1);
//        substructure.addAtom(carbon2);
//        substructure.addCrossingBond(crossingBond);
//        substructure.setCrossingVector(crossingBond, 0.1, 0.1);
//        substructure.setLabel("Ethyl");
//
//        return result;
//    }
//

    public static IMolecule createCyclohexane() throws CDKException {

        IMolecule result = DefaultChemObjectBuilder.getInstance().newMolecule();

        IAtom c1 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c1.setID("1");
        IAtom c2 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c2.setID("2");
        IAtom c3 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c3.setID("3");
        IAtom c4 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c4.setID("4");
        IAtom c5 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c5.setID("5");
        IAtom c6 = DefaultChemObjectBuilder.getInstance().newAtom("C");
        c6.setID("6");

        result.addAtom(c1);
        result.addAtom(c2);
        result.addAtom(c3);
        result.addAtom(c4);
        result.addAtom(c5);
        result.addAtom(c6);

        IBond bond1 = new Bond(c1, c2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(c2, c3, IBond.Order.SINGLE);
        IBond bond3 = new Bond(c3, c4, IBond.Order.SINGLE);
        IBond bond4 = new Bond(c4, c5, IBond.Order.SINGLE);
        IBond bond5 = new Bond(c5, c6, IBond.Order.SINGLE);
        IBond bond6 = new Bond(c6, c1, IBond.Order.SINGLE);


        result.addBond(bond1);
        result.addBond(bond2);
        result.addBond(bond3);
        result.addBond(bond4);
        result.addBond(bond5);
        result.addBond(bond6);


        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(result);



        return result;

    }
//
    //    public static Molecule createCyclopropane() {
//        Molecule result = new DefaultMolecule();
//        Atom c0 = result.addAtom("C");
//        Atom c1 = result.addAtom("C");
//        Atom c2 = result.addAtom("C");
//
//        result.connect(c0, c1, 1);
//        result.connect(c1, c2, 1);
//        result.connect(c2, c0, 1);
//
//        return result;
//    }
}
