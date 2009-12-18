/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.algorithm.vflib;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import junit.framework.Assert;
import org.junit.BeforeClass;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IMapper;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.INode;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IQuery;
import org.openscience.cdk.smsd.algorithm.vflib.interfaces.IState;
import org.openscience.cdk.smsd.algorithm.vflib.map.VFMapper;
import org.openscience.cdk.smsd.algorithm.vflib.map.VFState;
import org.openscience.cdk.smsd.algorithm.vflib.query.TemplateCompiler;
import org.openscience.cdk.smsd.algorithm.vflib.validator.VFMatch;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @cdk.module test-smsd
 */
public class VFLibTest extends CDKTestCase {

    private static IAtomContainer hexane;
    private static IQuery hexaneQuery;
    private static IAtomContainer benzene;
    private static IQuery benzeneQuery;

    @BeforeClass
    public static void setUp() throws CDKException {
        hexane = createHexane();
        System.out.println(" " + hexane.getAtomCount());
        hexaneQuery = TemplateCompiler.compile(hexane);
        System.out.println(" " + hexaneQuery.countNodes());
        benzene = createBenzene();
        benzeneQuery = TemplateCompiler.compile(benzene);
    }

    public void testItShouldFindAllMatchCandidatesInTheRootState() {
        IState state = new VFState(benzeneQuery, benzene);
        int count = 0;

        while (state.hasNextCandidate()) {
            state.nextCandidate();
            count++;
        }
        Assert.assertEquals(benzene.getAtomCount() * benzene.getAtomCount(), count);
    }

    public void testItShoudFindAllMatchCandidatesInThePrimaryState() {
        IState state = new VFState(benzeneQuery, benzene);
        VFMatch match = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState newState = state.nextState(match);
        List<VFMatch> candidates = new ArrayList<VFMatch>();

        while (newState.hasNextCandidate()) {
            candidates.add(newState.nextCandidate());
        }

        Assert.assertEquals(4, candidates.size());
    }

    public void testItShouldFindAllMatchCandidatesInTheSecondaryState() {
        IState state0 = new VFState(benzeneQuery, benzene);
        VFMatch match0 = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState state1 = state0.nextState(match0);
        VFMatch match1 = new VFMatch(benzeneQuery.getNode(1), benzene.getAtom(1));
        IState state2 = state1.nextState(match1);
        List<VFMatch> candidates = new ArrayList<VFMatch>();

        while (state2.hasNextCandidate()) {
            candidates.add(state2.nextCandidate());
        }

        Assert.assertEquals(1, candidates.size());
    }

    public void testItShouldMapAllAtomsInTheSecondaryState() {
        IState state0 = new VFState(benzeneQuery, benzene);
        VFMatch match0 = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState state1 = state0.nextState(match0);
        VFMatch match1 = new VFMatch(benzeneQuery.getNode(1), benzene.getAtom(1));
        IState state2 = state1.nextState(match1);

        Map<INode, IAtom> map = state2.getMap();

        Assert.assertEquals(2, map.size());
        Assert.assertEquals(benzene.getAtom(0), map.get(benzeneQuery.getNode(0)));
        Assert.assertEquals(benzene.getAtom(1), map.get(benzeneQuery.getNode(1)));
    }

    public void testItShouldFindAllMatchCandidatesFromTheTeriaryState() {
        IState state0 = new VFState(benzeneQuery, benzene);
        VFMatch match0 = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState state1 = state0.nextState(match0);
        VFMatch match1 = new VFMatch(benzeneQuery.getNode(1), benzene.getAtom(1));
        IState state2 = state1.nextState(match1);
        VFMatch match2 = new VFMatch(benzeneQuery.getNode(2), benzene.getAtom(2));
        IState state3 = state2.nextState(match2);
        List<VFMatch> candidates = new ArrayList<VFMatch>();

        while (state3.hasNextCandidate()) {
            candidates.add(state3.nextCandidate());
        }

        Assert.assertEquals(1, candidates.size());
    }

    public void testItShouldMapAllAtomsInTheTertiaryState() {
        IState state0 = new VFState(benzeneQuery, benzene);
        VFMatch match0 = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState state1 = state0.nextState(match0);
        VFMatch match1 = new VFMatch(benzeneQuery.getNode(1), benzene.getAtom(1));
        IState state2 = state1.nextState(match1);
        VFMatch match2 = new VFMatch(benzeneQuery.getNode(2), benzene.getAtom(2));
        IState state3 = state2.nextState(match2);
        Map<INode, IAtom> map = state3.getMap();

        Assert.assertEquals(3, map.size());
        Assert.assertEquals(benzene.getAtom(0), map.get(benzeneQuery.getNode(0)));
        Assert.assertEquals(benzene.getAtom(1), map.get(benzeneQuery.getNode(1)));
        Assert.assertEquals(benzene.getAtom(2), map.get(benzeneQuery.getNode(2)));
    }

    public void testItShouldReachGoalWhenAllAtomsAreMapped() {
        IState state0 = new VFState(benzeneQuery, benzene);
        VFMatch match0 = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));
        IState state1 = state0.nextState(match0);
        VFMatch match1 = new VFMatch(benzeneQuery.getNode(1), benzene.getAtom(1));
        IState state2 = state1.nextState(match1);
        VFMatch match2 = new VFMatch(benzeneQuery.getNode(2), benzene.getAtom(2));
        IState state3 = state2.nextState(match2);
        VFMatch match3 = new VFMatch(benzeneQuery.getNode(3), benzene.getAtom(3));
        IState state4 = state3.nextState(match3);
        VFMatch match4 = new VFMatch(benzeneQuery.getNode(4), benzene.getAtom(4));
        IState state5 = state4.nextState(match4);

        Assert.assertFalse(state5.isGoal());

        VFMatch match5 = new VFMatch(benzeneQuery.getNode(5), benzene.getAtom(5));
        IState state6 = state5.nextState(match5);

        Assert.assertTrue(state6.isGoal());
    }

    public void testItShouldHaveANextCandidateInTheSecondaryState() {
        IState state = new VFState(benzeneQuery, benzene);
        VFMatch match = new VFMatch(benzeneQuery.getNode(0), benzene.getAtom(0));

        IState nextState = state.nextState(match);

        Assert.assertTrue(nextState.hasNextCandidate());
    }

    public void testItShouldMatchHexaneToHexane() {
        IMapper mapper = new VFMapper(hexaneQuery);

        Assert.assertTrue(mapper.hasMap(hexane));
    }

    public void testItShouldMatchHexaneToHexaneWhenUsingMolecule() {
        IMapper mapper = new VFMapper(hexane);

        Assert.assertTrue(mapper.hasMap(hexane));
    }

    public void testItShouldFindTwoMapsFromHexaneToHexane() {
        IMapper mapper = new VFMapper(hexaneQuery);

        List<Map<INode, IAtom>> maps = mapper.getMaps(hexane);
        Assert.assertEquals(2, maps.size());
    }

    public static IMolecule createHexane() throws CDKException {
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
}
