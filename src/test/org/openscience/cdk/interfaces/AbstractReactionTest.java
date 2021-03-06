/* $RCSfile$
 * $Author$    
 * $Date$    
 * $Revision$
 * 
 * Copyright (C) 2003-2007  The Chemistry Development Kit (CDK) project
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
package org.openscience.cdk.interfaces;

import java.util.Iterator;

import org.junit.Assert;
import org.junit.Test;

/**
 * TestCase for {@link IReaction} implementations.
 *
 * @cdk.module test-interfaces
 */
public abstract class AbstractReactionTest extends AbstractChemObjectTest {

    @Test public void testGetReactantCount() {
        IReaction reaction = (IReaction)newChemObject();
        Assert.assertEquals(0, reaction.getReactantCount());
	reaction.addReactant(reaction.getBuilder().newMolecule());
        Assert.assertEquals(1, reaction.getReactantCount());
    }
    
    @Test public void testGetProductCount() {
        IReaction reaction = (IReaction)newChemObject();
        Assert.assertEquals(0, reaction.getProductCount());
	reaction.addProduct(reaction.getBuilder().newMolecule());
        Assert.assertEquals(1, reaction.getProductCount());
    }
    
    @Test public void testAddReactant_IMolecule() {
        IReaction reaction = (IReaction)newChemObject();
        IMolecule sodiumhydroxide = reaction.getBuilder().newMolecule();
        IMolecule aceticAcid = reaction.getBuilder().newMolecule();
        IMolecule water = reaction.getBuilder().newMolecule();
        IMolecule acetate = reaction.getBuilder().newMolecule();
        reaction.addReactant(sodiumhydroxide);
        reaction.addReactant(aceticAcid);
        reaction.addReactant(water);
        Assert.assertEquals(3, reaction.getReactantCount());
        // next one should trigger a growArray, if the grow
        // size is still 3.
        reaction.addReactant(acetate);
        Assert.assertEquals(4, reaction.getReactantCount());
        
        Assert.assertEquals(1.0, reaction.getReactantCoefficient(aceticAcid), 0.00001);
    }

    @Test public void testSetReactants_IMoleculeSet() {
        IReaction reaction = (IReaction)newChemObject();
        IMolecule sodiumhydroxide = reaction.getBuilder().newMolecule();
        IMolecule aceticAcid = reaction.getBuilder().newMolecule();
        IMolecule water = reaction.getBuilder().newMolecule();
        IMoleculeSet reactants = reaction.getBuilder().newMoleculeSet();
        reactants.addMolecule(sodiumhydroxide);
        reactants.addMolecule(aceticAcid);
        reactants.addMolecule(water);
        reaction.setReactants(reactants);
        Assert.assertEquals(3, reaction.getReactantCount());
        
        Assert.assertEquals(1.0, reaction.getReactantCoefficient(aceticAcid), 0.00001);
    }

    @Test public void testAddReactant_IMolecule_Double() {
        IReaction reaction = (IReaction)newChemObject();
        IMolecule proton = reaction.getBuilder().newMolecule();
        IMolecule sulfate = reaction.getBuilder().newMolecule();
        reaction.addReactant(proton, 2d);
        reaction.addReactant(sulfate, 1d);
        Assert.assertEquals(2.0, reaction.getReactantCoefficient(proton), 0.00001);
        Assert.assertEquals(1.0, reaction.getReactantCoefficient(sulfate), 0.00001);
    }
    
    @Test public void testAddProduct_IMolecule() {
        IReaction reaction = (IReaction)newChemObject();
        IMolecule sodiumhydroxide = reaction.getBuilder().newMolecule();
        IMolecule aceticAcid = reaction.getBuilder().newMolecule();
        IMolecule water = reaction.getBuilder().newMolecule();
        IMolecule acetate = reaction.getBuilder().newMolecule();
        reaction.addProduct(sodiumhydroxide);
        reaction.addProduct(aceticAcid);
        reaction.addProduct(water);
        Assert.assertEquals(3, reaction.getProductCount());
        // next one should trigger a growArray, if the grow
        // size is still 3.
        reaction.addProduct(acetate);
        Assert.assertEquals(4, reaction.getProductCount());
        
        Assert.assertEquals(1.0, reaction.getProductCoefficient(aceticAcid), 0.00001);
    }

    @Test public void testSetProducts_IMoleculeSet() {
        IReaction reaction = (IReaction)newChemObject();
        IMolecule sodiumhydroxide = reaction.getBuilder().newMolecule();
        IMolecule aceticAcid = reaction.getBuilder().newMolecule();
        IMolecule water = reaction.getBuilder().newMolecule();
        IMoleculeSet products = reaction.getBuilder().newMoleculeSet();
        products.addMolecule(sodiumhydroxide);
        products.addMolecule(aceticAcid);
        products.addMolecule(water);
        reaction.setProducts(products);
        Assert.assertEquals(3, reaction.getProductCount());
        
        Assert.assertEquals(1.0, reaction.getProductCoefficient(aceticAcid), 0.00001);
    }

    @Test public void testAddProduct_IMolecule_Double() {
        IReaction reaction = (IReaction)newChemObject();
        IMolecule proton = reaction.getBuilder().newMolecule();
        IMolecule sulfate = reaction.getBuilder().newMolecule();
        reaction.addProduct(proton, 2.0);
        reaction.addProduct(sulfate, 1.0);
        Assert.assertEquals(2.0, reaction.getProductCoefficient(proton), 0.00001);
        Assert.assertEquals(1.0, reaction.getProductCoefficient(sulfate), 0.00001);
    }
    
    @Test public void testAddAgent_IMolecule() {
        IReaction reaction = (IReaction)newChemObject();
        IMolecule proton = reaction.getBuilder().newMolecule();
        reaction.addAgent(proton);
        Assert.assertEquals(1, reaction.getAgents().getMoleculeCount());
    }

    @Test public void testGetReactantCoefficient_IMolecule() {
        IReaction reaction = (IReaction)newChemObject();
        IMolecule proton = reaction.getBuilder().newMolecule();
        reaction.addReactant(proton, 2.0);
        Assert.assertEquals(2.0, reaction.getReactantCoefficient(proton), 0.00001);
        
        Assert.assertEquals(-1.0, reaction.getReactantCoefficient(reaction.getBuilder().newMolecule()), 0.00001);
    }

    @Test public void testGetProductCoefficient_IMolecule() {
        IReaction reaction = (IReaction)newChemObject();
        IMolecule proton = reaction.getBuilder().newMolecule();
        reaction.addProduct(proton, 2.0);
        Assert.assertEquals(2.0, reaction.getProductCoefficient(proton), 0.00001);

        Assert.assertEquals(-1.0, reaction.getProductCoefficient(reaction.getBuilder().newMolecule()), 0.00001);
    }
    
	@Test public void testSetReactantCoefficient_IMolecule_Double() {
		IReaction reaction = (IReaction)newChemObject();
        IMolecule proton = reaction.getBuilder().newMolecule();
		reaction.addReactant(proton, 2.0);
		reaction.setReactantCoefficient(proton, 3.0);
		Assert.assertEquals(3.0, reaction.getReactantCoefficient(proton), 0.00001);
	}
	
	@Test public void testSetProductCoefficient_IMolecule_Double() {
		IReaction reaction = (IReaction)newChemObject();
        IMolecule proton = reaction.getBuilder().newMolecule();
		reaction.addProduct(proton, 2.0);
		reaction.setProductCoefficient(proton, 1.0);
		Assert.assertEquals(1.0, reaction.getProductCoefficient(proton), 0.00001);
	}
	
	@Test public void testGetReactantCoefficients() {
        IReaction reaction = (IReaction)newChemObject();
		IMolecule ed1 = reaction.getBuilder().newMolecule();
		IMolecule ed2 = reaction.getBuilder().newMolecule();
		reaction.addReactant(ed1, 2d);
		reaction.addReactant(ed2, 3d);
		Double[] ec = reaction.getReactantCoefficients();
		Assert.assertEquals(2.0, ec.length, 0.00001);
		Assert.assertEquals(reaction.getReactantCoefficient(ed1), ec[0], 0.00001);
		Assert.assertEquals(3.0, ec[1], 0.00001);
    }
	
	@Test public void testGetProductCoefficients() {
        IReaction reaction = (IReaction)newChemObject();
		IMolecule pr1 = reaction.getBuilder().newMolecule();
		IMolecule pr2 = reaction.getBuilder().newMolecule();
		reaction.addProduct(pr1, 1d);
		reaction.addProduct(pr2, 2d);
		Double[] pc = reaction.getProductCoefficients();
		Assert.assertEquals(2.0, pc.length, 0.00001);
		Assert.assertEquals(reaction.getProductCoefficient(pr1), pc[0], 0.00001);
		Assert.assertEquals(2.0, pc[1], 0.00001);
    }
	
	@Test public void testSetReactantCoefficients_arrayDouble() {
        IReaction reaction = (IReaction)newChemObject();
		IMolecule ed1 = reaction.getBuilder().newMolecule();
		IMolecule ed2 = reaction.getBuilder().newMolecule();
		reaction.addReactant(ed1, 2d);
		reaction.addReactant(ed2, 3d);
		Double[] ec = { 1.0, 2.0 };
		boolean coeffSet = reaction.setReactantCoefficients(ec);
		Assert.assertTrue(coeffSet);
		Assert.assertEquals(1.0, reaction.getReactantCoefficient(ed1), 0.00001);
		Assert.assertEquals(2.0, reaction.getReactantCoefficient(ed2), 0.00001);
		Double[] ecFalse = { 1.0 };
		Assert.assertFalse(reaction.setReactantCoefficients(ecFalse));
    }
	
	@Test public void testSetProductCoefficients_arrayDouble() {
        IReaction reaction = (IReaction)newChemObject();
		IMolecule pr1 = reaction.getBuilder().newMolecule();
		reaction.addProduct(pr1, 1d);
		Double[] pc = { 2.0 };
		boolean coeffSet = reaction.setProductCoefficients(pc);
		Assert.assertTrue(coeffSet);
		Assert.assertEquals(2.0, reaction.getProductCoefficient(pr1), 0.00001);
		Double[] pcFalse = { 1.0 , 2.0 };
		Assert.assertFalse(reaction.setProductCoefficients(pcFalse));
    }
	
    @Test public void testGetReactants() {
        IReaction reaction = (IReaction)newChemObject();
        IMolecule sodiumhydroxide = reaction.getBuilder().newMolecule();
        IMolecule aceticAcid = reaction.getBuilder().newMolecule();
        IMolecule water = reaction.getBuilder().newMolecule();
        reaction.addReactant(sodiumhydroxide);
        reaction.addReactant(aceticAcid);
        reaction.addReactant(water);
        Assert.assertEquals(3, reaction.getReactants().getMoleculeCount());
    }
	
    @Test public void testGetProducts() {    
        IReaction reaction = (IReaction)newChemObject();
        IMolecule sodiumhydroxide = reaction.getBuilder().newMolecule();
        IMolecule aceticAcid = reaction.getBuilder().newMolecule();
        IMolecule water = reaction.getBuilder().newMolecule();
        reaction.addProduct(sodiumhydroxide);
        reaction.addProduct(aceticAcid);
        reaction.addProduct(water);
        Assert.assertEquals(3, reaction.getProducts().getMoleculeCount());
    }
    
    @Test public void testGetAgents() {    
        IReaction reaction = (IReaction)newChemObject();
        IMolecule water = reaction.getBuilder().newMolecule();
        reaction.addAgent(water);
        Assert.assertEquals(1, reaction.getAgents().getMoleculeCount());
    }
    
    @Test public void testSetDirection_IReaction_Direction() {
        IReaction reaction = (IReaction)newChemObject();
        IReaction.Direction direction = IReaction.Direction.BIDIRECTIONAL;
        reaction.setDirection(direction);
        Assert.assertEquals(direction, reaction.getDirection());
    }

    @Test public void testGetDirection() {
        IReaction reaction = (IReaction)newChemObject();
        Assert.assertEquals(IReaction.Direction.FORWARD, reaction.getDirection());
    }

    /**
     * Method to test whether the class complies with RFC #9.
     */
    @Test public void testToString() {
        IReaction reaction = (IReaction)newChemObject();
        String description = reaction.toString();
        for (int i=0; i< description.length(); i++) {
            Assert.assertTrue(description.charAt(i) != '\n');
            Assert.assertTrue(description.charAt(i) != '\r');
        }
    }
    
    @Test public void testClone() throws Exception {
        IReaction reaction = (IReaction)newChemObject();
        Object clone = reaction.clone();
        Assert.assertNotNull(clone);
        Assert.assertTrue(clone instanceof IReaction);
    }
    
    @Test public void testClone_Mapping() throws Exception {
        IReaction reaction = (IReaction)newChemObject();
        IMapping mapping = reaction.getBuilder().newMapping(reaction.getBuilder().newAtom("C"), reaction.getBuilder().newAtom("C"));
        reaction.addMapping(mapping);
        IReaction clonedReaction = (IReaction)reaction.clone();
        Iterator<IMapping> mappings = reaction.mappings().iterator();
        Iterator<IMapping> clonedMappings = clonedReaction.mappings().iterator();
        Assert.assertNotNull(mappings);
        Assert.assertTrue(mappings.hasNext());
        Assert.assertNotNull(clonedMappings);
        Assert.assertTrue(clonedMappings.hasNext());
    }
    
    @Test public void testAddMapping_IMapping() {
        IReaction reaction = (IReaction)newChemObject();
        IMapping mapping = reaction.getBuilder().newMapping(reaction.getBuilder().newAtom("C"), reaction.getBuilder().newAtom("C"));
        reaction.addMapping(mapping);
        Iterator<IMapping> mappings = reaction.mappings().iterator();
        Assert.assertNotNull(mappings);
        Assert.assertTrue(mappings.hasNext());
        Assert.assertEquals(mapping, (IMapping)mappings.next());
    }
    
    @Test public void testRemoveMapping_int() {
        IReaction reaction = (IReaction)newChemObject();
        IMapping mapping = reaction.getBuilder().newMapping(reaction.getBuilder().newAtom("C"), reaction.getBuilder().newAtom("C"));
        reaction.addMapping(mapping);
        Assert.assertEquals(1, reaction.getMappingCount());
        reaction.removeMapping(0);
        Assert.assertEquals(0, reaction.getMappingCount());
    }
    
    @Test public void testGetMapping_int() {
    	IReaction reaction = (IReaction)newChemObject();
        IMapping mapping = reaction.getBuilder().newMapping(reaction.getBuilder().newAtom("C"), reaction.getBuilder().newAtom("C"));
        reaction.addMapping(mapping);
        IMapping gotIt = reaction.getMapping(0);
        Assert.assertEquals(mapping, gotIt);
    }
    
    @Test public void testGetMappingCount() {
    	IReaction reaction = (IReaction)newChemObject();
        IMapping mapping = reaction.getBuilder().newMapping(reaction.getBuilder().newAtom("C"), reaction.getBuilder().newAtom("C"));
        reaction.addMapping(mapping);
        Assert.assertEquals(1, reaction.getMappingCount());
    }
    
    @Test public void testMappings() {
    	IReaction reaction = (IReaction)newChemObject();
        IMapping mapping = reaction.getBuilder().newMapping(reaction.getBuilder().newAtom("C"), reaction.getBuilder().newAtom("C"));
        reaction.addMapping(mapping);
        Assert.assertEquals(1, reaction.getMappingCount());
    }
}
