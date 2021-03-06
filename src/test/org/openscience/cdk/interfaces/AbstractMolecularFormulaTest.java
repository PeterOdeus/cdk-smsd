/* $RCSfile$
 * $Author$    
 * $Date$    
 * $Revision$
 * 
 *  Copyright (C) 2007  Miguel Rojasch <miguelrojasch@users.sf.net>
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
 * 
 */
package org.openscience.cdk.interfaces;

import java.util.Iterator;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.CDKTestCase;

/**
 * Checks the functionality of {@link IMolecularFormula} implementations.
 *
 * @cdk.module test-interfaces
 */
public abstract class AbstractMolecularFormulaTest extends CDKTestCase {

    private static IChemObjectBuilder builder;

    public static IChemObjectBuilder getBuilder() {
        return builder;
    }

    public static void setBuilder( IChemObjectBuilder builder ) {
    	AbstractMolecularFormulaTest.builder = builder;
    }

    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetIsotopeCount0() {

        IMolecularFormula mf = getBuilder().newMolecularFormula();
    	
        Assert.assertEquals(0, mf.getIsotopeCount());
    }
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetIsotopeCount() {
    	
        IMolecularFormula mf = getBuilder().newMolecularFormula();
        mf.addIsotope( getBuilder().newIsotope("C") );
        mf.addIsotope( getBuilder().newIsotope("H") );
        mf.addIsotope( getBuilder().newIsotope("H") );
        mf.addIsotope( getBuilder().newIsotope("H") );
        mf.addIsotope( getBuilder().newIsotope("H") );
        
        Assert.assertEquals(2, mf.getIsotopeCount());
    }
    
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testAddIsotope_IIsotope() {
    	
    	IMolecularFormula mf = getBuilder().newMolecularFormula();
    	mf.addIsotope( getBuilder().newIsotope("C") );
        mf.addIsotope( getBuilder().newIsotope("H") );
        
    	IIsotope hy = getBuilder().newIsotope("C");
    	hy.setNaturalAbundance(2.00342342);
    	mf.addIsotope( hy);
        
        Assert.assertEquals(3, mf.getIsotopeCount());
    }
    
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetIsotopeCount_IIsotope() {
        IMolecularFormula mf = getBuilder().newMolecularFormula();
        
    	IIsotope carb = getBuilder().newIsotope("C");
    	IIsotope flu = getBuilder().newIsotope("F");
    	IIsotope h1 = getBuilder().newIsotope("H");
    	IIsotope h2 = getBuilder().newIsotope("H");
    	IIsotope h3 = getBuilder().newIsotope("H");
        mf.addIsotope( carb );
        mf.addIsotope( flu );
        mf.addIsotope( h1 );
        mf.addIsotope( h2 );
        mf.addIsotope( h3 );
        
        Assert.assertEquals(3, mf.getIsotopeCount());
        Assert.assertEquals(1, mf.getIsotopeCount(carb));
        Assert.assertEquals(1, mf.getIsotopeCount(flu));
        Assert.assertEquals(3, mf.getIsotopeCount(h1));
    }

    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetIsotopeCount_IIsotope2() {
        IMolecularFormula mf = getBuilder().newMolecularFormula();
        
        IIsotope carb = getBuilder().newIsotope("C");
        IIsotope flu = getBuilder().newIsotope("F");
        IIsotope h1 = getBuilder().newIsotope("H");
        mf.addIsotope( carb );
        mf.addIsotope( flu );
        mf.addIsotope( h1);
        mf.addIsotope( h1);
        mf.addIsotope( h1);
        
        Assert.assertEquals(3, mf.getIsotopeCount());
        Assert.assertEquals(1, mf.getIsotopeCount(carb));
        Assert.assertEquals(1, mf.getIsotopeCount(flu));
        Assert.assertEquals(3, mf.getIsotopeCount(h1));
    }
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testAddIsotope_IIsotope_int() {
        IMolecularFormula mf = getBuilder().newMolecularFormula();
        
        IIsotope carb = getBuilder().newIsotope("C");
        IIsotope flu = getBuilder().newIsotope("F");
        IIsotope h1 = getBuilder().newIsotope("H");
        mf.addIsotope( carb );
        mf.addIsotope( flu );
        mf.addIsotope( h1,3);
        
        Assert.assertEquals(3, mf.getIsotopeCount());
        Assert.assertEquals(1, mf.getIsotopeCount(carb));
        Assert.assertEquals(1, mf.getIsotopeCount(flu));
        Assert.assertEquals(3, mf.getIsotopeCount(h1));
        // In a List the objects are not stored in the same order than called
//        Assert.assertEquals("C", mf.getIsotope(0).getSymbol());
//        Assert.assertEquals("F", mf.getIsotope(1).getSymbol());
//        Assert.assertEquals("H", mf.getIsotope(2).getSymbol());
    }
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetIsotope_Number_Clone() throws Exception {
        IMolecularFormula mf = getBuilder().newMolecularFormula();
        
        IIsotope carb = getBuilder().newIsotope("C");
        IIsotope flu = getBuilder().newIsotope("F");
        IIsotope h1 = getBuilder().newIsotope("H");
        mf.addIsotope( carb );
        mf.addIsotope( flu );
        mf.addIsotope( h1,3);
        
        Object clone = mf.clone();
        Assert.assertTrue(clone instanceof IMolecularFormula);
        
        IMolecularFormula cloneFormula =(IMolecularFormula) clone;

        Assert.assertEquals(1, cloneFormula.getIsotopeCount(carb));
        Assert.assertEquals(1, cloneFormula.getIsotopeCount(flu));
        Assert.assertEquals(3, cloneFormula.getIsotopeCount(h1));
        // In a List the objects are not stored in the same order than called
//        Assert.assertEquals("C", cloneFormula.getIsotope(0).getSymbol());
//        Assert.assertEquals("F", cloneFormula.getIsotope(1).getSymbol());
//        Assert.assertEquals("H", cloneFormula.getIsotope(2).getSymbol());
    }
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetIsotopeCount_IIsotope_Occurr() {
        IMolecularFormula mf = getBuilder().newMolecularFormula();
        
        IIsotope carb = getBuilder().newIsotope("C");
        IIsotope flu = getBuilder().newIsotope("F");
        IIsotope h1 = getBuilder().newIsotope("H");
        mf.addIsotope( carb );
        mf.addIsotope( flu );
        mf.addIsotope( h1,3);
        
        Assert.assertEquals(3, mf.getIsotopeCount());
        Assert.assertEquals(1, mf.getIsotopeCount(carb));
        Assert.assertEquals(1, mf.getIsotopeCount(flu));
        Assert.assertEquals(3, mf.getIsotopeCount(h1));
    }
    
    
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testAdd_IMolecularFormula() {

        IMolecularFormula acetone = getBuilder().newMolecularFormula();
        acetone.addIsotope(getBuilder().newIsotope("C"),3);
        IIsotope oxig = getBuilder().newIsotope("O");
        acetone.addIsotope(oxig);
        
        
        Assert.assertEquals(2, acetone.getIsotopeCount());
        
        IMolecularFormula water = getBuilder().newMolecularFormula();
        water.addIsotope(getBuilder().newIsotope("H"),2);
        water.addIsotope(oxig);
        acetone.add(water);
        
        Assert.assertEquals(3, acetone.getIsotopeCount());
        
    }
    
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testIsotopes() {

        IMolecularFormula mf = getBuilder().newMolecularFormula();
        mf.addIsotope( getBuilder().newIsotope("C") );
        mf.addIsotope( getBuilder().newIsotope("F") );
        mf.addIsotope( getBuilder().newIsotope("H") ,3);
        
        Iterator<IIsotope> istoIter = mf.isotopes().iterator();
        int counter = 0;
        while (istoIter.hasNext()) {
        	istoIter.next();
            counter++;
        }
        Assert.assertEquals(3, counter);
    }
    
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testContains_IIsotope() {
    	IMolecularFormula mf = getBuilder().newMolecularFormula();
        
        IIsotope carb = getBuilder().newIsotope("C");
        IIsotope h1 = getBuilder().newIsotope("H");
        IIsotope h2 = getBuilder().newIsotope("H");
        h2.setExactMass(2.0004);
        
        mf.addIsotope( carb );
        mf.addIsotope( h1 );
    	
        Assert.assertTrue(mf.contains(carb));
        Assert.assertTrue(mf.contains(h1));
        Assert.assertFalse(mf.contains(h2));
    }

    
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testInstance_IIsotope() {

        IMolecularFormula mf = getBuilder().newMolecularFormula();
        
        IIsotope carb = getBuilder().newIsotope("C");
        IIsotope flu = getBuilder().newIsotope("F");
        IIsotope h1 = getBuilder().newIsotope("H");
        mf.addIsotope( carb );
        mf.addIsotope( flu );
        mf.addIsotope( h1 ,3);

        Iterator<IIsotope> istoIter = mf.isotopes().iterator();
        Assert.assertNotNull(istoIter);
        Assert.assertTrue(istoIter.hasNext());
        IIsotope next = istoIter.next();
        Assert.assertTrue(next instanceof IIsotope);
//        Assert.assertEquals(carb, next);
        
        Assert.assertTrue(istoIter.hasNext());
        next = istoIter.next();
        Assert.assertTrue(next instanceof IIsotope);
//        Assert.assertEquals(flu, next);
        
        Assert.assertTrue(istoIter.hasNext());
        next = istoIter.next();
        Assert.assertTrue(next instanceof IIsotope);
//        Assert.assertEquals(h1, next);
        
        Assert.assertFalse(istoIter.hasNext());
    }
    
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testGetCharge() {

        IMolecularFormula mf = getBuilder().newMolecularFormula();
        mf.setCharge(1);
        mf.addIsotope( getBuilder().newAtom("C") );
        mf.addIsotope( getBuilder().newAtom("F") );
        mf.addIsotope( getBuilder().newAtom("H"),3 );
        
        Assert.assertEquals(3, mf.getIsotopeCount());
        Assert.assertEquals(1.0,mf.getCharge(), 0.001);
        
    }
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testSetCharge_Double() {
    	

        IMolecularFormula mf = getBuilder().newMolecularFormula();
    	Assert.assertEquals(CDKConstants.UNSET, mf.getCharge());

    	
    	mf.setCharge(1);
        Assert.assertEquals(1.0, mf.getCharge(), 0.001);

        mf.add(mf);
        Assert.assertEquals(2.0, mf.getCharge(), 0.001);
    }
    
    @Test
    public void testSetCharge_Integer() {

        IMolecularFormula mf = getBuilder().newMolecularFormula();
        mf.setCharge(1);
        mf.addIsotope(getBuilder().newAtom("C"));
        mf.addIsotope(getBuilder().newAtom("F"));
        mf.addIsotope(getBuilder().newAtom("H"), 3);

        Assert.assertEquals(1.0, mf.getCharge(), 0.001);

    }

    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testCharge_rest() {
    	
        IMolecularFormula mf = getBuilder().newMolecularFormula();
    	Assert.assertEquals(CDKConstants.UNSET, mf.getCharge());
 	
    	mf.setCharge(1);
        Assert.assertEquals(1.0, mf.getCharge(), 0.001);

        IMolecularFormula mf2 = getBuilder().newMolecularFormula();
        mf2.setCharge(-1);
        mf.add(mf2);
        Assert.assertEquals(0.0, mf.getCharge(), 0.001);
    }
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testRemoveIsotope_IIsotope() {
    	
        IMolecularFormula mf = getBuilder().newMolecularFormula();
        IIsotope carb = getBuilder().newIsotope("C");
        IIsotope flu = getBuilder().newIsotope("F");
        IIsotope h1 = getBuilder().newIsotope("H");
        mf.addIsotope( carb );
        mf.addIsotope( flu );
        mf.addIsotope( h1,3);
        
        // remove the Fluorine 
        mf.removeIsotope(flu);
        
        Assert.assertEquals(2, mf.getIsotopeCount());
        
    }

    
    /**
	 * A unit test suite for JUnit.
	 *
	 * @return    The test suite
	 */
    @Test 
    public void testRemoveAllIsotopes() {
        IMolecularFormula mf = getBuilder().newMolecularFormula();
        IIsotope carb = getBuilder().newIsotope("C");
        IIsotope flu = getBuilder().newIsotope("F");
        IIsotope h1 = getBuilder().newIsotope("H");
        mf.addIsotope( carb );
        mf.addIsotope( flu );
        mf.addIsotope( h1,3);
        
        // remove the Fluorine 
        mf.removeAllIsotopes();
        
        Assert.assertEquals(0, mf.getIsotopeCount());
        
    }
    
    /**
	 * A unit test suite for JUnit. Only test whether the 
	 * MolecularFormula are correctly cloned.
	 *
	 * @return    The test suite
   	*/
    @Test 
	public void testClone() throws Exception {
        IMolecularFormula mf = getBuilder().newMolecularFormula();
        mf.setCharge(1);
        Object clone = mf.clone();
        Assert.assertTrue(clone instanceof IMolecularFormula);
        Assert.assertEquals(mf.getIsotopeCount(), ((IMolecularFormula) clone).getIsotopeCount());
        Assert.assertEquals(mf.getCharge(), ((IMolecularFormula) clone).getCharge());
        
	}   
	/**
	 * A unit test suite for JUnit. Only test whether 
	 * the MolecularFormula are correctly cloned.
   	*/
    @Test 
	public void testClone_Isotopes() throws Exception {
        IMolecularFormula mf = getBuilder().newMolecularFormula();
        IIsotope carb = getBuilder().newIsotope("C");
        IIsotope flu = getBuilder().newIsotope("F");
        IIsotope h1 = getBuilder().newIsotope("H");
        mf.addIsotope( carb );
        mf.addIsotope( flu );
        mf.addIsotope( h1 ,3);
        

        Assert.assertEquals(3, mf.getIsotopeCount());
        Assert.assertEquals(1, mf.getIsotopeCount(carb));
        Assert.assertEquals(1, mf.getIsotopeCount(flu));
        Assert.assertEquals(3, mf.getIsotopeCount(h1));
        
        Object clone = mf.clone();
        Assert.assertTrue(clone instanceof IMolecularFormula);
        Assert.assertEquals(mf.getIsotopeCount(), ((IMolecularFormula) clone).getIsotopeCount());

        Assert.assertEquals(3, ((IMolecularFormula) clone).getIsotopeCount());
	}  
    /**
	 * A unit test suite for JUnit.
   	*/
    @Test 
	public void testSetProperty_Object_Object() throws Exception {
        IMolecularFormula mf = getBuilder().newMolecularFormula();
        mf.setProperty("blabla", 2);
        Assert.assertNotNull(mf.getProperty("blabla"));
	} 
    /**
	 * A unit test suite for JUnit.
   	*/
    @Test 
	public void testRemoveProperty_Object() throws Exception {
    	 IMolecularFormula mf = getBuilder().newMolecularFormula();
    	 String blabla = "blabla";
    	 double number = 2;
         mf.setProperty(blabla,number);
         Assert.assertNotNull(mf.getProperty(blabla));
         
         mf.removeProperty("blabla");
         Assert.assertNull(mf.getProperty(blabla));
        
	} 
    /**
	 * A unit test suite for JUnit.
   	*/
    @Test 
	public void testGetProperty_Object() throws Exception {
    	testSetProperty_Object_Object();
        
	} 
    /**
	 * A unit test suite for JUnit.
   	*/
    @Test 
	public void testGetProperties() throws Exception {
    	 IMolecularFormula mf = getBuilder().newMolecularFormula();
         mf.setProperty("blabla", 2);
         mf.setProperty("blabla3", 3);
         Assert.assertEquals(2,mf.getProperties().size());
	} 
    /**
	 * A unit test suite for JUnit.
   	*/
    @Test 
	public void testSetProperties_Map() throws Exception {
    	testGetProperties();
        
	} 

    @Test public void testGetBuilder() {
    	IMolecularFormula add = getBuilder().newMolecularFormula();
    	IChemObjectBuilder builder = add.getBuilder();
    	Assert.assertNotNull(getBuilder());
    	Assert.assertEquals(getBuilder().getClass().getName(), builder.getClass().getName());
    }
}
