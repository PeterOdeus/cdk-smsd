/* $Revision: 5921 $ $Author: egonw $ $Date: 2006-04-12 11:16:35 +0200 (Wed, 12 Apr 2006) $    
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
package org.openscience.cdk.nonotify;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectChangeEvent;
import org.openscience.cdk.interfaces.IChemObjectListener;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.AtomContainerTest;

/**
 * Checks the functionality of the AtomContainer.
 *
 * @cdk.module test-nonotify
 */
public class NNAtomContainerTest extends AtomContainerTest {

    @BeforeClass public static void setUp() {
    	AtomContainerTest.builder = NoNotificationChemObjectBuilder.getInstance();
    }

    @Test public void testStateChanged_IChemObjectChangeEvent() {
        ChemObjectListenerImpl listener = new ChemObjectListenerImpl();
        IAtomContainer chemObject = builder.newAtomContainer();
        chemObject.addListener(listener);
        
        chemObject.addAtom(builder.newAtom());
        Assert.assertFalse(listener.changed);
        
        listener.reset();
        Assert.assertFalse(listener.changed);
        chemObject.addBond(builder.newBond(builder.newAtom(), builder.newAtom()));
        Assert.assertFalse(listener.changed);
    }

    private class ChemObjectListenerImpl implements IChemObjectListener {
        private boolean changed;
        
        private ChemObjectListenerImpl() {
            changed = false;
        }
        
        public void stateChanged(IChemObjectChangeEvent e) {
            changed = true;
        }
        
        public void reset() {
            changed = false;
        }
    }
    
}