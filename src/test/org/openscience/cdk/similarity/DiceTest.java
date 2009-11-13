/*
 * $RCSfile$    
 * $Author$    
 * $Date$    
 * $Revision$
 * 
 * Copyright (C) 1997-2007  The Chemistry Development Kit (CKD) project
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All I ask is that proper credit is given for my work, which includes
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
 * 
 */

package org.openscience.cdk.similarity;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.templates.MoleculeFactory;

import java.util.BitSet;

/**
 * @cdk.module test-fingerprint
 */
public class DiceTest extends CDKTestCase {

    boolean standAlone = false;

    @Test
    public void testDice1() throws Exception {
        Molecule mol1 = MoleculeFactory.makeIndole();
        Molecule mol2 = MoleculeFactory.makePyrrole();
        Fingerprinter fingerprinter = new Fingerprinter();
        BitSet bs1 = fingerprinter.getFingerprint(mol1);
        BitSet bs2 = fingerprinter.getFingerprint(mol2);
        float dice = Dice.calculate(bs1, bs2);
        if (standAlone) System.out.println("Dice: " + dice);
        if (!standAlone) Assert.assertEquals(0.5245, dice, 0.01);
    }

    @Test
    public void testDice2() throws Exception {
        Molecule mol1 = MoleculeFactory.makeIndole();
        Molecule mol2 = MoleculeFactory.makeIndole();
        Fingerprinter fingerprinter = new Fingerprinter();
        BitSet bs1 = fingerprinter.getFingerprint(mol1);
        BitSet bs2 = fingerprinter.getFingerprint(mol2);
        float dice = Dice.calculate(bs1, bs2);
        if (standAlone) System.out.println("Dice: " + dice);
        if (!standAlone) Assert.assertEquals(1.0, dice, 0.001);
    }
}