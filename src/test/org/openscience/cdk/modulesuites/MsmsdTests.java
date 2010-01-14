/* Copyright (C) 2009  Egon Willighagen <egonw@user.sf.net>
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
package org.openscience.cdk.modulesuites;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;
import org.junit.runners.Suite.SuiteClasses;
import org.openscience.cdk.coverage.SmsdCoverageTest;
import org.openscience.cdk.smsd.SMSDTest;
import org.openscience.cdk.smsd.algorithm.cdk.CDKMCSHandlerTest;
import org.openscience.cdk.smsd.algorithm.cdk.CDKMCSTest;
import org.openscience.cdk.smsd.algorithm.mcsplus.MCSPlusHandlerTest;
import org.openscience.cdk.smsd.algorithm.single.SingleMappingHandlerTest;
import org.openscience.cdk.smsd.algorithm.vflib.VFLibTest;
import org.openscience.cdk.smsd.algorithm.vflib.VFlibMCSHandlerTest;
import org.openscience.cdk.smsd.factory.FragmentMatcherTest;
import org.openscience.cdk.smsd.factory.SubStructureSearchAlgorithmsTest;
import org.openscience.cdk.smsd.global.TimeOutTest;
import org.openscience.cdk.smsd.helper.LabelContainerTest;
import org.openscience.cdk.smsd.tools.BondEnergiesTest;
import org.openscience.cdk.smsd.tools.TimeManagerTest;

/**
 * TestSuite that runs all the unit tests for the smsd module.
 *
 * @cdk.module test-smsd
 */
@RunWith(value=Suite.class)
@SuiteClasses(value={
    SmsdCoverageTest.class,

    BondEnergiesTest.class,
    CDKMCSTest.class,
    SMSDTest.class,
    TimeManagerTest.class,
    VFLibTest.class,
    SubStructureSearchAlgorithmsTest.class,
    FragmentMatcherTest.class,
    TimeOutTest.class,
    LabelContainerTest.class,
    MCSPlusHandlerTest.class,
    VFlibMCSHandlerTest.class,
    SingleMappingHandlerTest.class,
    CDKMCSHandlerTest.class
})
public class MsmsdTests {}
