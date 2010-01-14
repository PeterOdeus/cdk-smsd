/**
 *
 * Copyright (C) 2006-2009  Syed Asad Rahman {asad@ebi.ac.uk}
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
package org.openscience.cdk.smsd.global;

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;

/**
 * @cdk.module smsd
 */
@TestClass("org.openscience.cdk.smsd.global")
public class TimeOut {

    private static TimeOut instance = null;
    private double time = -1;

    /**
     *
     * @return
     */
    @TestMethod("testGetInstance")
    public static synchronized TimeOut getInstance() {
        if (instance == null) {
            // it's ok, we can call this constructor
            instance = new TimeOut();
        }
        return instance;
    }

    protected TimeOut() {
    }

    /**
     *
     * @param timeout
     */
    @TestMethod("testSetTimeOut")
    public void setTimeOut(double timeout) {
        this.time = timeout;
    }

    /**
     *
     * @return
     */
    @TestMethod("testSetTimeOut")
    public double getTimeOut() {
        return time;
    }
}
