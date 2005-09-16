/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 * 
 * Copyright (C) 2004-2005  The Chemistry Development Kit (CDK) project
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
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.openscience.cdk.isomorphism.matchers.smarts;

import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.isomorphism.matchers.QueryAtom;

/**
 * @cdk.module extra
 */
public class AromaticQueryBond extends SMARTSBond {

    public AromaticQueryBond() {
    }

    public AromaticQueryBond(QueryAtom atom1, QueryAtom atom2, double order) {
        super(atom1, atom2, order);
    }
    
	public boolean matches(Bond bond) {
        if (bond.getFlag(CDKConstants.ISAROMATIC)) {
            return true;
        }
        return false;
    };

    public String toString() {
		StringBuffer s = new StringBuffer();
		s.append("AromaticQueryBond(");
        s.append(this.hashCode() + ", ");
        org.openscience.cdk.interfaces.Atom[] atoms = getAtoms();
		s.append("#A:" + atoms.length);
		for (int i = 0; i < atoms.length; i++) {
			if (atoms[i] == null) {
				s.append(", null");
			} else {
				s.append(", " + atoms[i].toString());
			}
		}
		s.append(")");
		return s.toString();
    }
}

