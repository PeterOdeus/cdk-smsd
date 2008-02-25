/* $RCSfile$
 * $Author$
 * $Date$
 * $Revision$
 *
 * Copyright (C) 1997-2007  The JChemPaint project
 *
 * Contact: jchempaint-devel@lists.sourceforge.net
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
package org.openscience.cdk.applications.jchempaint.action;

import java.awt.event.ActionEvent;

import org.openscience.cdk.applications.jchempaint.JChemPaintEditorPanel;

/**
 * Opens a new empty JChemPaintFrame.
 *
 * @cdk.module jchempaint
 * @cdk.svnrev  $Revision$
 * @author     steinbeck
 * @cdk.created    1997
 */
public class ToolBarAction extends JCPAction {


    private static final long serialVersionUID = -7685519078261241187L;

    /**
     *  Opens an empty JChemPaint frame.
     *
     *@param  e  Description of the Parameter
     */
    public void actionPerformed(ActionEvent e) 
    {
	    if (((JChemPaintEditorPanel)jcpPanel).getShowToolBar()) ((JChemPaintEditorPanel)jcpPanel).setShowToolBar(false);
	    else ((JChemPaintEditorPanel)jcpPanel).setShowToolBar(true);
    }
}
