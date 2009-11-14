/* Copyright (C) 2006-2009  Syed Asad Rahman {asad@ebi.ac.uk}
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
package org.openscience.cdk.smsd.algorithm.mcsplus;

import java.io.IOException;
import java.util.List;
import java.util.Stack;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.smsd.algorithm.mcgregor.McGregor;
import org.openscience.cdk.smsd.core.tools.EBITimeManager;
import org.openscience.cdk.smsd.filters.ExactMapping;
import org.openscience.cdk.smsd.global.TimeOut;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @cdk.module smsd
 */
public class MCSPlus {

    public static double timeout = TimeOut.getInstance().getTimeOut();
    private static EBITimeManager TM;
    private static boolean timeoutFlag = false;

    /**
     * 
     * @param ac1
     * @param ac2
     * @param removeHydrogen
     * @return
     * @throws CDKException
     */
    public List<List<Integer>> getOverlaps(IAtomContainer ac1, IAtomContainer ac2, boolean removeHydrogen) throws CDKException {
        Stack<List<Integer>> Max_Cliques_Set = null;
        List<List<Integer>> _mappings = new Vector<List<Integer>>();


        try {

            GenerateCompatibilityGraph gcg = new GenerateCompatibilityGraph(ac1, ac2, removeHydrogen);
            List<Integer> comp_graph_nodes = gcg.get_comp_graph_nodes();

            List<Integer> C_edges = gcg.get_C_egdes();
            List<Integer> D_edges = gcg.get_D_egdes();

            if (D_edges.size() > 99999 && C_edges.size() > 2000) {
                System.err.println("D-edges Size " + D_edges.size() + " > : " + 99999);
                return null;
            }
            TM = new EBITimeManager();
            if (timeout > -1 && TM.getElapsedTimeInMinutes() > timeout) {
                //System.out.println("|Hello|");
                timeoutFlag = true;
                throw new CDKException("Timeout exceeded in getOverlaps");
            }
            BKKCKCF init = new BKKCKCF(comp_graph_nodes, C_edges, D_edges);
            Max_Cliques_Set = init.getMaxCliqueSet(); // build the CDKRGraph corresponding to this problem

            if (timeout > -1 && TM.getElapsedTimeInMinutes() > timeout) {
                //System.out.println("|Hello|");
                timeoutFlag = true;
                throw new CDKException("Timeout exceeded in getOverlaps");
            }
            //Clear all the compatibility graph content
            gcg.Clear();
            int clique_number = 1;
            while (!Max_Cliques_Set.empty()) {

                // build the CDKRGraph corresponding to this problem

                if (timeout > -1 && TM.getElapsedTimeInMinutes() > timeout) {
                    //System.out.println("|Hello|");
                    timeoutFlag = true;
                    throw new CDKException("Timeout exceeded in getOverlaps");
                }
                List<Integer> clique_vector = Max_Cliques_Set.peek();

                int clique_size = clique_vector.size();
                if (clique_size < ac1.getAtomCount() && clique_size < ac2.getAtomCount()) {
                    McGregor mgit = new McGregor(ac1, ac2, _mappings);
                    mgit.McGregor_IterationStart(mgit.getMCSSize(), clique_vector, comp_graph_nodes); //Start McGregorMCSPlus search

                    _mappings = mgit.getMappings();
                    mgit = null;

                } else {
                    _mappings = ExactMapping.extract_mapping(_mappings, comp_graph_nodes, clique_vector); //erfolgt KEINE McGregorMCSPlus-Suche->speichere L�ungen in final_MAPPINGS

                }
                Max_Cliques_Set.pop();
                clique_number++;


            }
        } catch (IOException ex) {
            Logger.getLogger(MCSPlus.class.getName()).log(Level.SEVERE, null, ex);
        }

        return _mappings;
    }

    /**
     *
     * @return
     */
    public static boolean getTimeOutFlag() {

        boolean flag = false;

        if (timeoutFlag) {
            flag = true;
        }
        System.out.println("Flag: " + flag);
        return flag;
    }
}
