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
package org.openscience.cdk.smsd.helper;

import java.util.HashMap;
import java.util.Map;

/**
 * @cdk.module smsd
 */
public class LabelContainer {

    private HashMap<String, Integer> labelMap = null;
    private int labelCounter = 1;
    private static LabelContainer instance = null;

    protected LabelContainer() {

        // System.err.println("List Initialized");
        labelMap = new HashMap<String, Integer>();

        labelMap.put("C", labelCounter++);
        labelMap.put("O", labelCounter++);
        labelMap.put("N", labelCounter++);
        labelMap.put("S", labelCounter++);

        labelMap.put("P", labelCounter++);
        labelMap.put("F", labelCounter++);
        labelMap.put("I", labelCounter++);
        labelMap.put("R", labelCounter++);

        labelMap.put("Br", labelCounter++);
        labelMap.put("Cl", labelCounter++);
        labelMap.put("Co", labelCounter++);
        labelMap.put("Fe", labelCounter++);
        labelMap.put("Na", labelCounter++);
        labelMap.put("Ca", labelCounter++);

        labelMap.put("K", labelCounter++);
        labelMap.put("Mg", labelCounter++);
        labelMap.put("Se", labelCounter++);
        labelMap.put("Cu", labelCounter++);
        labelMap.put("Hg", labelCounter++);

        labelMap.put("X", labelCounter++);
        labelMap.put("R", labelCounter++);
        labelMap.put("X1", labelCounter++);


    }

    /**
     *
     * @return
     */
    synchronized public static LabelContainer getInstance() {
        if (instance == null) {
            instance = new LabelContainer();
        }
        return instance;
    }

    /**
     *
     * @param label
     */
    synchronized public void addLabel(String label) {
        //System.err.println("List added");
        int lableNo = labelMap.size() + 1;
        labelMap.put(label, lableNo);
    }

    /**
     *
     * @param label
     * @return
     */
    synchronized public Integer getLabelID(String label) {

        int labelID = -1;

        if (!labelMap.containsKey(label)) {
            int lableNo = labelMap.size() + 1;
            labelMap.put(label, lableNo);

        }

        labelID = labelMap.get(label);


        return labelID;
    }

    /**
     *
     * @param labelID
     * @return
     */
    synchronized public String getLabel(Integer labelID) {

        String indexLabel = null;
        boolean flag = false;

        for (Map.Entry<String, Integer> map : labelMap.entrySet()) {

            indexLabel = map.getKey();

            if (labelID.equals(map.getValue())) {
                flag = true;
                break;
            }
        }

        if (!flag) {
            indexLabel = null;
        }
        return indexLabel;
    }

    /**
     *
     * @return
     */
    synchronized public int getSize() {
        return labelMap.size();
    }
}
