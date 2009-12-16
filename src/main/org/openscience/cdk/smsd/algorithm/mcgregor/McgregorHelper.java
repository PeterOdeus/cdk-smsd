/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.algorithm.mcgregor;

import java.util.List;

/**
 *
 * @author Asad
 */
public class McgregorHelper {

    private final List<String> c_bond_setA;
    private final List<String> c_bond_setB;
    private final boolean mappingCheckFlag;
    private final int mappedAtomCount;
    private final List<Integer> mappedAtomsOrg;
    private final int neighborBondNumA;
    private final int neighborBondNumB;
    private final List<Integer> iBondNeighborAtomsA;
    private final List<Integer> iBondNeighborAtomsB;
    private final List<String> cBondNeighborsA;
    private final List<String> cBondNeighborsB;
    private final int setNumA;
    private final int setNumB;
    private final List<Integer> i_bond_setA;
    private final List<Integer> i_bond_setB;

    public McgregorHelper(boolean mappingCheckFlag,
            int mappedAtomCount,
            List<Integer> mappedAtomsOrg,
            int neighborBondNumA,
            int neighborBondNumB,
            List<Integer> iBondNeighborAtomsA,
            List<Integer> iBondNeighborAtomsB,
            List<String> cBondNeighborsA,
            List<String> cBondNeighborsB,
            int setNumA,
            int setNumB,
            List<Integer> i_bond_setA,
            List<Integer> i_bond_setB,
            List<String> c_bond_setA,
            List<String> c_bond_setB) {
        this.c_bond_setA = c_bond_setA;
        this.c_bond_setB = c_bond_setB;
        this.mappingCheckFlag = mappingCheckFlag;
        this.mappedAtomCount = mappedAtomCount;
        this.mappedAtomsOrg = mappedAtomsOrg;
        this.neighborBondNumA = neighborBondNumA;
        this.neighborBondNumB = neighborBondNumB;
        this.iBondNeighborAtomsA = iBondNeighborAtomsA;
        this.iBondNeighborAtomsB = iBondNeighborAtomsB;
        this.cBondNeighborsA = cBondNeighborsA;
        this.cBondNeighborsB = cBondNeighborsB;
        this.setNumA = setNumA;
        this.setNumB = setNumB;
        this.i_bond_setA = i_bond_setA;
        this.i_bond_setB = i_bond_setB;

    }

    /**
     * @return the c_bond_setA
     */
    public List<String> getC_bond_setA() {
        return c_bond_setA;
    }

    /**
     * @return the c_bond_setB
     */
    public List<String> getC_bond_setB() {
        return c_bond_setB;
    }

    /**
     * @return the mappingCheckFlag
     */
    public boolean isMappingCheckFlag() {
        return mappingCheckFlag;
    }

    /**
     * @return the mappedAtomCount
     */
    public int getMappedAtomCount() {
        return mappedAtomCount;
    }

    /**
     * @return the mappedAtomsOrg
     */
    public List<Integer> getMappedAtomsOrg() {
        return mappedAtomsOrg;
    }

    /**
     * @return the neighborBondNumA
     */
    public int getNeighborBondNumA() {
        return neighborBondNumA;
    }

    /**
     * @return the neighborBondNumB
     */
    public int getNeighborBondNumB() {
        return neighborBondNumB;
    }

    /**
     * @return the iBondNeighborAtomsA
     */
    public List<Integer> getiBondNeighborAtomsA() {
        return iBondNeighborAtomsA;
    }

    /**
     * @return the iBondNeighborAtomsB
     */
    public List<Integer> getiBondNeighborAtomsB() {
        return iBondNeighborAtomsB;
    }

    /**
     * @return the cBondNeighborsA
     */
    public List<String> getcBondNeighborsA() {
        return cBondNeighborsA;
    }

    /**
     * @return the cBondNeighborsB
     */
    public List<String> getcBondNeighborsB() {
        return cBondNeighborsB;
    }

    /**
     * @return the setNumA
     */
    public int getSetNumA() {
        return setNumA;
    }

    /**
     * @return the i_bond_setA
     */
    public List<Integer> getI_bond_setA() {
        return i_bond_setA;
    }

    /**
     * @return the i_bond_setB
     */
    public List<Integer> getI_bond_setB() {
        return i_bond_setB;
    }

    int getsetNumB() {
        return setNumB;
    }
}
