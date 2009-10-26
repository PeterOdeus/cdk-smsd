/*
 * SearchCliques.java
 *
 * Created on January 17, 2007, 11:25 PM
 *
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package org.openscience.cdk.smsd.algorithm.mcsplus;


import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.smsd.core.tools.EBIException;
import org.openscience.cdk.smsd.filters.PostFilter;
import org.openscience.cdk.smsd.helper.FinalMappings;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.IMCS;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * 
 * @author Syed Asad Rahman 
 * @version March 2007
 * This is the main interface of the MCS algorithm. 
 */
public class MCSPlusHandler implements IMCS {

    private static List<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static List<TreeMap<Integer, Integer>> allMCS = null;
    private IAtomContainer ac1 = null;
    private IAtomContainer ac2 = null;
    private boolean flagExchange = false;

    public MCSPlusHandler() {


        allAtomMCS = new Vector<Map<IAtom, IAtom>>();
        atomsMCS = new HashMap<IAtom, IAtom>();
        firstMCS = new TreeMap<Integer, Integer>();
        allMCS = new Vector<TreeMap<Integer, Integer>>();

    }

    /** Creates a new instance of SearchCliques
     * @param Reactant
     * @param Product
     * @throws java.io.IOException
     *
     *
     */
    @Override
    public void set(MolHandler Reactant, MolHandler Product) throws IOException {

        this.ac1 = Reactant.getMolecule();
        this.ac2 = Product.getMolecule();


    }

    /** Creates a new instance of SearchCliques
     * @param ReactantMolFileName
     * @param ProductMolFileName
     * @throws java.io.IOException
     */
    @Override
    public void set(String ReactantMolFileName, String ProductMolFileName) throws IOException {


        String mol1 = ReactantMolFileName;
        String mol2 = ProductMolFileName;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);

        this.set(Reactant, Product);
        // System.exit(1);

    }

    /** Creates a new instance of SearchCliques
     * @param ReactantMol
     * @param ProductMol
     * @throws java.io.IOException
     */
    @Override
    public void set(IAtomContainer ReactantMol, IAtomContainer ProductMol) throws IOException {

        IAtomContainer mol1 = ReactantMol;
        IAtomContainer mol2 = ProductMol;

        MolHandler Reactant = new MolHandler(mol1, false);
        MolHandler Product = new MolHandler(mol2, false);

        this.set(Reactant, Product);
        // System.exit(1);

    }

    //Function is called by the main program and serves as a starting point for the comparision procedure.
    /**
     *
     * @param removeHydrogen
     * @return
     * @throws java.io.IOException
     */
    @Override
    public int search_MCS(boolean removeHydrogen) throws IOException {
        List<List<Integer>> _mappings = null;
        try {
            if (ac1.getAtomCount() > ac2.getAtomCount()) {
                _mappings = new MCSPlus().getOverlaps(ac1, ac2, removeHydrogen);
            } else {
                flagExchange = true;
                _mappings = new MCSPlus().getOverlaps(ac2, ac1, removeHydrogen);

            }

            PostFilter.filter(_mappings);
            setAllMapping();
            setAllAtomMapping();
            setFirstMapping();
            setFirstAtomMapping();
        } catch (EBIException e) {
            _mappings = null;
//            System.err.println("WARNING: MCSPlus: most probably time out error ");
        }

//        System.out.println("Mappings" + _mappings);

        return 0;
    }

    public final void setAllMapping() {

//        int count_final_sol = 1;
        //System.out.println("Output of the final FinalMappings: ");
        try {

            List<TreeMap<Integer, Integer>> final_solution = FinalMappings.getInstance().getFinalMapping();
            int counter = 0;
            for (TreeMap<Integer, Integer> solution : final_solution) {
//                System.out.println("Number of MCS solution: " + solution);
                TreeMap<Integer, Integer> validSolution = new TreeMap<Integer, Integer>();

                if (!flagExchange) {
                    for (Map.Entry<Integer, Integer> map : solution.entrySet()) {
                        validSolution.put(map.getKey(), map.getValue());
                    }
                } else {
                    for (Map.Entry<Integer, Integer> map : solution.entrySet()) {
                        validSolution.put(map.getValue(), map.getKey());
                    }
                }
                allMCS.add(counter++, validSolution);
            }

        } catch (Exception ex) {
            ex.printStackTrace();
        }

//        System.out.println("Number of MCS solution: " + allMCS.size());


    }

    public final synchronized void setAllAtomMapping() {

        try {
//            int count_final_sol = 1;
            //System.out.println("Output of the final FinalMappings: ");
            List<TreeMap<Integer, Integer>> final_solution = FinalMappings.getInstance().getFinalMapping();

            int counter = 0;
            for (TreeMap<Integer, Integer> solution : final_solution) {


                Map<IAtom, IAtom> atomMappings = new HashMap<IAtom, IAtom>();

//                System.out.println("Sol size " + solution.size());


                for (Map.Entry<Integer, Integer> map : solution.entrySet()) {

                    //  int IIndex = Solutions.get(i) -1;
                    //  int JIndex = Solutions.get(i + 1) -1;

                    int IIndex = map.getKey();
                    int JIndex = map.getValue();

                    //System.out.println("I " + IIndex + " J " + JIndex);

                    //IAtom A = ac1.getAtom(IIndex);
                    //IAtom B = ac2.getAtom(JIndex);

                    IAtom A = null;
                    IAtom B = null;

                    if (!flagExchange) {
                        A = ac1.getAtom(IIndex);
                        B = ac2.getAtom(JIndex);
                    } else {
                        A = ac2.getAtom(IIndex);
                        B = ac1.getAtom(JIndex);
                    }

                    atomMappings.put(A, B);

//                    System.out.println("I " + A.getSymbol() + " J " + B.getSymbol());
                }

                allAtomMCS.add(counter++, atomMappings);
            }
        } catch (Exception I) {
            I.getCause();
        }

//        System.out.println("Number of Atom MCS solution: " + allAtomMCS.size());


    }

    public final synchronized void setFirstMapping() {

        if (allMCS.size() > 0) {
            firstMCS = new TreeMap<Integer, Integer>(allMCS.get(0));
        }

    }

    public final synchronized void setFirstAtomMapping() {
        if (allAtomMCS.size() > 0) {
//            System.out.println("In MCS+handle First Atom MCS: " + allAtomMCS.get(0));
            atomsMCS = new HashMap<IAtom, IAtom>(allAtomMCS.get(0));
        }

    }

    private IAtomContainer getSubstructreBasedOnAtomUID(
            IAtomContainer RefMol, IAtomContainer Substructure) {

        IAtomContainer needle = new AtomContainer();

        Vector<IAtom> idlist = new Vector<IAtom>();

        // get the ID's (corresponding to the serial number of the Bond object in
        // the AtomContainer for the supplied molecule) of the matching bonds
        // (there will be repeats)

        for (int i = 0; i < RefMol.getAtomCount(); i++) {

            for (int j = 0; j < Substructure.getAtomCount(); j++) {

                if ((RefMol.getAtom(i).getID()).equals(Substructure.getAtom(j).getID())) {

                    idlist.add(RefMol.getAtom(i));

                }

            }

        }

        // get a unique list of bond ID's and add them to an AtomContainer
        HashSet<IAtom> hs = new HashSet<IAtom>(idlist);
        for (IAtom h : hs) {
            needle.addAtom(h);
        }

        return needle;

    }

    @Override
    public List<TreeMap<Integer, Integer>> getAllMapping() {
        return allMCS;
    }

    @Override
    public TreeMap<Integer, Integer> getFirstMapping() {

        return firstMCS;
    }

    @Override
    public List<Map<IAtom, IAtom>> getAllAtomMapping() {
        return allAtomMCS;
    }

    @Override
    public Map<IAtom, IAtom> getFirstAtomMapping() {
//        System.out.println("In MCS+handle First Atom MCS: " + atomsMCS);
        return atomsMCS;
    }
}
