
/*
 * MappingHandler.java
 *
 * Created on March 1, 2006, 7:57 PM
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
package org.openscience.cdk.smsd.algorithm.cdk;

//~--- classes ----------------------------------------------------------------
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import org.openscience.cdk.smsd.core.tools.EBIException;
import org.openscience.cdk.smsd.helper.FinalMappings;
import org.openscience.cdk.smsd.helper.MolHandler;
import org.openscience.cdk.smsd.interfaces.IMCS;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMoleculeSet;

/**
 *
 * @Copyright (C)   2009  Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class CDKMCSHandler implements IMCS {

//    //~--- fields -------------------------------------------------------------
    private IAtomContainer Reactant;
    private IAtomContainer Product;
    private boolean RonPFlag = false;
    private static List<Map<IAtom, IAtom>> allAtomMCS = null;
    private static Map<IAtom, IAtom> atomsMCS = null;
    private static TreeMap<Integer, Integer> firstMCS = null;
    private static List<TreeMap<Integer, Integer>> allMCS = null;

    //~--- constructors -------------------------------------------------------
    /** Creates a new instance of MappingHandler
     */
    public CDKMCSHandler() {

        allAtomMCS = new Vector<Map<IAtom, IAtom>>();
        atomsMCS = new HashMap<IAtom, IAtom>();
        firstMCS = new TreeMap<Integer, Integer>();
        allMCS = new Vector<TreeMap<Integer, Integer>>();


    }

    private int checkForH(IAtomContainer mol) {
        int hCount = 0;


        for (int i = 0; i < mol.getAtomCount(); i++) {


            if (mol.getAtom(i).getSymbol().equals("H")) {
                hCount++;
            }

        }


        return hCount;
    }

    /**
     * Set the JMCS software
     * 
     * @param reactant
     * @param product
     */
    @Override
    public void set(IAtomContainer reactant, IAtomContainer product) {

        this.Reactant = reactant;
        this.Product = product;

    }

    /**
     * Set the JMCS software
     * 
     * @param reactant
     * @param product
     */
    @Override
    public void set(MolHandler reactant, MolHandler product) {


        this.Reactant = reactant.getMolecule();
        this.Product = product.getMolecule();

    }

    /**
     * Creates a new instance of SearchCliques
     * @param ReactantMolFileName
     * @param ProductMolFileName
     */
    @Override
    public void set(String ReactantMolFileName, String ProductMolFileName) {


        String mol1 = ReactantMolFileName;
        String mol2 = ProductMolFileName;

        this.Reactant = new MolHandler(mol1, false).getMolecule();
        this.Product = new MolHandler(mol2, false).getMolecule();


    }

    /**
     * 
     * @return
     */
    @Override
    public int search_MCS(boolean removeHydrogen) throws EBIException {


        CDKRMapHandler rmap = new CDKRMapHandler();

        try {

            if ((Reactant.getAtomCount() == Product.getAtomCount()) && Reactant.getBondCount() == Product.getBondCount()) {
                RonPFlag = true;
                //System.err.print("True");

                rmap.calculateOverlapsAndReduceExactMatch(Reactant, Product);

            } else if (Reactant.getAtomCount() > Product.getAtomCount() && Reactant.getBondCount() != Product.getBondCount()) {
                RonPFlag = true;
                //System.err.print("True");

                rmap.calculateOverlapsAndReduce(Reactant, Product);

            } else {
                RonPFlag = false;
                rmap.calculateOverlapsAndReduce(Product, Reactant);


            }

            setAllMapping();
            setAllAtomMapping();
            setFirstMapping();
            setFirstAtomMapping();

        } catch (EBIException e) {
            rmap = null;
//            System.err.println("WARNING: CDKMCS: most probably time out error ");
        }

        return 0;
    }

    /**
     *
     * @param ac1
     * @param ac2
     * 
     */
    private IAtomContainer CalculateMCSS(IAtomContainer ac1, IAtomContainer ac2) throws EBIException {



        IAtomContainer MaxCommonSubgraphContainer = null;
        int maxac = -999;
        List overlaps = CDKMCS.getOverlaps(ac1, ac2);

        for (Object aSsList : overlaps) {
            IAtomContainer ss = (IAtomContainer) aSsList;
            if (ss.getAtomCount() > maxac) {
                MaxCommonSubgraphContainer = ss;
                maxac = ss.getAtomCount();
            }
        }

        return MaxCommonSubgraphContainer;

    }

    /**
     * 
     * @param mol
     * @param mcss
     * @return
     * @throws EBIException 
     */
    protected IMoleculeSet getUncommon(IAtomContainer mol, IAtomContainer mcss) throws EBIException {
        ArrayList<Integer> atomSerialsToDelete = new ArrayList<Integer>();

        List matches = CDKMCS.getSubgraphAtomsMaps(mol, mcss);
        List mapList = (List) matches.get(0);
        for (Object o : mapList) {
            CDKRMap rmap = (CDKRMap) o;
            atomSerialsToDelete.add(rmap.getId1());
        }

        // at this point we have the serial numbers of the bonds to delete
        // we should get the actual bonds rather than delete by serial numbers
        ArrayList<IAtom> atomsToDelete = new ArrayList<IAtom>();
        for (Integer serial : atomSerialsToDelete) {
            atomsToDelete.add(mol.getAtom(serial));
        }

        // now lets get rid of the bonds themselves
        for (IAtom atom : atomsToDelete) {
            mol.removeAtomAndConnectedElectronContainers(atom);
        }

        // now we probably have a set of disconnected components
        // so lets get a set of individual atom containers for
        // corresponding to each component
        return ConnectivityChecker.partitionIntoMolecules(mol);
    }

    /*public void setStereoScore() {
    int size = atomsMCS.size();
    for (int i = 0; i < size - 3; i += 4) {
    int RBondNumber = Reactant.getBondCount(atomsMCS.get(i), atomsMCS.get(i + 2));
    int PBondNumber = Product.getBondCount(atomsMCS.get(i + 1), atomsMCS.get(i + 3));
    if (RBondNumber > -1 && PBondNumber > -1) {
    IBond RBond = Reactant.getBond(RBondNumber);
    IBond PBond = Product.getBond(PBondNumber);
    if (RBond.getStereo() == PBond.getStereo()) {
    stereoScore++;
    }
    }
    }
    }*/
    private IAtomContainer getSubstructreBasedOnAtomUID(IAtomContainer RefMol, IAtomContainer Substructure) {

        IAtomContainer needle = new AtomContainer();

        Vector<IAtom> idlist = new Vector<IAtom>();

        // get the ID's (corresponding to the serial number of the Bond object in
        // the AtomContainer for the supplied molecule) of the matching bonds
        // (there will be repeats)

        //IAtomContainer IM=DefaultChemObjectBuilder.getInstance().newAtomContainer();

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

    //~--- get methods --------------------------------------------------------
    public final synchronized void setAllMapping() {

        //int count_final_sol = 1;
        //System.out.println("Output of the final FinalMappings: ");
        try {
            List<TreeMap<Integer, Integer>> sol = FinalMappings.getInstance().getFinalMapping();
            int counter = 0;
            for (TreeMap<Integer, Integer> final_solution : sol) {
                TreeMap<Integer, Integer> atomMappings = new TreeMap<Integer, Integer>();
                for (Map.Entry<Integer, Integer> Solutions : final_solution.entrySet()) {

                    int IIndex = Solutions.getKey().intValue();
                    int JIndex = Solutions.getValue().intValue();

                    if (RonPFlag) {
                        atomMappings.put(IIndex, JIndex);
                    } else {
                        atomMappings.put(JIndex, IIndex);
                    }
                }


                if (!allMCS.contains(atomMappings)) {
//                    System.out.println("atomMappings: " + atomMappings);
                    allMCS.add(counter++, atomMappings);
                }
            }

        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    public final synchronized void setAllAtomMapping() {
        //int count_final_sol = 1;
        //System.out.println("Output of the final FinalMappings: ");
        List<TreeMap<Integer, Integer>> sol = allMCS;

        int counter = 0;
        for (TreeMap<Integer, Integer> final_solution : sol) {

            Map<IAtom, IAtom> atomMappings = new HashMap<IAtom, IAtom>();


            for (Map.Entry<Integer, Integer> Solutions : final_solution.entrySet()) {

                int IIndex = Solutions.getKey().intValue();
                int JIndex = Solutions.getValue().intValue();

                IAtom A = null;

                IAtom B = null;

                A = Reactant.getAtom(IIndex);
                B = Product.getAtom(JIndex);
//                A.setID(Solutions.getKey().toString());
//                B.setID(Solutions.getValue().toString());
                atomMappings.put(A, B);


            }

            allAtomMCS.add(counter++, atomMappings);
        }



    }

    public final synchronized void setFirstMapping() {

        if (allMCS.size() > 0) {
            firstMCS = new TreeMap<Integer, Integer>(allMCS.get(0));
        }

    }

    public final synchronized void setFirstAtomMapping() {
        if (atomsMCS.size() > 0) {
            atomsMCS = new TreeMap<IAtom, IAtom>(allAtomMCS.get(0));
        }

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
        return atomsMCS;
    }
}
