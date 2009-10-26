/*
 * JMCSHandler.java
 *
 * Created on January 28, 2007, 11:37 AM
 *
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package org.openscience.cdk.smsd.helper;

//~--- JDK imports ------------------------------------------------------------
import org.openscience.cdk.smsd.core.tools.EBIAtomContainerManipulator;
import org.openscience.cdk.smsd.core.tools.MoleculeSanityCheck;
import java.io.FileInputStream;
import java.io.IOException;
//~--- non-JDK imports --------------------------------------------------------

import java.io.InputStreamReader;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.BondTools;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.MDLReader;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class MolHandler {

    private IAtomContainer Mol = null;
    private IAtomContainerSet FragmentMolSet = null;
    private boolean removeHydrogen = false;
    private boolean FragmentFlag = false;

    private void checkFragmentation() {

        if (Mol.getAtomCount() > 0) {
            FragmentFlag = ConnectivityChecker.isConnected(Mol);
        }
        FragmentMolSet = DefaultChemObjectBuilder.getInstance().newMoleculeSet();


        //System.out.println("isConnected : " + FragmentFlag);

        if (!FragmentFlag) {

            /*System.err.println("The molecule is not connected, " +
            "Fragment Matcher Will Handle this _molecule");*/

            FragmentMolSet.add(ConnectivityChecker.partitionIntoMolecules(Mol));
            FragmentMolSet.setID(Mol.getID());

        } else {

            FragmentMolSet.addAtomContainer(Mol);
            FragmentMolSet.setID(Mol.getID());
        }
    }

    /** Creates a new instance of JMCSHandler
     * @param MolFile Mol file name
     * @param cleanMolecule
     * @param removeHydrogen 
     *  
     */
    public MolHandler(String MolFile, boolean cleanMolecule, boolean removeHydrogen) {

        MDLReader MolRead;

        //Mol2Reader Mol2Read;
        this.removeHydrogen = removeHydrogen;


        try {


            FileInputStream ReadMolecule;

            ReadMolecule = new FileInputStream(MolFile);
            MolRead = new MDLReader(new InputStreamReader(ReadMolecule));
            this.Mol = (IMolecule) MolRead.read(new Molecule());


            if (cleanMolecule) {
                MoleculeSanityCheck.fixAromaticity((IMolecule) Mol);
            }


            BondTools.makeUpDownBonds(Mol);



//            /*Remove Hydrogen by Asad*/
            if (removeHydrogen) {
//                removeHydrogenAtoms(Mol);
                Mol = EBIAtomContainerManipulator.removeHydrogens(Mol);

            }

            checkFragmentation();



        } catch (IOException ex) {
            Logger.getLogger(MolHandler.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException e) {
            System.err.println(e);
        }
    }

    /**
     *
     * @param MolFile
     * @param cleanMolecule
     */
    public MolHandler(String MolFile, boolean cleanMolecule) {

        MDLReader MolRead;
        this.removeHydrogen = false;


        try {


            FileInputStream ReadMolecule;

            ReadMolecule = new FileInputStream(MolFile);
            MolRead = new MDLReader(new InputStreamReader(ReadMolecule));
            this.Mol = (IMolecule) MolRead.read(new Molecule());


            if (cleanMolecule) {
                MoleculeSanityCheck.fixAromaticity((IMolecule) Mol);
            }


            BondTools.makeUpDownBonds(Mol);


//            /*Remove Hydrogen by Asad*/
            if (removeHydrogen) {
//                removeHydrogenAtoms(Mol);
                Mol = EBIAtomContainerManipulator.removeHydrogens(Mol);

            }

            checkFragmentation();

        } catch (IOException ex) {
            Logger.getLogger(MolHandler.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException e) {
            System.err.println(e);
        }
    }

    /**
     * 
     * @param _molecule Molecule AtomContainer
     * @param cleanMolecule
     * @param removeHydrogen
     */
    public MolHandler(IAtomContainer _molecule, boolean cleanMolecule, boolean removeHydrogen) {
        try {
            CDKHueckelAromaticityDetector.detectAromaticity(_molecule);
        } catch (CDKException ex) {
            Logger.getLogger(MolHandler.class.getName()).log(Level.SEVERE, null, ex);
        }

        String ID = _molecule.getID();
        this.removeHydrogen = removeHydrogen;


        this.Mol = _molecule;


        if (cleanMolecule) {
            MoleculeSanityCheck.fixAromaticity((IMolecule) Mol);
        }
        /*Hydrogen are always removed for this molecule before mapping*/

        if (removeHydrogen) {
            try {
                this.Mol = (IMolecule) EBIAtomContainerManipulator.removeHydrogensAndPreserveAtomID(Mol);
                Mol.setID(ID);
            } catch (Exception ex) {
                Logger.getLogger(MolHandler.class.getName()).log(Level.SEVERE, null, ex);
            }

        } else {
            this.Mol = DefaultChemObjectBuilder.getInstance().newMolecule(Mol);
            Mol.setID(ID);

        }

        checkFragmentation();


    }

    /**
     *
     * @param _molecule
     * @param cleanMolecule
     */
    public MolHandler(IAtomContainer _molecule, boolean cleanMolecule) {

        String ID = _molecule.getID();
        this.removeHydrogen = false;
        this.Mol = _molecule;
        if (cleanMolecule) {
            MoleculeSanityCheck.fixAromaticity((IMolecule) Mol);
        }

//         /*Hydrogen are always removed for this molecule before mapping*/

        if (removeHydrogen) {
            try {
                this.Mol = (IMolecule) EBIAtomContainerManipulator.removeHydrogensAndPreserveAtomID(Mol);
                Mol.setID(ID);
            } catch (Exception ex) {
                Logger.getLogger(MolHandler.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else {
            this.Mol = DefaultChemObjectBuilder.getInstance().newMolecule(Mol);
            Mol.setID(ID);
        }

        checkFragmentation();


    }

    public IAtomContainer getMolecule() {

//        System.out.println("MCSMolHandle Mol:" + Mol.getID() + "  " + Mol.getAtomCount());
        return Mol;
    }

    public boolean getRemoveHydrogenFlag() {

        return removeHydrogen;
    }

    /**
     *
     * @return
     */
    public IAtomContainerSet getFragmentedMolecule() {

        return this.FragmentMolSet;
    }

    public boolean getFragmentFlag() {
        return this.FragmentFlag;
    }
}
