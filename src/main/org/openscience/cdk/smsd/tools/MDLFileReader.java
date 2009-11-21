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
 * This program is distributed inputStream the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.smsd.tools;

import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.ChemModel;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.IChemObjectReader.Mode;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.tools.manipulator.MoleculeSetManipulator;

/**
 * @cdk.module smsd
 */
public class MDLFileReader {

    private static IMolecule Mol = null;

    /**
     *
     * @param inputStream
     * @param mode
     * @throws IOException
     */
    public MDLFileReader(InputStream inputStream, Mode mode) throws IOException {

        try {
            MDLV2000Reader reader2 = new MDLV2000Reader(inputStream, mode);
            Mol = (IMolecule) reader2.read(new Molecule());
            reader2.close();
        } catch (CDKException e) {
            String string = e.toString();
            if (string.contains("This file must be read with the MDLV3000Reader.")) {
                MDLV3000Reader reader2 = new MDLV3000Reader(inputStream, mode);
                try {
                    Mol = (IMolecule) reader2.read(new Molecule());
                } catch (CDKException ex) {
                    Logger.getLogger(MDLFileReader.class.getName()).log(Level.SEVERE, null, ex);
                }
                reader2.close();

            } else if (string.contains("This file must be read with the MDLReader.")) {
                try {
                    MDLReader reader2 = new MDLReader(inputStream, mode);
                    Mol = (IMolecule) reader2.read(new Molecule());
                    reader2.close();
                } catch (CDKException ex) {
                    Logger.getLogger(MDLFileReader.class.getName()).log(Level.SEVERE, null, ex);
                }


            }


        }

    }

    /**
     * 
     * @param inputStream
     * @throws IOException
     * @throws CDKException
     */
    public MDLFileReader(InputStream inputStream) throws IOException, CDKException {
        this(inputStream, Mode.RELAXED);
    }

    /**
     *
     * @param reader
     * @throws IOException
     * @throws CDKException
     */
    public MDLFileReader(Reader reader) throws IOException, CDKException {
        this(reader, Mode.RELAXED);
    }

    /**
     *
     * @param reader
     * @param mode
     * @throws IOException
     */
    public MDLFileReader(Reader reader, Mode mode) throws IOException {

        try {
            MDLV2000Reader reader2 = new MDLV2000Reader(reader, mode);
            Mol = (IMolecule) reader2.read(new Molecule());
            reader2.close();
        } catch (CDKException e) {
            String string = e.toString();
            if (string.contains("This file must be read with the MDLV3000Reader.")) {
                MDLV3000Reader reader2 = new MDLV3000Reader(reader, mode);
                try {
                    Mol = (IMolecule) reader2.read(new Molecule());
                } catch (CDKException ex) {
                    Logger.getLogger(MDLFileReader.class.getName()).log(Level.SEVERE, null, ex);
                }
                reader2.close();

            } else if (string.contains("This file must be read with the MDLReader.")) {
                try {
                    MDLReader reader2 = new MDLReader(reader, mode);
                    Mol = (IMolecule) reader2.read(new Molecule());
                    reader2.close();
                } catch (CDKException ex) {
                    Logger.getLogger(MDLFileReader.class.getName()).log(Level.SEVERE, null, ex);
                }


            }


        }

    }

    /**
     *
     * @return
     */
    public IMolecule getMolecule() {
        return Mol;
    }

    /**
     *
     * @return
     */
    public IMolecule getMoleculeWithLayoutCheck() {
        if (GeometryTools.has2DCoordinatesNew(Mol) != 2) {
            try {
                StructureDiagramGenerator sdg = new StructureDiagramGenerator(Mol);
                sdg.generateCoordinates();
                Mol = sdg.getMolecule();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        return Mol;
    }

    /**
     *
     * @return
     */
    public IChemModel getChemModelWithMoleculeWithLayoutCheck() {
        IChemModel chemModel = new ChemModel();
        chemModel.setMoleculeSet(ConnectivityChecker.partitionIntoMolecules(Mol));
        for (IAtomContainer molecule : MoleculeSetManipulator.getAllAtomContainers(chemModel.getMoleculeSet())) {
            if (GeometryTools.has2DCoordinatesNew(molecule) != 2) {
                try {
                    StructureDiagramGenerator sdg = new StructureDiagramGenerator((IMolecule) molecule);
                    sdg.generateCoordinates();
                    molecule = sdg.getMolecule();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        return chemModel;
    }
}


