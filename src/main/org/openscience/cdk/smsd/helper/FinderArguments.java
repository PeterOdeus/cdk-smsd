/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.cdk.smsd.helper;

import java.util.ArrayList;
import java.util.List;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author      Syed Asad Rahman {asad@ebi.ac.uk}
 */

public class FinderArguments {

    //define the size of the fingerprint
    //NOTE: this should be a multiple of 64 and preferably not 1024 or 2048
    //as for these values we often get the random numbers for one-atom or
    //two-atom paths the same!
    protected boolean applyHAdding = false;
    protected boolean applyHRemoval = false;
    protected boolean applyTest = false;
    protected boolean applySuffix = false;
    protected boolean appendMode = false;
    protected boolean matchBondType = false;
    protected static final String matchFile = "mcs";
    protected static final String fingerFile = "finger";
    protected static final String graphFile = "graph";
    protected static final String descriptorFile = "molDescriptors";
    protected String suffix = null;
    protected List<String> fileName = new ArrayList<String>();

    public FinderArguments() {
        super();
    }

    protected void printError(String errmessg) {
        System.err.println(errmessg);
        System.exit(-1);
    }

    /**
     * Parses the options in the command line arguments and returns
     * an array of strings corresponding to the filenames given as arguments only
     * @param args
     * @return
     * @throws org.apache.commons.cli.ParseException 
     */
    protected String[] parseCommandLineOptions(String[] args) throws org.apache.commons.cli.ParseException {

        Options options = new Options();


        options.addOption("h", "help", false, "Help page for command usage");

        options.addOption("a", "addHydrogen", false, "Add Hydrogen");

        options.addOption("r", "removeHydrogen", false, "Remove Hydrogen");

        options.addOption("b", "bondType", false, "Match Bond types (Single, Double etc)");

        options.addOption("q", "Query", true, "Reactant/Query");

        options.addOption("t", "Target", true, "Product/Target");

        options.addOption("s", "addSuffix", true, "Add suffix to the files");

        options.addOption("m", "appendMode", false, "Appends the output to the exsisting files else" +
                " creates new files");


        CommandLine line = null;
        //try {
        CommandLineParser parser = new PosixParser();
        line = parser.parse(options, args);

        //System.out.println(line.getArgList().size());

        //} catch (ParseException exception) {
        //    System.err.println("Unexpected exception: " + exception.toString());
        //}


        if (line.hasOption("a")) {
            this.applyHAdding = true;
        }
        if (line.hasOption("r")) {
            this.applyHRemoval = true;
        }
        if (line.hasOption("s")) {
            this.applySuffix = true;
        }

        if (line.hasOption("b")) {
            this.matchBondType = true;
        }


        String[] files = line.getArgs();

        if (line.hasOption("h")) {
            System.out.println("Hello1");
            printHelp(options);
        }

        if (line.hasOption("s")) {
            //assert();
            String[] suffix_reader = line.getOptionValues("s");
            if (suffix_reader.length < 1) {
                System.out.println("Suffix required!");
                printHelp(options);
            }
            suffix = suffix_reader[0];
        //System.out.println(suffix);
        }



        if (!(line.hasOption("q") && line.hasOption("t"))) {
            //System.out.println("Hello");
            printHelp(options);
        }

        if (line.hasOption("q")) {

            String[] temp = line.getOptionValues("q");
            //System.out.println(temp[0].toString());
            fileName.add(temp[0]);
        }

        if (line.hasOption("t")) {
            //System.out.println(line.getOptionValue("t"));
            String[] temp = line.getOptionValues("t");
            fileName.add(temp[0]);
        }

        if (line.hasOption("m")) {
            this.appendMode = true;
        }

        //System.out.println(fileName);

        return files;
    }

    private void printHelp(Options options) {
        HelpFormatter formatter = new HelpFormatter();
        /*options.addOption("a", false, "Add Hydrogen");
        options.addOption("r", false, "Remove Hydrogen");
        options.addOption("s", false, "Add suffix to the files");
        options.addOption("t", false, "test case, presently disabled");*/

        options.addOption("h", "help", false, "Help page for command usage");

        options.addOption("a", "addHydrogen", false, "Add explicit Hydrogens where " +
                "\nthey are missing (default = false)");

        options.addOption("r", "removeHydrogen", false, "Remove all explicit " +
                "\nHydrogens (default=false)");

        options.addOption("b", "bondType", false, "Match Bond types (Single, Double etc)");

        options.addOption("q", "Query", true, "Reactant/Query");

        options.addOption("t", "Target", true, "Product/Target");

        options.addOption("s", "addSuffix", true, "Add suffix to the files (optional)");

        options.addOption("m", "appendMode", false, "Appends the output to the exsisting files else" +
                " creates new files");

        System.out.println("\n++++++++++++++++++++++++++++++++++++++++++++++" +
                "\nSMSD (Small Molecule Similarity Detector)" +
                "\n++++++++++++++++++++++++++++++++++++++++++++++" +
                "\nContact: Syed Asad Rahman," +
                "\n\t EMBL-EBI, Hinxton " +
                "\n\t Cambridge CB10 1SD" +
                "\n\t United Kingdom " +
                "\ne-mail: asad@ebi.ac.uk" +
                "\n++++++++++++++++++++++++++++++++++++++++++++++\n" +
                "\nSMSD software can calculate the similarity between" +
                "\ntwo small molecules by using an inhouse algorithm" +
                "\ndeveloped at EMBL-EBI. " +
                "It also uses CDK based" +
                "\nFingerprints and Maximum Common Subgraph (MCS) " +
                "\nmatching for finding similarity. " +
                "It create four" +
                "\nfiles which contains similarity scores between" +
                "\nmatched atoms.\n" +
                "\n++++++++++++++++++++++++++++++++++++++++++++++" +
                "\nNOTE: The graph matching is performed by removing" +
                "\nthe Hydrogens" +
                "\n++++++++++++++++++++++++++++++++++++++++++++++\n");

        formatter.printHelp("\n", options);
        System.out.println("\n++++++++++++++++++++++++++++++++++++++++++++++\n");

        System.exit(0);
    }
}



