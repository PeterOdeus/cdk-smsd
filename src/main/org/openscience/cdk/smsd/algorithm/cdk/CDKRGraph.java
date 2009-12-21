
/* 
 *
 * Copyright (C) 2002-2007  Stephane Werner <mail@ixelis.net>
 *               2007-2009  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * This code has been kindly provided by Stephane Werner
 * and Thierry Hanser from IXELIS mail@ixelis.net.
 *
 * IXELIS sarl - Semantic Information Systems
 *               17 rue des C?dres 67200 Strasbourg, France
 *               Tel/Fax : +33(0)3 88 27 81 39 Email: mail@ixelis.net
 *
 * CDK Contact: cdk-devel@lists.sf.net
 *
 * This program is free software; you can redistribute maxIterator and/or
 * modify maxIterator under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that maxIterator will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR sourceBitSet PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.smsd.algorithm.cdk;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.smsd.tools.TimeManager;
import org.openscience.cdk.smsd.global.TimeOut;

/**
 * This class implements the Resolution Graph (CDKRGraph).
 * The CDKRGraph is a graph based representation of the search problem.
 * An CDKRGraph is constructed from the two compared graphs (G1 and G2).
 * Each vertex (node) in the CDKRGraph represents a possible association
 * from an edge in G1 with an edge in G2. Thus two compatible bonds
 * in two molecular graphs are represented by a vertex in the CDKRGraph.
 * Each edge in the CDKRGraph corresponds to a common adjacency relationship
 * between the 2 couple of compatible edges associated to the 2 CDKRGraph nodes
 * forming this edge.
 *
 * <p>Example:
 * <pre>
 *    G1 : C-C=O  and G2 : C-C-C=0
 *         1 2 3           1 2 3 4
 * </pre>
 *
 *  <p>The resulting CDKRGraph(G1,G2) will contain 3 nodes:
 *  <ul>
 *    <li>Node sourceBitSet : association between bond C-C :  1-2 in G1 and 1-2 in G2
 *    <li>Node targetBitSet : association between bond C-C :  1-2 in G1 and 2-3 in G2
 *    <li>Node C : association between bond C=0 :  2-3 in G1 and 3-4 in G2
 *  </ul>
 *  The CDKRGraph will also contain one edge representing the
 *  adjacency between node targetBitSet and C  that is : bonds 1-2 and 2-3 in G1
 *  and bonds 2-3 and 3-4 in G2.
 *
 *  <p>Once the CDKRGraph has been built from the two compared graphs
 *  maxIterator becomes a very interesting tool to perform all kinds of
 *  structural search (isomorphism, substructure search, maximal common
 *  substructure,....).
 *
 *  <p>The  search may be constrained by mandatory elements (e.g. bonds that
 *  have to be present in the mapped common substructures).
 *
 *  <p>Performing a query on an CDKRGraph requires simply to set the constrains
 *  (if any) and to invoke the parsing method (parse())
 *
 *  <p>The CDKRGraph has been designed to be a generic tool. It may be constructed
 *  from any kind of source graphs, thus maxIterator is not restricted to a chemical
 *  context.
 *
 *  <p>The CDKRGraph model is indendant from the CDK model and the link between
 *  both model is performed by the RTools class. In this way the CDKRGraph
 *  class may be reused in other graph context (conceptual graphs,....)
 *
 *  <p><bitSet>Important note</bitSet>: This implementation of the algorithm has not been
 *                      optimized for speed at this stage. It has been
 *                      written with the goal to clearly retrace the
 *                      principle of the underlined search method. There is
 *                      room for optimization in many ways including the
 *                      the algorithm itself.
 *
 *  <p>This algorithm derives from the algorithm described in
 *  {@cdk.cite HAN90} and modified in the thesis of T. Hanser {@cdk.cite HAN93}.
 */
/**
 * @cdk.module smsd
 */
public class CDKRGraph {
    // an CDKRGraph is a list of CDKRGraph nodes
    // each node keeping track of its
    // neighbors.

    List<CDKRNode> graph = null;
    // maximal number of iterations before
    // search break
    int maxIteration = -1;
    // dimensions of the compared graphs
    int firstGraphSize = 0;
    int secondGraphSize = 0;
    // constrains
    BitSet sourceBitSet = null;
    BitSet targetBitSet = null;
    // current solution list
    List<BitSet> solutionList = null;
    // flag to define if we want to get all possible 'mappings'
    boolean findAllMap = false;
    // flag to define if we want to get all possible 'structures'
    boolean findAllStructure = true;
    // working variables
    boolean stop = false;
    int nbIteration = 0;
    BitSet graphBitSet = null;
    private static double timeout = TimeOut.getInstance().getTimeOut(); // -1 for infinite search and one min is 1
    private TimeManager timeManager;

    /**
     * Constructor for the CDKRGraph object and creates an empty CDKRGraph.
     */
    public CDKRGraph() {
        graph = new ArrayList<CDKRNode>();
        solutionList = new ArrayList<BitSet>();
        graphBitSet = new BitSet();
    }

    private void checkTimeOut() throws CDKException {
        if (timeout > -1 && timeManager.getElapsedTimeInMinutes() > timeout) {
            throw new CDKException("Timeout exceeded in getOverlaps");
        }
    }

    /**
     *  Returns the size of the first of the two
     *  compared graphs.
     * @return The size of the first of the two compared graphs
     */
    public int getFirstGraphSize() {
        return firstGraphSize;
    }

    /**
     *  Returns the size of the second of the two
     *  compared graphs.
     * @return The size of the second of the two compared graphs
     */
    public int getSecondGraphSize() {
        return secondGraphSize;
    }

    /**
     *  Sets the size of the first of the two
     *  compared graphs.
     * @param graphSize The size of the second of the two compared graphs
     */
    public void setFirstGraphSize(int graphSize) {
        firstGraphSize = graphSize;
    }

    /**
     *  Returns the size of the second of the two
     *  compared graphs.
     * @param graphSize The size of the second of the two compared graphs
     */
    public void setSecondGraphSize(int graphSize) {
        secondGraphSize = graphSize;
    }

    /**
     *  Reinitialisation of the TGraph.
     */
    public void clear() {
        graph.clear();
        graphBitSet.clear();
    }

    /**
     *  Returns the graph object of this CDKRGraph.
     * @return      The graph object, a list
     */
    public List<CDKRNode> getGraph() {
        return this.graph;
    }

    /**
     *  Adds a new node to the CDKRGraph.
     * @param  newNode  The node to add to the graph
     */
    public void addNode(CDKRNode newNode) {
        graph.add(newNode);
        graphBitSet.set(graph.size() - 1);
    }

    /**
     *  Parsing of the CDKRGraph. This is the main method
     *  to perform a query. Given the constrains sourceBitSet and targetBitSet
     *  defining mandatory elements in G1 and G2 and given
     *  the search options, this method builds an initial set
     *  of starting nodes (targetBitSet) and parses recursively the
     *  CDKRGraph to find a list of solution according to
     *  these parameters.
     *
     * @param  sourceBitSet  constrain on the graph G1
     * @param  targetBitSet  constrain on the graph G2
     * @param  findAllStructure true if we want all results to be generated
     * @param  findAllMap true is we want all possible 'mappings'
     * @param timeManager
     * @throws CDKException
     */
    public void parse(BitSet sourceBitSet, BitSet targetBitSet, boolean findAllStructure, boolean findAllMap, TimeManager timeManager) throws CDKException {
        // initialize the list of solution

        this.timeManager = timeManager;
        // initialize the list of solution
        solutionList.clear();

        // builds the set of starting nodes
        // according to the constrains
        BitSet bitSet = buildB(sourceBitSet, targetBitSet);

        // setup options
        setAllStructure(findAllStructure);
        setAllMap(findAllMap);

        // parse recursively the CDKRGraph
        parseRec(new BitSet(bitSet.size()), bitSet, new BitSet(bitSet.size()));
    }

    /**
     *  Parsing of the CDKRGraph. This is the recursive method
     *  to perform a query. The method will recursively
     *  parse the CDKRGraph thru connected nodes and visiting the
     *  CDKRGraph using allowed adjacency relationship.
     *
     * @param  traversed  node already parsed
     * @param  extension  possible extension node (allowed neighbors)
     * @param  forbiden   node forbidden (set of node incompatible with the current solution)
     */
    private void parseRec(BitSet traversed, BitSet extension, BitSet forbidden) throws CDKException {
        BitSet newTraversed = null;
        BitSet newExtension = null;
        BitSet newForbidden = null;
        BitSet potentialNode = null;

        checkTimeOut();

        // if there is no more extension possible we
        // have reached a potential new solution
        if (extension.isEmpty()) {
            solution(traversed);
        } // carry on with each possible extension
        else {
            // calculates the set of nodes that may still
            // be reached at this stage (not forbidden)
            potentialNode = ((BitSet) graphBitSet.clone());
            potentialNode.andNot(forbidden);
            potentialNode.or(traversed);

            // checks if we must continue the search
            // according to the potential node set
            if (mustContinue(potentialNode)) {
                // carry on research and update iteration count
                nbIteration++;

                // for each node in the set of possible extension (neighbors of
                // the current partial solution, include the node to the solution
                // and parse recursively the CDKRGraph with the new context.
                for (int x = extension.nextSetBit(0); x >= 0 && !stop; x = extension.nextSetBit(x + 1)) {
                    // evaluates the new set of forbidden nodes
                    // by including the nodes not compatible with the
                    // newly accepted node.
                    newForbidden = (BitSet) forbidden.clone();
                    newForbidden.or((graph.get(x)).forbidden);

                    // if maxIterator is the first time we are here then
                    // traversed is empty and we initialize the set of
                    // possible extensions to the extension of the first
                    // accepted node in the solution.
                    if (traversed.isEmpty()) {
                        newExtension = (BitSet) ((graph.get(x)).extension.clone());
                    } // else we simply update the set of solution by
                    // including the neighbors of the newly accepted node
                    else {
                        newExtension = (BitSet) extension.clone();
                        newExtension.or((graph.get(x)).extension);
                    }

                    // extension my not contain forbidden nodes
                    newExtension.andNot(newForbidden);

                    // create the new set of traversed node
                    // (update current partial solution)
                    // and add x to the set of forbidden node
                    // (a node may only appear once in a solution)
                    newTraversed = (BitSet) traversed.clone();
                    newTraversed.set(x);
                    forbidden.set(x);

                    // parse recursively the CDKRGraph
                    parseRec(newTraversed, newExtension, newForbidden);
                }
            }
        }
    }

    /**
     * Checks if a potential solution is a real one
     * (not included in a previous solution)
     *  and add this solution to the solution list
     * in case of success.
     *
     * @param  traversed  new potential solution
     */
    private void solution(BitSet traversed) throws CDKException {
        boolean included = false;
        BitSet projG1 = projectG1(traversed);
        BitSet projG2 = projectG2(traversed);

        // the solution must follow the search constrains
        // (must contain the mandatory elements in G1 an G2)
        if (isContainedIn(sourceBitSet, projG1) && isContainedIn(targetBitSet, projG2)) {
            findNewSolutions(projG1, projG2, traversed, included);
            addNewSolution(included, traversed);
            isSubstructureFound();
        }
    }

    /**
     * Determine if there are potential solution remaining.
     * @param       potentialNode  set of remaining potential nodes
     * @return      true if maxIterator is worse to continue the search
     */
    private boolean mustContinue(BitSet potentialNode) throws CDKException {
        boolean cancel = false;
        BitSet projG1 = projectG1(potentialNode);
        BitSet projG2 = projectG2(potentialNode);

        // if we reached the maximum number of
        // search iterations than do not continue
        if (maxIteration != -1 && nbIteration >= maxIteration) {
            return false;
        }

        // if constrains may no more be fulfilled then stop.
        if (!isContainedIn(sourceBitSet, projG1) || !isContainedIn(targetBitSet, projG2)) {
            return false;
        }

        return checkPotentialSolution(projG1, projG2, cancel);
    }

    /**
     *  Builds the initial extension set. This is the
     *  set of node that may be used as seed for the
     *  CDKRGraph parsing. This set depends on the constrains
     *  defined by the user.
     * @param  sourceBitSet  constraint in the graph G1
     * @param  targetBitSet  constraint in the graph G2
     * @return
     */
    private BitSet buildB(BitSet sourceBitSet, BitSet targetBitSet) throws CDKException {
        this.sourceBitSet = sourceBitSet;
        this.targetBitSet = targetBitSet;

        BitSet bistSet = new BitSet();

        // only nodes that fulfill the initial constrains
        // are allowed in the initial extension set : targetBitSet
        for (Iterator<CDKRNode> i = graph.iterator(); i.hasNext();) {
            CDKRNode rNode = i.next();
            checkTimeOut();
            if ((sourceBitSet.get(rNode.getRMap().id1) || sourceBitSet.isEmpty()) && (targetBitSet.get(rNode.getRMap().id2) || targetBitSet.isEmpty())) {
                bistSet.set(graph.indexOf(rNode));
            }
        }
        return bistSet;
    }

    /**
     *  Returns the list of solutions.
     *
     * @return    The solution list
     */
    public List<BitSet> getSolutions() {
        return solutionList;
    }

    /**
     *  Converts a CDKRGraph bitset (set of CDKRNode)
     * to a list of CDKRMap that represents the
     * mapping between to substructures in G1 and G2
     * (the projection of the CDKRGraph bitset on G1
     * and G2).
     *
     * @param  set  the BitSet
     * @return      the CDKRMap list
     */
    public List<CDKRMap> bitSetToRMap(BitSet set) {
        List<CDKRMap> rMapList = new ArrayList<CDKRMap>();

        for (int x = set.nextSetBit(0); x >= 0; x = set.nextSetBit(x + 1)) {
            CDKRNode xNode = graph.get(x);
            rMapList.add(xNode.getRMap());
        }
        return rMapList;
    }

    /**
     *  Sets the 'AllStructres' option. If true
     * all possible solutions will be generated. If false
     * the search will stop as soon as a solution is found.
     * (e.g. when we just want to know if a G2 is
     *  a substructure of G1 or not).
     *
     * @param  findAllStructure
     */
    public void setAllStructure(boolean findAllStructure) {
        this.findAllStructure = findAllStructure;
    }

    /**
     *  Sets the 'finAllMap' option. If true
     * all possible 'mappings' will be generated. If false
     * the search will keep only one 'mapping' per structure
     * association.
     *
     * @param  findAllMap
     */
    public void setAllMap(boolean findAllMap) {
        this.findAllMap = findAllMap;
    }

    /**
     * Sets the maxIteration for the CDKRGraph parsing. If set to -1,
     * then no iteration maximum is taken into account.
     *
     * @param  maxIterator  The new maxIteration value
     */
    public void setMaxIteration(int maxIterator) {
        this.maxIteration = maxIterator;
    }

    /**
     *  Returns a string representation of the CDKRGraph.
     * @return the string representation of the CDKRGraph
     */
    @Override
    public String toString() {
        String message = "";
        int jIndex = 0;

        for (Iterator<CDKRNode> i = graph.iterator(); i.hasNext();) {
            CDKRNode rNode = i.next();
            message += "-------------\n" + "CDKRNode " + jIndex + "\n" + rNode.toString() + "\n";
            jIndex++;
        }
        return message;
    }

    /////////////////////////////////
    // BitSet tools
    /**
     *  Projects a CDKRGraph bitset on the source graph G1.
     * @param  set  CDKRGraph BitSet to project
     * @return      The associate BitSet in G1
     */
    public BitSet projectG1(BitSet set) {
        BitSet projection = new BitSet(firstGraphSize);
        CDKRNode xNode = null;

        for (int x = set.nextSetBit(0); x >= 0; x = set.nextSetBit(x + 1)) {
            xNode = graph.get(x);
            projection.set(xNode.getRMap().id1);
        }
        return projection;
    }

    /**
     *  Projects a CDKRGraph bitset on the source graph G2.
     * @param  set  CDKRGraph BitSet to project
     * @return      The associate BitSet in G2
     */
    public BitSet projectG2(BitSet set) {
        BitSet projection = new BitSet(secondGraphSize);
        CDKRNode xNode = null;

        for (int x = set.nextSetBit(0); x >= 0; x = set.nextSetBit(x + 1)) {
            xNode = graph.get(x);
            projection.set(xNode.getRMap().id2);
        }
        return projection;
    }

    /**
     *  Test if set sourceBitSet is contained in  set targetBitSet.
     * @param  sourceBitSet  a bitSet
     * @param  targetBitSet  a bitSet
     * @return    true if  sourceBitSet is contained in  targetBitSet
     */
    private boolean isContainedIn(BitSet sourceBitSet, BitSet targetBitSet) {
        boolean result = false;

        if (sourceBitSet.isEmpty()) {
            return true;
        }

        BitSet setA = (BitSet) sourceBitSet.clone();
        setA.and(targetBitSet);

        if (setA.equals(sourceBitSet)) {
            result = true;
        }

        return result;
    }

    private void isSubstructureFound() {
        if (!findAllStructure) {
            // if we need only one solution
            // stop the search process
            // (e.g. substructure search)
            stop = true;
        }
    }

    private void addNewSolution(boolean included, BitSet traversed) {
        if (included == false) {
            // if maxIterator is really a new solution add maxIterator to the
            // list of current solution
            solutionList.add(traversed);
        }
    }

    private void findNewSolutions(BitSet projG1, BitSet projG2, BitSet traversed, boolean included) throws CDKException {
        // the solution should not be included in a previous solution
        // at the CDKRGraph level. So we check against all previous solution
        // On the other hand if a previous solution is included in the
        // new one, the previous solution is removed.
        for (Iterator<BitSet> i = solutionList.listIterator(); i.hasNext() && !included;) {
            BitSet sol = i.next();
            checkTimeOut();
            if (!sol.equals(traversed)) {
                // if we asked to save all 'mappings' then keep this mapping
                if (findAllMap && (projG1.equals(projectG1(sol)) || projG2.equals(projectG2(sol)))) {
                    // do nothing
                } // if the new solution is included mark maxIterator as included
                else if (isContainedIn(projG1, projectG1(sol)) || isContainedIn(projG2, projectG2(sol))) {
                    included = true;
                } // if the previous solution is contained in the new one, remove the previous solution
                else if (isContainedIn(projectG1(sol), projG1) || isContainedIn(projectG2(sol), projG2)) {
                    i.remove();
                }
            } else {
                // solution already exists
                included = true;
            }
        }
    }

    private boolean checkPotentialSolution(BitSet projG1, BitSet projG2, boolean cancel) throws CDKException {

        boolean result = true;
        // check if the potential solution is not included in an already
        // existing solution
        for (Iterator<BitSet> i = solutionList.iterator(); i.hasNext() && !cancel;) {
            BitSet sol = i.next();
            checkTimeOut();
            // if we want every 'mappings' do not stop
            if (findAllMap && (projG1.equals(projectG1(sol)) || projG2.equals(projectG2(sol)))) {
                // do nothing
            } // if maxIterator is not possible to do better than an already existing solution than stop.
            else if (isContainedIn(projG1, projectG1(sol)) || isContainedIn(projG2, projectG2(sol))) {
                result = false;
                cancel = true;
            }
        }

        return result;
    }
}
