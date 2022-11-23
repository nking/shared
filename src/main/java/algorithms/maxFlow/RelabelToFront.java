package algorithms.maxFlow;

import algorithms.DoublyLinkedList;
import algorithms.VertexNode;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectDoubleIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;

/**
 * <pre>
 * reference:
 * chap 26, and specifically 26.5 of Cormen, Leiserson, Rivest, & Stein "Introduction to Algorithms".
 * </pre>
 * 
 * we can interpret a directed graph as a “flow network” and use it to answer 
 * questions about material flows.  a material coursing through a system from 
 * a source, where the material is produced, to a sink, where it is consumed
 * is produced and consumed at the same steady rate.
 * 
 * In the maximum-flow problem, we wish to compute the greatest rate at which 
 * we can ship material from the source to the sink without violating any 
 * capacity constraints.
 * 
 * Relabel-To-Front is a maximum flow algorithm that begins with an
 * initialization of creating
 * a unchanging height of |V| for the source vertex and a unchanging height of 
 * 0 for the sink vertex, and modifiable heights of 0 for the vertexes in between
 * the source and sink.  The algorithm starts by sending flow from the 
 * source to each of its neighboring vertexes - the amount of flow is equal
 * to the capacity for the edge.  The goal for each vertex (except source and sink)
 * is to have 0 excess
 * flow (where the excess flow is the sum of incoming flow to a vertex minus the
 * sum of outgoing flow from the same vertex).  
 * The method for a vertex to reduce its excess flow to 0 is
 * called discharge and it's responsibilities are to invoke push of flow to
 * eligible neighbors, else relabel its height and continue to push and relabel
 * until the excess flow is 0.
 * At completion, the source vertex excess flow is -|flow| as it has given
 * that flow to the network, and the sink excess flow = the flow.  All other
 * vertexes will have excess flow=0.
 * 
 * runtime complexity is O(|V|^3).
 * 
 * @author nichole
 */
public class RelabelToFront {
    
    /**
     * represents a generous machine precision in comparisons to 0.
     * this could be reduced to 1e-11 or be a number the user can set.
     */
    protected final double eps = 1e-7;
    
    /** list of vertexes of G (minus source and sink) in any order, minus srcIdx and sinkIdx 
    modifications such as move node to front are O(1) for a doubly-linked list
    */
    protected final DoublyLinkedList<VertexNode> ell;
    
    /** rewrite graph adjacency into a map of a key which is the index of vertex u and the values which 
     are the neighbors of vertex u in ordered lists to be accessed by indexes.
     the neighbors are the edges into u and out of u aggregated into a single 
     sequential list. 
    */
    protected final TIntObjectMap<TIntList> uNMap;
    
    /**
     * a reference to the adjacency map representing the original graph G.
     * It is the edges out of a vertex.
     */
    protected final TIntObjectMap<TIntSet> adj;
    
    /**
     * an adjacency map of the edges into a vertex
     */
    protected final TIntObjectMap<TIntSet> revAdj;
     
    /** capacity of edges.  only positive entries are present.  no entries have .lte. 0 */
    protected final TObjectDoubleMap<PairInt> c;
    
    /** flow of edges*/
    protected final TObjectDoubleMap<PairInt> f;
    
    /** height of nodes */
    protected final int[] h;
    
    /** excess flow for a node: flow in exceeds flow out.*/
    protected final double[] eF;
    
    protected final int srcIdx;
    
    protected final int sinkIdx;

    public RelabelToFront(TIntObjectMap<TIntSet> adj, 
        TObjectDoubleMap<PairInt> cap, int srcIdx, int sinkIdx) {
        
        this.srcIdx = srcIdx;
        this.sinkIdx = sinkIdx;
        
        this.adj = adj;
        
        this.uNMap = createAdjLists(adj);
        
        this.revAdj = createIncomingEdgesVertexMap(adj);
        
        // contains V - {src, sink}
        this.ell = constructL();
        
        assert(uNMap.size() == (ell.size() + 2));
        
        this.c = MatrixUtil.copy(cap);
        
        if (srcIdx < 0 || !uNMap.containsKey(srcIdx)) {
            throw new IllegalArgumentException("srcIdx must be a key in adj");
        }
        if (sinkIdx < 0 || !uNMap.containsKey(sinkIdx)) {
            throw new IllegalArgumentException("sinkIdx must be a value in adj");
        }
        
        int nV = uNMap.size();
        
        // initalized to 0
        this.h = new int[nV];
        this.eF = new double[nV];
        
        // initalized to 0
        this.f = initFlow();
        
        initPreFlow();
        
        //System.out.printf("after initPreFlow excess[srcIdx]=%.3e = -|f*|%n", eF[srcIdx]);                
    }
    
    /**
     * 
     * @param adj graph as an adjacency matrix.  the elements in the matrix are
     * the edge capacities.
     * @param srcIdx
     * @param sinkIdx 
     */
    public RelabelToFront(double[][] adj, int srcIdx, int sinkIdx) {
        
        this.srcIdx = srcIdx;
        this.sinkIdx = sinkIdx;
        
        this.adj = createOutgoingEdgesVertexMap(adj);
        
        this.uNMap = createAdjLists(this.adj);
        
        this.revAdj = createIncomingEdgesVertexMap(this.adj);
        
        // contains V - {src, sink}
        this.ell = constructL();
        
        assert(uNMap.size() == (ell.size() + 2));
        
        this.c = createCapacityMap(adj);
        
        if (srcIdx < 0 || !uNMap.containsKey(srcIdx)) {
            throw new IllegalArgumentException("srcIdx must be a key in adj");
        }
        if (sinkIdx < 0 || !uNMap.containsKey(sinkIdx)) {
            throw new IllegalArgumentException("sinkIdx must be a value in adj");
        }
        
        int nV = uNMap.size();
        
        // initalized to 0
        this.h = new int[nV];
        this.eF = new double[nV];
        
        // initalized to 0
        this.f = initFlow();
        
        initPreFlow();
        
        //System.out.printf("after initPreFlow excess[srcIdx]=%.3e = -|f*|%n", eF[srcIdx]);                
    }
    
    public MaxFlowResults findMaxFlow() {
        
        VertexNode u = ell.peekFirst();
        
        while (u != null) {
            
            dischargeLoop(u);
            
            // u may now be at the front of the list (modified in dischargeLoop), 
            //    and in that case, the next u is the 2nd in the list ell.
            
            u = (VertexNode) u.next;
        }
        
        MaxFlowResults r = new MaxFlowResults();
        r.srcIdx = this.srcIdx;
        r.sinkIdx = this.sinkIdx;
        r.flow = this.eF[sinkIdx];
        r.edgeFlows = MatrixUtil.copy(this.f);
        
        //print();
        
        assertZeroExcess();
        
        return r;
    }
    
    /**
     * if (u,v) in E: c_f(u,v)=c(u,v)-f(u,v); else if (v,u) in E c_f(u,v)=f(v,u), else=0
     * @param u
     * @param v
     * @return the residual capacity
     */
    protected double calculateResidualCapacity(int u, int v) {
        return calculateResidualCapacity(new PairInt(u, v));
    }
    
    /**
     * if (u,v) in E: c_f(u,v)=c(u,v)-f(u,v); else if (v,u) in E c_f(u,v)=f(v,u), else=0
     * @param p edge (u, v)
     * @return the residual capacity
     */
    protected double calculateResidualCapacity(PairInt p) {
        
        if (this.adj.containsKey(p.getX()) &&
            this.adj.get(p.getX()).contains(p.getY())) {
            // is an edge in G.E
            double cUV = this.c.get(p);
            double fUV = this.f.get(p);
            return cUV - fUV;
        } else if (this.revAdj.containsKey(p.getY()) &&
            this.revAdj.get(p.getY()).contains(p.getX())) {
            return this.f.get(new PairInt(p.getY(), p.getX()));
        } else {
            // is not an edge in G.E
            return 0;
        }
    }
    
    protected void printF() {
        printF(f);
    }
    
    protected static void printF(TObjectDoubleMap<PairInt> flow) {
         TObjectDoubleIterator<PairInt> iter = flow.iterator();
         PairInt p;
         double fP;
         while (iter.hasNext()) {
             iter.advance();
             p = iter.key();
             fP = iter.value();
             System.out.printf("f%s = %.3e%n", p.toString(), fP);
         }
    }
    
    protected void print() {
        printF();
        printC();
        printE();
        printH();
    }
    
    protected void printC() {
         TObjectDoubleIterator<PairInt> iter = c.iterator();
         PairInt p;
         double cP;
         while (iter.hasNext()) {
             iter.advance();
             p = iter.key();
             cP = iter.value();
             System.out.printf("c%s = %.3e%n", p.toString(), cP);
         }
    }
    
    protected void printE() {
         System.out.printf("e=%s%n", FormatArray.toString(eF, "%.3e"));
    }
    
    protected void printH() {
         System.out.printf("h=%s%n", Arrays.toString(h));
    }

    /**
     * sum(flow into u) - sum(flow out of u).
      capacity constraint: for all (u,v) in E, 0 .lte. f(u,v) .lte. c(u,v).
      residual network: E_f={(u,v) pairs of V in which c_f(u,v) > 0;  
          the pairs can be existing edges in E or their reversals.
      f_p(u,v): { if (u,v) is on path p, f_p(u,v) = c_f(p), else f_p(u,v):=0.  
      note that p is an augmenting path in G_f
     * @param u
     * @return 
     */
    protected double calculateExcessFlow(int u) {
        double sumToU = calculateSumOfFlow(u, this.revAdj.get(u));
        double sumFromU = calculateSumOfFlow(u, this.adj.get(u));
        return sumToU - sumFromU;
    }
    
    /**
     * sum(flow into u) - sum(flow out of u).
      capacity constraint: for all (u,v) in E, 0 .lte. f(u,v) .lte. c(u,v).
      residual network: E_f={(u,v) pairs of V in which c_f(u,v) > 0;  
          the pairs can be existing edges in E or their reversals.
      f_p(u,v): { if (u,v) is on path p, f_p(u,v) = c_f(p), else f_p(u,v):=0.  
      note that p is an augmenting path in G_f
     * @param u
     * @param v the neighbors of u
     * @return 
     */
    protected double calculateSumOfFlow(int u, TIntSet v) {
        if (u == srcIdx) {
            return 0;
        }
        PairInt p;
        double sum = 0;
        TIntIterator iter = v.iterator();
        while (iter.hasNext()) {
            p = new PairInt(u, iter.next());
            assert(this.f.containsKey(p));
            sum += this.f.get(p);
        }
        return sum;
    }

    /**
     * @param uNode
     */
    protected void dischargeLoop(VertexNode uNode) {
        /*
        old-height = u.h
        discharge(u);
        if (u.h > old-height)
            move u to the front of the list L
        */
        
        int u = uNode.vertex;
        
        double oH = h[u];
        discharge(u);
        if (h[u] > oH) {
            ell.unlink(uNode);
            ell.addFirst(uNode);
        } 
    }

    /**
     * construct a linked list of all vertexes except src and sink
     * @return 
     */
    private DoublyLinkedList<VertexNode> constructL() {
        
        if (uNMap == null) {
            throw new IllegalStateException("uNMap cannot be null");
        }
        
        DoublyLinkedList<VertexNode> vList = new DoublyLinkedList<VertexNode>();
        
        int u;
        VertexNode vNode;
        
        // uNMap has all vertexes as keys
        TIntObjectIterator<TIntList> iter = uNMap.iterator();
        while (iter.hasNext()) {
            iter.advance();
            u = iter.key();
            if (u == srcIdx || u == sinkIdx) {
                continue;
            }
            vNode = new VertexNode(u);
            vList.addFirst(vNode);
        }
        
        return vList;
    }
    
    private TIntObjectMap<TIntList> createAdjLists(TIntObjectMap<TIntSet> adj) {
        
        adj = MatrixUtil.copyToSymmetricMap(adj);
        
        TIntObjectMap<TIntList> adjLists = new TIntObjectHashMap<TIntList>();
        
        TIntIterator iter2;
        TIntSet uSet;
        int u;
        int v;
        TIntList nhbrList;
        
        TIntObjectIterator<TIntSet> iter = adj.iterator();
        while (iter.hasNext()) {
            iter.advance();
            u = iter.key();
            
            nhbrList = new TIntArrayList();
            adjLists.put(u, nhbrList);
            
            uSet = iter.value();
            if (uSet == null || uSet.isEmpty()) {
                continue;
            }
            iter2 = uSet.iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                nhbrList.add(v);
            }
        }
        
        //TODO: use a design pattern to inject this only for tests.
        //    e.g. aspects can follow the method with byte code enhancement
        //    that adds this after method completes.
        //this loop is only used in tests to reverse the order in each nhbrList
        TIntObjectIterator<TIntList> iter3 = adjLists.iterator();
        TIntObjectMap<TIntList> out = new TIntObjectHashMap<TIntList>();
        while (iter3.hasNext()) {
            iter3.advance();
            u = iter3.key();
            nhbrList = iter3.value();
            nhbrList.reverse();
            out.put(u, nhbrList);
        }
        adjLists = out;
        
        return adjLists;
    }

    private void initPreFlow() {
        if (uNMap == null) {
            throw new IllegalStateException("uNMap cannot be null");
        }
        /*for each vertex v in G.V   (already handled in constructor)
             v.h = 0;
             v.e = 0;
          for each edge (u, v) in G.E (already handled in constructor)
             (u,v).f = 0
          s.h = |G.V|
          for each vertex v in s.Adj
            (s,v).f = c(s,v);
            v.e = c(s,v)
            s.e = s.e - c(s,v)
        */
        
        h[srcIdx] = uNMap.size();
        
        TIntSet nhbrSet = adj.get(srcIdx);
        if (nhbrSet == null) {
            throw new IllegalStateException("adj must have an entry for the srcIdx");
        }
        PairInt p;
        double cSV;
        int v;
        TIntIterator iter3 = nhbrSet.iterator();
        while (iter3.hasNext()) {
            v = iter3.next();
            p = new PairInt(srcIdx, v);
            cSV = this.c.get(p);
            
            //(s,v).f = c(s,v);
            this.f.put(p, cSV);
            
            //v.e = c(s,v)
            this.eF[v] = cSV;
            //s.e = s.e - c(s,v)
            this.eF[srcIdx] -= cSV;
        }
    
    }
    
    /**
     * initialize flow to 0 for every edge in capacity map
     * @return 
     */
    private TObjectDoubleMap<PairInt> initFlow() {
        
        if (c == null) {
            throw new IllegalStateException("c cannot be null");
        }
        
        TObjectDoubleMap<PairInt> out = new TObjectDoubleHashMap<PairInt>();
        
        TObjectDoubleIterator<PairInt> iter = c.iterator();
        PairInt p;
        while (iter.hasNext()) {
            iter.advance();
            p = iter.key();
            // initalize to 0
            out.put(p.copy(), 0);
        }
        
        return out;
    }

    protected TIntObjectMap<TIntSet> createIncomingEdgesVertexMap(TIntObjectMap<TIntSet> adj) {
        return MatrixUtil.createReverseMap(adj);
    }
    
    /**
     * if u is overflowing and all neighbors in E_f are not downhill from u,
     * relabel increases the height of u to be 1 plus the minimum of the
     * heights of the neighbors of u that are in E_f.
     * @throws java.lang.IllegalArgumentException throws exception if a neighbor
     * of u is eligible for a push operation.
     * @param u
     */
    protected void relabel(int u) {
        
        if (u == srcIdx || u == sinkIdx) {
            throw new IllegalArgumentException("u cannot be srcIdx or sinkIdx");
        }
        
        TIntList vs = uNMap.get(u);
        
        /*System.out.printf("relabel(%d).  N=%s%n", u, Arrays.toString(vs.toArray()));
        print();
        System.out.printf("  for V_f: ");*/
        
        int v;
        double cF;
        boolean inE;
        int minH = Integer.MAX_VALUE;
        for (int i = 0; i < vs.size(); ++i) {
            v = vs.get(i);
            
            inE = this.adj.containsKey(u) && this.adj.get(u).contains(v);
            
            //System.out.printf("  inE=%b h[%d]=%d  h[%d]=%d%n", inE, u, h[u], v, h[v]);
            
            // Lemma 26.16
            if (this.h[v] < this.h[u]) {
                continue;
            }
                                    
            if (inE) {
                cF = calculateResidualCapacity(u, v);
                                
                if (Math.abs(cF) < 1e-7) {
                    // not an edge in G.E
                    continue;
                }
            }
            
            // edge is in E_f
            
            //System.out.printf("h[%d]=%d, ", v, this.h[v]);
            
            // calculate minimum of heights of u neighbors in E_f
            if (this.h[v] < minH) {
                minH = this.h[v];
            }
        }
        
        this.h[u] = minH + 1;
        
        //Lemma 26.20
        assert(h[u] <= (2*uNMap.size() - 1));
        
        //DEBUG
        /*System.out.printf("%n");        
        System.out.printf("  minH=%d,  h=%s%n", minH, u, Arrays.toString(h));        
        if (!(minH < Integer.MAX_VALUE)) {
            System.out.printf("ERROR in relabel(%d).  N=%s%n",
                u, Arrays.toString(vs.toArray()));
        }*/
        
        assert(minH < Integer.MAX_VALUE);     
    }
    
    /**
     * The basic operation PUSH(u,v) applies if u is an overflowing vertex, 
     * c_f(u,v) > 0, and h.u = h.v + 1.
     * @param u
     * @param v 
     * @throws java.lang.IllegalArgumentException throws exception cf_(u,v) is
     * not a positive number or h.u != (h.v + 1).
     */
    protected void push(int u, int v) {
        
        //System.out.printf("push(%d, %d)%n", u, v);
        
        boolean inE = this.adj.containsKey(u) && this.adj.get(u).contains(v);
        
        if (this.h[u] != (this.h[v] + 1)) {
            throw new IllegalStateException("cannot push because "+u+".h != ("+v+
                ".h + 1)"); 
        }
        
        double cf = calculateResidualCapacity(u, v);
        if (inE && cf <= 0) {
            throw new IllegalStateException("cannot push because c_f("+u+","+v+
                ")=" + cf + " which is not positive"); 
        }
        
        double delta;
        if (inE) {
            delta = Math.min(eF[u], cf);
            PairInt p = new PairInt(u, v);
            double fUV = f.get(p) + delta;
            f.put(p, fUV);
            
            //DEBUG
            /*System.out.printf("  ==> f(%d,%d)=%.3e, c(%d,%d)=%.3e, delta=%.3e%n", 
                u, v, fUV, u, v, c.get(p), delta);
            if (Math.abs(c.get(p) - fUV) < 1e-7) {
                System.out.printf("     saturated%n");
            }*/
            
        } else {
            double fVU = 0;
            PairInt p = new PairInt(v, u);
            if (f.containsKey(p)) {
                fVU = f.get(p);
            }
            // edit to algorithm, similar to residual capacity for back-edge
            delta = Math.min(eF[u], fVU);
            
            fVU -= delta;
            
            assert(fVU <= c.get(p));
            
            f.put(p, fVU);
                        
            //DEBUG
            //System.out.printf("  ==> f(%d,%d)=%.3e, c(%d,%d)=%.3e delta=%.3e for v,u%n", 
            //    v, u, fVU, v, u, c.get(p), delta);
        }
        this.eF[u] -= delta;
        this.eF[v] += delta;
        
        /*
        //DEBUG
        System.out.printf("  eF[%d]=%.3e, eF[%d]=%.3e%n", u, eF[u], v, eF[v]);
        
        if (Math.abs(this.eF[u]) <1e-7) {
            System.out.printf("      removes %d from E_f%n", u);
        }
        if (Math.abs(this.eF[v]) <1e-7) {
            System.out.printf("      removes %d from E_f%n", v);
        }*/
    }
    
    protected void discharge(int u) {

        /*//DEBUG        
        print();                
        System.out.printf("discharge(%d)  excess=%.3e%n", u, eF[u]);
        */
        
        TIntList uNList = uNMap.get(u);
        if (uNList == null) {
            throw new IllegalStateException("vertex " + u + " has no neighbors");
        }
        int v;
        int nIter = 1;
        double cF;
        boolean inE;
                
        TIntObjectMap<TIntList> neighborHeightMap = createNeighborHeightMap(uNList);
        
        int currNIdx = 0;
        TIntList currNList;
        
        DoublyLinkedList<VertexNode> backEdgeIdxs = new DoublyLinkedList<VertexNode>();
        
        while (eF[u] > 0) {
                        
            currNList = neighborHeightMap.get(h[u] - 1);
            
            /*
            //DEBUG
            System.out.printf("  nIter=%d eF[%d]=%.3e%n", nIter, u, eF[u]);
            if (currNList != null) {
                System.out.printf("  [%d].N=%s u.h=%d [h-1].currNList=%s%n", u, Arrays.toString(uNList.toArray()),
                    h[u], Arrays.toString(currNList.toArray()));
            }*/
                        
            if (currNList == null || currNIdx >= currNList.size()) {
                
                if (!backEdgeIdxs.isEmpty()) {
                    
                    // for a height, push to the back-edges.
                    currNIdx = backEdgeIdxs.removeFirst().vertex;
                    
                    v = currNList.get(currNIdx);
                    
                    //System.out.printf("    <-*cf(%d,%d)=%.3e, h[%d]=%d, h[%d]=%d  inE=false%n", 
                    //    u, v, calculateResidualCapacity(u, v), u, h[u], v, h[v]);
                    
                    push(u, v);
                    
                    // set to a value that results in the v=nil condition above
                    currNIdx = currNList.size();
                    
                    nIter++;
                    
                    continue;
                }
                
                relabel(u);
                
                currNIdx = 0;
                
                backEdgeIdxs.clear();
                
                //System.out.printf("    after relabel returning to discharge(%d).  reset u.currN%n", u);
                
                nIter++;
                
                continue;
            }
            
            v = currNList.get(currNIdx);
            
            inE = this.adj.containsKey(u) && this.adj.get(u).contains(v);
                        
            cF = calculateResidualCapacity(u, v);
            
            if (cF > 0) {
                
                //System.out.printf("    *cf(%d,%d)=%.3e, h[%d]=%d, h[%d]=%d  inE=%b%n", u, v, cF, u, h[u], v, h[v], inE);
                
                push(u, v);
                
            } else if (!inE && cF <= 0) {
                
                backEdgeIdxs.add(new VertexNode(currNIdx));
                
                //System.out.printf("did not push to reverse edge for v=%d cF=%.3e"
                //    + " currNList=%s%n",
                //    v, cF, Arrays.toString(currNList.toArray()));
                
                currNIdx++;
                
            } else {
                //System.out.printf("    cf(%d,%d)=%.3e, h[%d]=%d, h[%d]=%d inE=%b. incr u.currN%n", u, v, cF, u, h[u], v, h[v], inE);
                currNIdx++;
            }
            
            nIter++;
        }
    }

    /**
     * assert each  vertex in V - {s,t} has an excess of 0
     */
    private void assertZeroExcess() {
        
        VertexNode u = ell.peekFirst();
        
        while (u != null) {
            
            assert(Math.abs(this.eF[u.vertex]) < 1e-7);
            
            u = (VertexNode) u.next;
        }
    }

    private TIntObjectMap<TIntList> createNeighborHeightMap(TIntList list) {
        
        TIntObjectMap<TIntList> map = new TIntObjectHashMap<TIntList>();
        
        int v;
        int h;
        TIntList mList;
        for (int i = 0; i < list.size(); ++i) {
            v = list.get(i);
            h = this.h[v];
            mList = map.get(h);
            if (mList == null) {
                mList = new TIntArrayList();
                map.put(h, mList);
            }
            mList.add(v);
        }
        return map;
    }

    private TIntObjectMap<TIntSet> createOutgoingEdgesVertexMap(double[][] adj) {
        
        TIntObjectMap<TIntSet> map = new TIntObjectHashMap<TIntSet>();
        int u;
        int v;
        TIntSet set;
        for (u = 0; u < adj.length; ++u) {
            
            set = map.get(u);
            if (set == null) {
                set = new TIntHashSet();
                map.put(u, set);
            }
            
            for (v = 0; v < adj[u].length; ++v) {
                if (Math.abs(adj[u][v]) > eps) {
                    set.add(v);
                }
            }
            
            if (set.isEmpty()) {
                map.remove(u);
            }
        }
        
        return map;
    }

    private TObjectDoubleMap<PairInt> createCapacityMap(double[][] adjM) {
        
        TObjectDoubleMap<PairInt> map = new TObjectDoubleHashMap<PairInt>();
        
        int u;
        int v;
        PairInt p;
        for (u = 0; u < adjM.length; ++u) {
            for (v = 0; v < adjM[u].length; ++v) {
                if (Math.abs(adjM[u][v]) > eps) {
                    map.put(new PairInt(u, v), adjM[u][v]);
                }
            }
        }
        
        return map;
    }

    public static class MaxFlowResults {
        public TObjectDoubleMap<PairInt> edgeFlows;
        public int srcIdx;
        public int sinkIdx;
        public double flow;
        public void print() {
            if (edgeFlows != null) {
                System.out.printf("max flow=%.3e from src=%d to sink=%d%n", 
                    flow, srcIdx, sinkIdx);
                printF(edgeFlows);
            }
        }
    }
}
