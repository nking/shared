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
 * Relabel-To-Front is a maximum flow algorithm.
 * 
 * runtime complexity is O(|V|^3).
 * 
 * @author nichole
 */
public class RelabelToFront {
    
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
    
    /** the indexes are the U vertex numbers.  the values are the current neighbor 
    index (which is vertex v index).
    if there are no neighbors for vertex U, the value is -1.
    */
    protected final int[] currNU;
    
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
        
        this.revAdj = createMapIntoVertices(adj);
        
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
        this.currNU = new int[nV];
        Arrays.fill(this.currNU, -1);
        
        this.f = initFlow();
        
        this.h = new int[nV];
        
        this.eF = new double[nV];
        Arrays.fill(this.eF, Double.NEGATIVE_INFINITY);
        
        initPreFlow();
                
        initUNeighborListPointers();
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
        r.flow = this.f.get(new PairInt(srcIdx, sinkIdx));
        r.edgeFlows = MatrixUtil.copy(this.f);
        
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
         TObjectDoubleIterator<PairInt> iter = f.iterator();
         PairInt p;
         double fP;
         while (iter.hasNext()) {
             iter.advance();
             p = iter.key();
             fP = iter.value();
             System.out.printf("f%s = %.3e\n", p.toString(), fP);
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
             System.out.printf("c%s = %.3e\n", p.toString(), cP);
         }
    }
    
    protected void printE() {
         System.out.printf("e=%s\n", FormatArray.toString(eF, "%.3e"));
    }
    protected void printH() {
         System.out.printf("h=%s\n", Arrays.toString(h));
    }

    /*
     * min(c_f(u, v)) where (u,v) are in path p
     * @param p
     * @return 
    double calculateMinResidualCapacity(PairInt[] p) {
        double min = Double.POSITIVE_INFINITY;
        double r;
        for (PairInt pi : p) {
            r = calculateResidualCapacity(pi);
            if (r < min) {
                min = r;
            }
        }
        return min;
    }
    */

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
        /*for each vertex v in G.V
             v.h = 0;
             v.e = 0;
          for each edge (u, v) in G.E
             (u,v).f = 0
          s.h = |G.V|
          for each vertex v in s.Adj
            (s,v).f = c(s,v);
            v.e = c(s,v)
            s.e = s.e - c(s,v)
        */
        
        TIntObjectIterator<TIntList> iter = uNMap.iterator();
        int v, v2;
        int i;
        //for each vertex v in G.V, set h and e to 0
        while (iter.hasNext()) {
            iter.advance();
            v = iter.key();
            h[v] = 0;
            eF[v] = 0;
        }
            
        // for each edge (u, v) in G.E, set f to 0
        //TODO: consider whether to set (v, u) to 0 also
        TIntObjectIterator<TIntSet> iter2 = adj.iterator();
        TIntSet nhbrSet;
        TIntIterator iter3;
        while (iter2.hasNext()) {
            iter2.advance();
            v = iter2.key();
            
            nhbrSet = iter2.value();
            if (nhbrSet == null || nhbrSet.isEmpty()) {
                continue;
            }
            iter3 = nhbrSet.iterator();
            while (iter3.hasNext()) {
                v2 = iter3.next();
                eF[v2] = 0;
            }
        }
        
        h[srcIdx] = uNMap.size();
        
        nhbrSet = adj.get(srcIdx);
        if (nhbrSet == null) {
            throw new IllegalStateException("adj must have an entry for the srcIdx");
        }
        PairInt p;
        double cSV;
        iter3 = nhbrSet.iterator();
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
    
    private void initUNeighborListPointers() {
        VertexNode u = ell.peekFirst();
        int uIdx;
        while (u != null) {
            uIdx = u.vertex;
            if (uNMap.containsKey(uIdx) && !uNMap.get(uIdx).isEmpty()) {
                this.currNU[uIdx] = 0;
            } else {
                this.currNU[uIdx] = -1;
            }
            u = (VertexNode) u.next;
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
            out.put(p.copy(), 0);
        }
        
        return out;
    }

    protected TIntObjectMap<TIntSet> createMapIntoVertices(TIntObjectMap<TIntSet> adj) {
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
        
        System.out.printf("relabel(%d).  N=%s\n", u, Arrays.toString(vs.toArray()));
        
        //DEBUG
            printF();
            printC();
            printE();
            printH();
        
        System.out.printf("u=%d\n", u);
        
        System.out.printf("  for V_f=");
        int v;
        double cF;
        int minH = Integer.MAX_VALUE;
        for (int i = 0; i < vs.size(); ++i) {
            v = vs.get(i);
            cF = calculateResidualCapacity(u, v);
            if (cF == 0) {
                // not an edge in G.E
                continue;
            }
            // edge is in E_f
            if (this.h[u] > this.h[v]) { // or more specifically if h.u == (h.v + 1) ?
                // v is downhill from u so can receive a push, making relabel invalid
                throw new IllegalStateException("cannot relabel because there is a "
                    + "neighboring vertex eligible for a push");
            }
            System.out.printf("h[%d]=%d, ", v, this.h[v]);
            // calculate minimum of heights of u neighbors in E_f
            if (this.h[v] < minH) {
                minH = this.h[v];
            }
        }
        System.out.printf("\n");        
        System.out.printf("  minH=%d\n", minH);
        
        //DEBUG
        if (!(minH < Integer.MAX_VALUE)) {
            System.out.printf("ERROR in relabel(%d).  N=%s\n",
                u, Arrays.toString(vs.toArray()));
        }
        
        assert(minH < Integer.MAX_VALUE);
        
        this.h[u] = minH + 1;        
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
        
        System.out.printf("push(%d, %d)\n", u, v);
        
        assert(uNMap.containsKey(u) && uNMap.get(u).contains(v));
        
        double cf = calculateResidualCapacity(u, v);
        if (cf <= 0) {
            throw new IllegalStateException("cannot push because c_f("+u+","+v+
                ")=" + cf + " which is not positive"); 
        }
        if (this.h[u] != (this.h[v] + 1)) {
            throw new IllegalStateException("cannot push because "+u+".h != ("+v+
                ".h + 1)"); 
        }
        double delta = Math.min(eF[u], cf);
        if (adj.containsKey(u) && adj.get(u).contains(v)) {
            PairInt p = new PairInt(u, v);
            double fUV = f.get(p) + delta;
            System.out.printf("  ==> %.3e\n", fUV);
            f.put(p, fUV);
        } else {
            double fVU = 0;
            PairInt p = new PairInt(v, u);
            if (f.containsKey(p)) {
                fVU = f.get(p);
            }
            fVU -= delta;
            f.put(p, fVU);
            System.out.printf("  ==> %.3e\n", fVU);
        }
        this.eF[u] -= delta;
        this.eF[v] += delta;
    }
    
    protected void discharge(int u) {
        
        System.out.printf("discharge(%d)\n", u);
        
        TIntList uNList = uNMap.get(u);
        if (uNList == null) {
            throw new IllegalStateException("vertex " + u + " has no neighbors");
        }
        int uNListIdx, v;
        while (eF[u] > 0) {
            uNListIdx = this.currNU[u];
            if (uNListIdx >= uNList.size()) {
                relabel(u);
                this.currNU[u] = 0;
                continue;
            }
            v = uNList.get(uNListIdx);
            double cf = calculateResidualCapacity(u, v);
            if (cf > 0 && (h[u] == (h[v] + 1))) {
                push(u, v);
            } else {
                this.currNU[u]++;
            }
        }
    }

    public static class MaxFlowResults {
        public TObjectDoubleMap<PairInt> edgeFlows;
        public int srcIdx;
        public int sinkIdx;
        public double flow;
    }
}
