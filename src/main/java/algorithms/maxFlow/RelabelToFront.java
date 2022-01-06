package algorithms.maxFlow;

import algorithms.DoublyLinkedList;
import algorithms.VertexNode;
import algorithms.matrix.MatrixUtil;
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
     the neighbors are the edges into u and out of u.  
    */
    protected final TIntObjectMap<TIntList> uNMap;
    
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
    
    /** excess flow for a node: flow in exceeds flow out.  vertex s should have -inf.*/
    protected final double[] eF;
    
    protected final int srcIdx;
    
    protected final int sinkIdx;

    public RelabelToFront(TIntObjectMap<TIntSet> adj, 
        TObjectDoubleMap<PairInt> cap, int srcIdx, int sinkIdx) {
        
        this.uNMap = createAdjLists(adj);
        
        this.ell = constructL();
        
        assert(uNMap.size() == ell.size());
        
        this.c = MatrixUtil.copy(cap);
        
        if (srcIdx < 0 || !uNMap.containsKey(srcIdx)) {
            throw new IllegalArgumentException("srcIdx must be an index in adj");
        }
        
        this.srcIdx = srcIdx;
        this.sinkIdx = sinkIdx;
        
        int nV = uNMap.size();
        this.currNU = new int[nV];
        Arrays.fill(this.currNU, -1);
        
        this.f = initFlow();
        
        this.h = new int[nV];
        
        this.eF = new double[nV];
        Arrays.fill(this.eF, Double.NEGATIVE_INFINITY);
    }
    
    public MaxFlowResults findMaxFlow() {
        
        initPreFlow();
                
        initUNeighborListPointers();
                
        VertexNode u = ell.peekFirst();
        
        while (u != null) {
            
            dischargeLoop(u.vertex);
            
            u = (VertexNode) u.next;
        }
    }
    
    /**
     * if (u,v) in E: c_f(u,v)=c(u,v)-f(u,v); else if (v,u) in E c_f(u,v)=f(v,u), else=0
     * @param u
     * @param v
     * @return 
     */
    double calculateResidualCapacity(int u, int v) {
        
    }

    /**
     * min(c_f(u, v)) where (u,v) are in path p
     * @param p
     * @return 
     */
    double calculateResidualCapacity(PairInt[] p) {
        
    }

    /**
     * sum flow into u - sum of flow out of u.
      capacity constraint: for all (u,v) in E, 0 .lte. f(u,v) .lte. c(u,v).
      resid network: E_f={(u,v) pairs of V in which c_f(u,v) > 0;  
          the pairs can be existing edges in E or their reversals.
      f_p(u,v): { if (u,v) is on path p, f_p(u,v) = c_f(p), else f_p(u,v):=0.  
      note that p is an augmenting path in G_f
     * @param u
     * @param v
     * @return 
     */
    double calculateExcessFlow(int u, int v) {
        
    }

    void dischargeLoop(int vertex) {
        
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
            vList.add(vNode);
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
        
        return adjLists;
    }

    void initPreFlow() {
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
        TIntList nhbrList;
        int v, v2;
        int i;
        //for each vertex v in G.V, set h and e to 0
        while (iter.hasNext()) {
            iter.advance();
            v = iter.key();
            h[v] = 0;
            eF[v] = 0;
            
            // for each edge (u, v) in G.E, set f to 0
            nhbrList = iter.value();
            if (nhbrList == null || nhbrList.isEmpty()) {
                continue;
            }
            for (i = 0; i < nhbrList.size(); ++i) {
                v2 = nhbrList.get(i);
                eF[v2] = 0;
            }
        }
        
        h[srcIdx] = uNMap.size();
        nhbrList = uNMap.get(srcIdx);
        if (nhbrList == null) {
            throw new IllegalStateException("uNMap must have an entry for the srcIdx");
        }
        PairInt p;
        double cSV;
        for (i = 0; i < nhbrList.size(); ++i) {
            v = nhbrList.get(i);
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

    public static class MaxFlowResults {
        
    }
}
