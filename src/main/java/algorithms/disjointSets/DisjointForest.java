package algorithms.disjointSets;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * An "all-in-one" class to handle the disjoint set methods and the nodes
 * in a disjoint forest of trees.
 *
 * The disjointSet representation has a path compression operation internal to
 * findSet method which makes subsequent membership disjointSet queries
 * faster.
 * 
 * implementation follows that in Cormen et al. Introduction to Algorithms.
 * 
 * @author nichole
 */
public class DisjointForest<T> {

    // key = parent nodes
    // value = disjoint sets
    private Map<DisjointSet2Node<T>, RootedTreeDisjointSet<T>> trees =
        new HashMap<>();
    
    /**
     * make a set out of the given node and add it to the internal forest.
     * runtime complexity is O(1).
     * @param x
     * @return
     */
    public RootedTreeDisjointSet<T> makeSet(DisjointSet2Node<T> x) {
        x.setParent(x);
        x.setRank(0);
        
        RootedTreeDisjointSet<T> set = new RootedTreeDisjointSet<T>(x);
        addToForest(set);
        return set;
    }
    
     /**
     * <pre>
     * find the set representative for the given node.  As a side effect,
     * also updates x and all of it's ancestors with the found result so that
     * they directly point to the top most parent and subsequent lookups are
     * faster (this is path compression).
     * runtime complexity:
     *     the method uses iteration.  
     *     if we represent x_height as the number of nodes between x's tree 
     *     root and x, we have 2 times x_height iterations of statements,
     * so O(x_height).   With path compression here, the amoritzed
     * worst-case running time is Θ(α(n)) where α is the inverse Ackermann 
     * function.  
     * 
     * The inverse Ackermann function is α(n) = min{k : A_k(1) ≥ n}.
     * For most purposes, α(n) = O(1) so then the amoritized running time is O(1).
     * </pre>
     * 
     * @param x
     * @return
     */
    public DisjointSet2Node<T> findSet(DisjointSet2Node<T> x) {
                
        // iterative
        if (!x.equals(x.getParent())) {
            
            List<DisjointSet2Node<T>> update = new ArrayList<DisjointSet2Node<T>>();
            
            DisjointSet2Node<T> parent = x;
            while (!parent.equals(parent.getParent())) {
                update.add(parent);
                parent = parent.getParent();
            }
            
            // update the nodes with parent
            for (DisjointSet2Node<T> node : update) {
                node.setParent(parent);
            }
 
        }
        
        return x.getParent();
    }
    
    /**
      assigns to x and y the same parent from both of their parents, choosing
      the one with largest rank.
          
       Runtime complexity is O(1).
       
     * @param x
     * @param y
     * @return the root found to be the one with equal number of nodes or more nodes
     */
    private DisjointSet2Node<T> link(DisjointSet2Node<T> x, DisjointSet2Node<T> y) {
        
        if (x.equals(y)) {
            return x;
        }
        
        DisjointSet2Node<T> parent;
        if (x.getRank() >= y.getRank()) {
            parent = x;
            y.setParent(parent);
            if (x.getRank() == y.getRank()) {
                parent.setRank(parent.getRank() + 1);
            }
        } else {
            parent = y;
            x.setParent(parent);
        }
        return parent;
    }
    
    /**
      assigns to x and y the same parent from both of their parents, choosing
      the y as the parent.
          
       Runtime complexity is O(1).
       
     * @param x
     * @param y
     * @return the root found to be the one with equal number of nodes or more nodes
     */
    private DisjointSet2Node<T> linkChooseY(DisjointSet2Node<T> x, DisjointSet2Node<T> y) {
        
        if (x.equals(y)) {
            return y;
        }
        
        DisjointSet2Node<T> parent = y;
        x.setParent(parent);
        if (x.getRank() == y.getRank()) {
            parent.setRank(parent.getRank() + 1);
        }
        
        return parent;
    }
    
    private void addToForest(RootedTreeDisjointSet<T> t) {
        getTrees().put(t.parent, t);
    }
    private RootedTreeDisjointSet<T> removeFromForest(DisjointSet2Node<T> key) {
        return getTrees().remove(key);
    }
    
    /**
     * append the shorter list onto the end of the longer list.
     * The amortized runtime complexity is O(1) due to the path compression
     * in findSet.
     * The method is also known as "union-find" because it uses findSet 
     * as the first steps before link.
     * 
     * @param x
     * @param y
     * @return
     */
    public DisjointSet2Node<T> union(DisjointSet2Node<T> x, DisjointSet2Node<T> y) {
        
        x = findSet(x);
        
        y = findSet(y);
        
        if (x.equals(y)) {
            return x;
        }
        
        RootedTreeDisjointSet<T> tX = removeFromForest(x.parent);
        RootedTreeDisjointSet<T> tY = removeFromForest(y.parent);
        
        DisjointSet2Node<T> parent = link(x, y);
        
        tX.nodes.addAll(tY.nodes);
        addToForest(tX);

        return parent;
    }
    
    /**
     * Compress the paths in x and y and set the parent to be y.
\     * The amortized runtime complexity is O(1) due to the path compression
     * in findSet.
     * The method is also known as "union-find" because it uses findSet 
     * as the first steps before link.
     * 
     * @param x
     * @param y
     * @return
     */
    public DisjointSet2Node<T> unionChooseY(DisjointSet2Node<T> x, DisjointSet2Node<T> y) {
        
        x = findSet(x);
        
        y = findSet(y);
        
        if (x.equals(y)) {
            return y;
        }
        
        RootedTreeDisjointSet<T> tX = removeFromForest(x.parent);
        RootedTreeDisjointSet<T> tY = removeFromForest(y.parent);
        
        DisjointSet2Node<T> parent = linkChooseY(x, y);
        
        tY.nodes.addAll(tX.nodes);
        
        addToForest(tY);
        
        return parent;
    }
  
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("trees.size=");
        if (getTrees() != null) {
            sb.append(getTrees().size());
        } else {
            sb.append("0");
        }
        sb.append(", [trees=");
        if (getTrees() != null) {
            Iterator<Entry<DisjointSet2Node<T>, RootedTreeDisjointSet<T>>> iter 
                = getTrees().entrySet().iterator();
            Entry<DisjointSet2Node<T>, RootedTreeDisjointSet<T>> entry;
            DisjointSet2Node<T> p;
            RootedTreeDisjointSet<T> t;
            while (iter.hasNext()) {
                entry = iter.next();
                p = entry.getKey();
                t = entry.getValue();
                sb.append("\n tree p=");
                if (p.member != null) {
                    sb.append(p.member);
                } else {
                    sb.append(p.hashCode());
                }
                sb.append(", t=").append(t.toString());
                sb.append(", ");
            }
        }
        sb.append(", ");
        sb.append("] ");

        return sb.toString();
    }

    /**
     * each disjoint set (== instance of this class), has a representative which
     * is a member of the nodes of this set
     *
     * @author nichole
     */
    public static class RootedTreeDisjointSet<T> {

        /**
         * the representative of the tree. it's the root node.
         */
        protected DisjointSet2Node<T> repr = null;

        /**
         * pointer to the tree parent in the forest
         */
        protected DisjointSet2Node<T> parent = null;

        private Set<DisjointSet2Node<T>> nodes = null;

        /**
         * create a new tree with the representative and parent being the node.
         * the node is added to the internal set holding tree nodes too.
         * @param node 
         */
        public RootedTreeDisjointSet(DisjointSet2Node<T> node) {
            this.parent = node;
            this.repr = node;
            nodes = new HashSet<DisjointSet2Node<T>>();
            nodes.add(node);
        }

        public RootedTreeDisjointSet() {
            nodes = new HashSet<DisjointSet2Node<T>>();
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("nodes.size=");
            if (nodes != null) {
                sb.append(nodes.size());
            } else {
                sb.append("0");
            }
            sb.append(", [repr=");
            if (repr != null) {
                if (repr.member != null) {
                    sb.append(repr.member);
                } else {
                    sb.append(repr.hashCode());
                }
            }
            sb.append(", ");
            sb.append("parent=");
            if (parent != null) {
                if (parent.member != null) {
                    sb.append(parent.member);
                } else {
                    sb.append(parent.hashCode());
                }
            }
            sb.append(", ");
            sb.append("\n  nodes=");
            if (nodes != null) {
                for (DisjointSet2Node<T> node : nodes){
                    sb.append("\n   node=").append(node.toString());
                }
            }
            sb.append(", ");
            sb.append("] ");

            return sb.toString();
        }

    }

    /**
     * @return the trees
     */
    public Map<DisjointSet2Node<T>, RootedTreeDisjointSet<T>> getTrees() {
        return trees;
    }
}
