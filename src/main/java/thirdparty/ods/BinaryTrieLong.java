package thirdparty.ods;

/*
The class is adapted from the open datastructures source code
http://opendatastructures.org/ods-java.pdf

"The source code available there is released under a Creative Commons
Attribution license, meaning that anyone is free to share: to copy, distribute
and transmit the work; and to remix: to adapt the work, including
the right to make commercial use of the work. The only condition on
these rights is attribution: you must acknowledge that the derived work
contains code and/or text from opendatastructures.org.
http://github.com/patmorin/ods

Edits were made to the code in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet, The MIT License (MIT)
and then moved to this project

*/
import algorithms.util.ObjectSpaceEstimator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.Deque;

/**
 *
 * @author nichole
 @param <S> parameter type for node
 @param <T> parameter type for value held by node
 */
@SuppressWarnings({"unchecked"})
public class BinaryTrieLong<S extends BinaryTrieNode<T>, T> 
    implements SSet<T> {
	
    /**
     *
     */
    protected static final int prev = 0;

    /**
     *
     */
    protected static final int next = 1;

    /**
     *
     */
    protected static final int left = 0;

    /**
     *
     */
    protected static final int right = 1;
	
    /**
     *
     */
    protected int w = 62;

    /**
     *
     */
    protected final long maxC;
    
	/**
	* The root node
	*/
	protected final S r;
	
	/**
	 * The dummy node in the doubly-linked list
	 */
	protected final S dummy;
	
	/**
	 * For converting elements of type T into integers
	 */
	protected final Longizer<T> it;

	/**
	 * The number of elements stored in the trie
	 */
	int n;

	/**
	 * To make a node factory
	 */
	protected final S sampleNode;
	
	/**
	 * Allocate a new node.  if S is an extension
     * of BinaryTrieNode,
     * make sure the class is public and if it's an
     * inner class, it has to be static too.
	 @return
	 */
	protected S newNode() {
		try {
            S u = (S)sampleNode.getClass().getDeclaredConstructor().newInstance();
			u.parent = u.child[0] = u.child[1] = null;
			return u;
		} catch (Exception e) {
            throw new UnsupportedOperationException(
               "sampleNode constructor is not reachable "
               + "as a public merhod with a no arguments constructor");
		}
	}

    /**
     * 
     @param sampleNode a node instance that is used for 
     * its class type when creating other nodes, such as the 
     * root and linked-list sentinel nodes.
     @param it class to provide the inner node key which
     * prefixes are extracted from
     */
	public BinaryTrieLong(S sampleNode, Longizer<T> it) {
		this.sampleNode = sampleNode;
		this.dummy = newNode();
		dummy.child[prev] = dummy.child[next] = dummy;
		this.r = newNode();
		r.jump = dummy;
		this.it = it;
		n = 0;
        maxC = (1L << w) - 1;
	}
	
    /**
     * constructor that accepts a sample node type for
     * internal nodes and accepts the maximum number of bits
     * that a value added will have.
     @param sampleNode a node instance that is used for 
     * its class type when creating other nodes, such as the 
     * root and linked-list sentinel nodes.
     @param it class to provide inner node keys which the
     * prefixes are extracted from
     @param maxNumBits maximum number of bits a value that is
     * added to the trie will have when it is known to be less than
     * 32 (else, can use default constructor);
     */
	public BinaryTrieLong(S sampleNode, Longizer<T> it,
        int maxNumBits) {
        if (maxNumBits < 63 && maxNumBits > 1) {
            this.w = maxNumBits;
        } else {
            throw new IllegalStateException("maxNumBits "
                + " should be greater than 1 and less than 63");
        }
        maxC = (1L << w) - 1;
		this.sampleNode = sampleNode;
		this.dummy = newNode();
		dummy.child[prev] = dummy.child[next] = dummy;
		this.r = newNode();
		r.jump = dummy;
		this.it = it;
		n = 0;
	}
	
    /**
     * add x to the tree.  runtime complexity is O(w)
     * where w is the number of
     * bits set in the constructor, else is 32.
     @param x
     @return 
     */
    @Override
	public boolean add(T x) {
		int i, c = 0;
        long ix = it.longValue(x);
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value can add is " + maxC);
        }
		BinaryTrieNode<T> u = r;
		// 1 - search for ix until falling out of the trie
		for (i = 0; i < w; i++) {
			c = (int)((ix >>> w-i-1) & 1);
            assert(c == 0 || c == 1);
			if (u.child[c] == null) break;
			u = u.child[c];
		}
		if (i == w) return false; // already contains x - abort
		BinaryTrieNode<T> pred = (c == right) ? u.jump : u.jump.child[0];
		u.jump = null;  // u will have two children shortly
		// 2 - add path to ix
		for (; i < w; i++) {
			c = (int)((ix >>> w-i-1) & 1);
            assert(c == 0 || c == 1);
			u.child[c] = newNode();
			u.child[c].parent = u;
			u = u.child[c];
		}
		u.x = x;
		// 3 - add u to linked list
		u.child[prev] = pred;
		u.child[next] = pred.child[next];
		u.child[prev].child[next] = u;
		u.child[next].child[prev] = u;
		// 4 - walk back up, updating jump pointers
		BinaryTrieNode<T> v = u.parent;
		while (v != null) {
			if ((v.child[left] == null 
	        	&& (v.jump == null ||
                it.longValue(v.jump.x) > ix))
			|| (v.child[right] == null 
	    		&& (v.jump == null || 
                (v.jump.x != null && it.longValue(v.jump.x) < ix))
                )) {
				v.jump = u;
            }
			v = v.parent;
		}
        
		n++;
        
		return true;
	}

    /**
     * find the node if it exists in the trie.
     * runtime complexity is O(w)
     * where w is the number of
     * bits set in the constructor, else is 32.
     @param x
     @return 
     */
    @Override
    public T find(T x) {
        long ix = it.longValue(x);
        T v = findValue(ix);
        if (v == null) {
            return v;
        } else if (it.longValue(v) == ix) {
            return v;
        }
        return null;
    }
    
    /**
     *
     @param ix
     @return
     */
    protected T findValue(long ix) {
		int i, c = 0;
		if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value can find is " + maxC);
        }
        BinaryTrieNode<T> u = r;
		for (i = 0; i < w; i++) {
			c = (int)((ix >>> w-i-1) & 1);
            assert(c == 0 || c == 1);
			if (u.child[c] == null) break;
			u = u.child[c];
		}
		if (i == w) return u.x;  // found it
        if (c == 0) {
            u = u.jump;
        } else {
            u = u.jump.child[next]; 
        }
		return u == dummy ? null : u.x;
	}

    /**
     * remove node from the trie.
     * runtime complexity is O(w) 
     * where w is the number of
     * bits set in the constructor, else is 32.
     @param x
     @return 
     */
    @Override
	public boolean remove(T x) {
		// 1 - find leaf, u, containing x
		int i, c;
        long ix = it.longValue(x);
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value can remove is " + maxC);
        }
		BinaryTrieNode<T> u = r;
		for (i = 0; i < w; i++) {
			c = (int)((ix >>> w-i-1) & 1);
            assert(c == 0 || c == 1);
			if (u.child[c] == null) return false;
			u = u.child[c];
		}
		// 2 - remove u from linked list
		u.child[prev].child[next] = u.child[next];
		u.child[next].child[prev] = u.child[prev];
		BinaryTrieNode<T> v = u;
		// 3 - delete nodes on path to u
		for (i = w-1; i >= 0; i--) {
            c = (int)((ix >>> w-i-1) & 1);
            assert(c == 0 || c == 1);
			v = v.parent;
			v.child[c] = null;
			if (v.child[1-c] != null) break;
		}
		// 4 - update jump pointers
		c = (int)((ix >>> w-i-1) & 1);
        assert(c == 0 || c == 1);
		v.jump = u.child[1-c];
		v = v.parent;
		i--;
		for (; i >= 0; i--) {
			c = (int)((ix >>> w-i-1) & 1);
            assert(c == 0 || c == 1);
			if (v.jump == u) {
				v.jump = u.child[1-c];
            }
			v = v.parent;
		}
		n--;
		return true;
	}
    
    /**
     *
     @param ix
     @return
     */
    protected T predecessor(long ix) {
        S q = predecessorNode(ix);
        if (q != null) {
            return q.x;
        }
        return null;
    }

    /**
	 * Find the key of the node that contains the predecessor of x.
	 * runtime complexity is O(w) 
     * where w is the number of
     * bits set in the constructor, else is 32.
     @param x
	 @return The node before the node that contains x w.r.t. 
     * nodes in the internal the linked list.
	 */
    @Override
	public T predecessor(T x) {
        S q = predecessorNode(x);
        if (q != null) {
            return q.x;
        }
        return null;
    }
    
    /**
     *
     @param x
     @return
     */
    protected S predecessorNode(T x) {
		int i, c = 0;
        long ix = it.longValue(x);
        return predecessorNode(ix); 
   }
    
    /**
     *
     @param ix
     @return
     */
    protected S predecessorNode(long ix) {
		int i, c = 0;
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value argument is " + maxC);
        }
		BinaryTrieNode<T> u = r;
		for (i = 0; i < w; i++) {
			c = (int)((ix >>> w-i-1) & 1);
            assert(c == 0 || c == 1);
			if (u.child[c] == null) break;
			u = u.child[c];
		}
		BinaryTrieNode<T> pred;
		if (i == w) pred = u.child[prev];
		else pred = (c == 1) ? u.jump : u.jump.child[0]; 
		return (pred != null) ? (S)pred : null;
	}
    
    /**
     *
     @param ix
     @return
     */
    protected T successor(long ix) {
        S q = successorNode(ix);
        if (q != null) {
            return q.x;
        }
        return null;
    }
    
    /**
	 * Find the key of the node that contains the successor of x.
	 * runtime complexity is O(w)
     * where w is the number of
     * bits set in the constructor, else is 32.
     @param x
	 @return The key of the node before the node that contains x w.r.t. 
     * nodes in the internal the linked list.
	 */
    @Override
	public T successor(T x) {
        S q = successorNode(x);
        if (q != null) {
            return q.x;
        }
        return null;
    }
    
    /**
     *
     @param x
     @return
     */
    protected S successorNode(T x) {
		int i, c = 0;
        long ix = it.longValue(x);
        return successorNode(ix);
    }
    
    /**
     *
     @param ix
     @return
     */
    protected S successorNode(long ix) {
		int i, c = 0;
        if (ix > maxC) {
            throw new IllegalArgumentException("w=" + w
               + " so max value argument is " + maxC);
        }
		BinaryTrieNode<T> u = r;
		for (i = 0; i < w; i++) {
			c = (int)((ix >>> w-i-1) & 1);
            assert(c == 0 || c == 1);
			if (u.child[c] == null) break;
			u = u.child[c];
		}
		BinaryTrieNode<T> successor;
		if (i == w) successor = u.child[next];
		else successor = (c == 0) ? u.jump : u.jump.child[1]; 
		return (successor != null) ? (S)successor : null;
	}

    /**
     * find the maximum key within the nodes. 
     * runtime complexity is O(w) where w is the number of
     * bits set in the constructor, else is 32.
     @return 
     */
    @Override
    public T maximum() {
        // O(w)
        T q = findValue(maxC);
        if (q != null) {
            return q;
        }
        return predecessor(maxC);
    }
    
    /**
     * find the minimum key within the nodes. 
     * runtime complexity is O(w) where w is the maximum
     * word size set in the constructor, else 32.
     @return 
     */
    @Override
    public T minimum() {
        // O(w)
        T q = findValue(0);
        if (q != null) {
            return q;
        }
        return successor(0);
    }
    
	public int size() {
		return n;
	}
    
    /**
     *
     */
    public void clear() {
		n = 0;
		r.child[0] = r.child[1] = null;
		r.jump = dummy;
		dummy.child[0] = dummy.child[1] = dummy;
	}
	
    /**
     * print out the dummy link keys then print the
     * root tree in level  traversal
     */
    void debugNodes() {
        TIntSet dummyHashCodes = new TIntHashSet();
        S node = dummy;
        //System.out.println("dummy.hashCode=" + dummy.hashCode());
        System.out.print("\ndummy=");
        do {
            int dhc = node.hashCode();
            System.out.print(node.x + ", ");
            dummyHashCodes.add(dhc);
            node = (S)node.child[1];
        } while (!node.equals(dummy));
        
        System.out.println();
        
        if (r == null) {
            return;
        }
        
        int dhc = dummy.hashCode();
        
        Deque<S> q0 = new ArrayDeque<S>();
        Deque<S> q1 = new ArrayDeque<S>();
        q0.offer(r);
        
        int level = 0;
        int count = 0;
        boolean skip = true;
        while(!q0.isEmpty() && level < w) {
            System.out.println("count=" + level);
            while(!q0.isEmpty()) {
                node = q0.poll();
                if (!node.equals(r)) {
                    System.out.println(node.toString2());
                    count++;
                }
                if (node.child[0] != null) {
                    int hc = node.child[0].hashCode();
                    if (level < w) {
                        q1.offer((S) node.child[0]);
                    }
                }
                if (node.child[1] != null) {
                    int hc = node.child[1].hashCode();
                    if (level < w) {
                        q1.offer((S) node.child[1]);
                    }
                }
            }
            if (!skip) {
                level++;
            } else {
                skip = false;
            }
            q0.addAll(q1);
            q1.clear();
        }
        System.out.println("number of trie nodes=" + count);
    }
    
    /**
     * NOTE: there are prefix entries in the trie, created as needed. total
     * number of trie prefix nodes for sequential data is 2 * n + padding to
     * next power of 2. The number of prefix nodes is due to the pattern of
     * numbers, so not predictable, but a few tests show range of factor 2 to 5
     * times the number of added nodes.
     * A factor of 5 is used here.
     @param numberOfEntries
     @return 
     */
    public static long estimateSizeOfTriePrefixNodes(int numberOfEntries) {
        
        long factor = 5;
        
        long nodeSz = BinaryTrieNode.estimateSizeOnHeap();
        
        long total = factor * numberOfEntries * nodeSz;
        
        return total;
    }
    
    /**
     * estimate the size of an instance of this class on the heap with
     * n number of inserted entries.
     * NOTE there are prefix entries in the trie, created as
       needed and the separate method should be
       used for those: estimateSizeOfTrieNodes()
        
     @param numberOfEntries
     @return 
     */
    public static long estimateSizeOnHeap(int numberOfEntries) {
       
        long nodeSz = BinaryTrieNode.estimateSizeOnHeap();
        
        long total = numberOfEntries*nodeSz;
        
        // one more for this instance superclass
        total += nodeSz;
       
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNIntFields(6);
        est.setNObjRefsFields(4);
        est.setNLongFields(1);
       
        total += est.estimateSizeOnHeap();
        
        return total;
    }
    
}
