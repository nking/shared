package algorithms.tsp;

import algorithms.Permutations;
import algorithms.SubsetChooser;
import algorithms.misc.MiscMath0;
import gnu.trove.iterator.TLongDoubleIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.TLongDoubleMap;
import gnu.trove.map.hash.TLongDoubleHashMap;
import java.math.BigInteger;
import java.util.Arrays;

/**
 *
 * The dynamic approach needs to store the partial sums as the path lengthens
     (without the start node).
    
  data structure requirements for memo:
     - need set(key, value), getValue(key), and contains(key), and prefer O(1) or O(small const) for each.
     - a key for one partial sum would have to hold information on the nodes in
       the partial path and their order.
       could consider many ways to do that:
         The number of bits needed to represent the number of a node would be math.log(n-1)/math.log(2).
         Let w = math.log(n-1)/math.log(2).  then a key as a bit-string representing node and order
         would be a variable bit length from 3 nodes = 3*w bits in length up to (n-1)*w bits in length.
         For example, a 3-node path sum from nodes 2, 3, 1 would be
              2       3       1
           ------  ------  ------
           w bits  w bits  w bits

         Compression options for the key could be considered too.

       Note that the start node isn't part of those keys as it will be added into the path at
       evaluation time.
     - A worse case number of keys to store is the permutation of all but one node: n!/n.
     - A long bit-string is limited to 63 bits.
     - A java array length is limited to signed integer length, 1 &lt&lt 31 -1.
     - so to fit more than 1 &lt&lt 31-1 permutations would need to use the bit-string key
       as a concatenation of more than one path bit-string.
       the java BigInteger or the VeryLongBitString could do so.
       - also, the values for those keys would need to be referenced in an O(1) manner also if possible.
         one could use a 2 dimensional array where the 2nd dimension
         holds the individual values for the keys in the concatenated key.
     - would need a wrapper class to manage the integrity of the array of BigInteger or
       VeryLongBitString and the 2 dimensional array of path sums.
     - what amount of concatenation of keys would improve the number of storable
       path sums significantly?

 * @author nichole
 */
public abstract class AbstractTSP {

    protected double minCost = Double.POSITIVE_INFINITY;
    protected final TLongList minPath = new TLongArrayList();
    protected final int startNode = 0;
    protected final double[][] dist;
    protected final TLongDoubleMap memo;
    protected final double sentinel = Double.POSITIVE_INFINITY;
    protected final long totalNPerm;
    protected final long totalNSubSet;
    protected final long totalNSubSeq;
    protected final int w; // number of bits a city takes in a path where path is a bitstring of type long

    public AbstractTSP(double[][] dist) {
        
        this.dist = dist;
        int n = dist.length;
        long nPerm = MiscMath0.factorial(n); // max for n=13 for limit of array length
        totalNPerm = nPerm/n;
        
        //TODO: add in the number of permutations for those not in a 3-set, that is,
        //   the 2 node and 1-node permutations
        totalNSubSet = countTotalNumSubSetInvocations(n - 1); // max for n=338 for limit of array length
        totalNSubSeq = countTotalNumSubSeqInvocations(n - 1); 
        
        System.out.printf("nPerm=%d, totalNSubSet=%d  totalNSubSeq=%d\n", 
            totalNPerm, totalNSubSet, totalNSubSeq);
        
        int sz = (int)MiscMath0.computeNDivNMinusK(dist.length-1, 3);
        
        /*if (totalNSubSeq > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("this class can solve for 13 cities at most."
                + " one could design a datastructure to hold more than (1<<31)-1 items, "
                + " but the algorithms in this class have exponential runtime complexity,"
                + " so it would be better to choose a good approximate TSP.");
        }*/
        n = (int)totalNPerm;
                        
        w = (int)(Math.ceil(Math.log(dist.length)/Math.log(2)));
        
        memo = new TLongDoubleHashMap(sz);
        
        /*
        //find max n cities for limit of array length
        int nc = 100;
        long wc, nb;
        while (true) {
            long c = countTotalNumSubSeqInvocations(nc);
            if (c > Integer.MAX_VALUE) {
                break;
            }
            wc = (long)Math.ceil(Math.log(nc-1)/Math.log(2)); // 1 city encoding in bits
            nb = (Long.MAX_VALUE/wc); // the number of cities one could fit into a long if using wc
            System.out.printf("nc=%d, ns=%d wc=%d nb=%d ns*1e-9=%e\n", nc, c, wc, nb, c*1.e-9);
            nc = (int)Math.round(nc*1.5);
        }
        */
    }
    
    protected void reset() {
        minCost = sentinel;
        minPath.clear();
        memo.clear();
    }

    /**
     * roughly counting k-permutations for a dynamic approach where k is
     * increased by a factor of 2 each time and begins with k=3.
     * @param n
     * @return
     */
    protected static BigInteger count0(int n) {
        n = n - 1;
        BigInteger c1;
        BigInteger c2;
        BigInteger c3;
        BigInteger c0 = BigInteger.ZERO;
        int k = 3;
        while ((n - k) > k) {
            c1 = MiscMath0.computeNDivKTimesNMinusKBigInteger(n, k);
            c2 = MiscMath0.factorialBigInteger(k);
            c3 = c1.multiply(c2);
            c0 = c0.add(c3);
            k *= 2;
        }
        return c0;
    }
  
    // total number of subsetchooser invocations.  max n = 507 for count < Integer.MAX_VALUE
    protected long countTotalNumSubSeqInvocations(int n) {
        int k = 3;
        int n2 = n;
        long c = 0;
        while (n2 > k) {
            c += MiscMath0.computeNDivKTimesNMinusK(n2, k);
            n2 -= k;
        }
        return c;
    }

    // total number of k-permutations, that is, n!/(n-k)!.  using the subset chooser
    //   to find the sets of size 3, then permuting each of those into
    //   distinct ordered sequences
    protected long countTotalNumSubSetInvocations(int n) {
        int k = 3;
        int n2 = n;
        long c = 0;
        while (n2 > k) {
            c += MiscMath0.computeNDivNMinusK(n2, k);
            n2 -= k;
        }
        return c;
    }

    /**
     * given a pathNodeNumber such as 0 being the first node in the path bit-string,
     * uses path node bit-string length w to read the set bits and
     * return the base 10 number (which can be used with the distance matrix
     * or in writing the final solution of ordered nodes).
     * @param path
     * @param pathNodeNumber number of the node within the path.  NOTE: this excludes
     * start node 0 and end node 0 (pathNodeNumber=0 corresponds to the
     * 2nd node in the final solution for the completed path for the algorithm instance).
     * @return
     */
    protected int getBase10NodeIndex(final long pathNodeNumber, final long path) {
        // read bits pathNodeNumber*w to pathNodeNumber*w + w
        long sum = 0;
        long b = pathNodeNumber * w;
        long end = b + w;
        int s = 0;
        for (b = pathNodeNumber * w, s = 0; b < end; ++b, s++) {
            if ((path & (1L << b)) != 0) {
                // test bit b in path
                sum += (1L << s);
            }
        }
        return (int) sum;
    }

    /**
     * write the base 10 indexes s into a bit-string in the encoding used by the
     * memo.
     * @param s base 10 node indexes in the order to set into the bit-string path.
     * NOTE that s should not contain the startNode.
     * @return
     */
    protected long createThe3NodeBitstring(int[] s) {
        assert (s.length == 3);
        long path = concatenate(0, 0, s);
        return path;
    }

    /**
     * for the pathNodeNumber which is the order number of the node in the path
     * (e.g. second node in the path), set the node number to be the base10Node.
     * @param base10Node
     * @param path
     * @param pathNodeNumber number of the node within the path.  NOTE: this excludes
     * start node 0 and end node 0 (pathNodeNumber=0 corresponds to the
     * 2nd node in the final solution for the completed path for the algorithm instance).
     * @return
     */
    protected long setBits(final int base10Node, final long path, final int pathNodeNumber) {
        long bitstring = path; // 11 10 01
        long b = pathNodeNumber * w;
        long end = b + w;
        int b0;
        for (b = pathNodeNumber * w, b0 = 0; b < end; ++b, b0++) {
            if ((base10Node & (1L << b0)) != 0) {
                // test bit b0 in pathNodeNumber
                bitstring |= (1L << b); // set bit b in bitstring
            } else {
                bitstring &= ~(1L << b); //clear bit
            }
        }
        return bitstring;
    }

    /**
     * read the given bit-string encoded for use with memo, to find the
     * set bits and return the nodes as a base10 bit-string (without the path order information).
     * @param bitstring
     * @return
     */
    protected long findSetBitsBase10(long bitstring) {
        long base10nodes = 0;
        int i;
        int b;
        int node;
        int j;
        int bf = w * (dist.length - 1);
        for (b = 0; b < bf; b += w) {
            node = 0;
            for (i = b, j = 0; i < (b + w); ++i, ++j) {
                if ((bitstring & (1L << i)) != 0) {
                    node += (1 << j);
                }
            }
            if (node > 0) {
                base10nodes |= (1L << node);
            }
        }
        return base10nodes;
    }

    /**
     * read the given bit-string encoded for use with memo, to find the
     * unset bits and return the nodes as a base10 bit-string.
     * Note that bit 0 is not read as that is the implicit startNode
     * which is excluded from bit-string operations.
     * @param bitstring
     * @return
     */
    protected long findUnsetBitsBase10(long bitstring) {
        long base10nodesSet = findSetBitsBase10(bitstring);
        long base10NotSet = 0;
        // find bits not set from bit 1 to dist.length
        int b;
        for (b = 1; b < dist.length; ++b) {
            if ((base10nodesSet & (1 << b)) == 0) {
                base10NotSet |= (1L << b);
            }
        }
        return base10NotSet;
    }

    /**
     * read the given bit-string encoded for use with memo, to find the
     * unset bits and return the nodes as a base10 bit-string.
     * Note that bit 0 is not read as that is the implicit startNode
     * which is excluded from bit-string operations.
     * @param bitstring
     * @param out list to hold the output node indexes
     */
    protected void findUnsetBitsBase10(long bitstring, TIntList out) {
        out.clear();
        long base10nodesSet = findSetBitsBase10(bitstring);
        // find bits not set from bit 1 to dist.length
        int b;
        for (b = 1; b < dist.length; ++b) {
            if ((base10nodesSet & (1 << b)) == 0) {
                out.add(b);
            }
        }
    }

    /**
     * read the given bit-string encoded for use with memo, to find the
     * unset bits and return the nodes as a base10 bit-string that has lost
     * information about the path node order.
     * Note that bit 0 is not read as that is the implicit startNode
     * which is excluded from bit-string operations.
     * @param bitstring
     * @param out list to hold the output node indexes
     */
    protected void findSetBitsBase10(long bitstring, TIntList out) {
        out.clear();
        long base10nodesSet = findSetBitsBase10(bitstring);
        // find bits not set from bit 1 to dist.length
        int b;
        for (b = 1; b < dist.length; ++b) {
            if ((base10nodesSet & (1 << b)) != 0) {
                out.add(b);
            }
        }
    }

    /**
     * read path into base10 node numbers, preserving order of path
     * @param bitstring
     * @param out
     */
    protected void readPathIntoBase10(long bitstring, TIntList out) {
        out.clear();
        int i;
        int b;
        int node;
        int j;
        int bf = w * (dist.length - 1);
        for (b = 0; b < bf; b += w) {
            node = 0;
            for (i = b, j = 0; i < (b + w); ++i, ++j) {
                if ((bitstring & (1L << i)) != 0) {
                    node += (1 << j);
                }
            }
            if (node > 0) {
                out.add(node);
            }
        }
    }

    /**
     * @param path encoded bit-string of ordered path nodes used in the memo.
     * @param nPathNodesSet the number of nodes currently set in the path
     * @param base10Nodes base 10 node indexes in the order to set into the bit-string path.
     * NOTE that s should not contain the startNode.
     */
    protected long concatenate(long path, int nPathNodesSet, int[] base10Nodes) {
        assert (numberOfSetNodes(path) == nPathNodesSet);
        long path2 = path;
        int i;
        int si;
        for (i = 0; i < base10Nodes.length; ++i) {
            si = base10Nodes[i];
            assert (si != startNode);
            path2 = setBits(si, path2, nPathNodesSet + i);
        }
        return path2;
    }
    
    /**
     * initialize memo with permutations for all unset path nodes, where the
     * number of unset path nodes < 4.
     */
    protected void initNodePaths() {
        assert(memo.isEmpty());
        
        int nUnset = dist.length - 1;
        
        int i;
        int[] sel = new int[nUnset];
        for (i = 1; i <= sel.length; ++i) {
            sel[i-1] = i;
        }
        
        int nPerm = (int)MiscMath0.factorial(nUnset);
        
        final int[][] selPerm = new int[nPerm][nUnset];
        for (i = 0; i < selPerm.length; ++i) {
            selPerm[i] = new int[nUnset];
        }
        
        Permutations.permute(sel, selPerm);
        
        double sum;
        long path;
        int j, i0, i1;
        for (i = 0; i < selPerm.length; ++i) {
            sum = 0;
            
            //System.out.println("    selPerm=" + Arrays.toString(selPerm[i]));

            path = createThe3NodeBitstring(selPerm[i]);

            for (j = 1; j < selPerm[i].length; ++j) {
                i0 = selPerm[i][j - 1];
                i1 = selPerm[i][j];
                sum += dist[i0][i1];
            }
            memo.put(path, sum);
        }
    }

    /**
     * initialize memo with permutations for all subsets of 3 path nodes, where the
     * number of unset path nodes is > 3.
     */
    protected void init3NodePaths() {        
        assert(memo.isEmpty());
        
        int k = 3;
        final int[] sel = new int[k];
        final int[] sel2 = new int[k];
        int s, i;
        final int[][] selPerm = new int[6][k];
        for (i = 0; i < selPerm.length; ++i) {
            selPerm[i] = new int[k];
        }
        
        TIntList remaining = new TIntArrayList();
        findUnsetBitsBase10(0, remaining);
        //System.out.println("remaining unset=" + remaining.toString());

        int j, i0, i1;
        long path, sum;
        SubsetChooser chooser = new SubsetChooser(dist.length-1, k);
        while (true) {
            s = chooser.getNextSubset(sel);
            if (s == -1) {
                break;
            }
            //System.out.println("sel=" + Arrays.toString(sel));
            //transform sel to the bitstring unset indexes
            for (i = 0; i < k; ++i) {
                sel2[i] = remaining.get(sel[i]);
            }
            //System.out.println("    sel2=" + Arrays.toString(sel2));

            Permutations.permute(sel2, selPerm);
            
            for (i = 0; i < selPerm.length; ++i) {
                sum = 0;
                path = createThe3NodeBitstring(selPerm[i]);
                
                for (j = 1; j < k; ++j) {
                    i0 = selPerm[i][j-1];
                    i1 = selPerm[i][j];
                    sum += dist[i0][i1];
                }
                memo.put(path, sum);
            }
        }
    }

    protected void compareToMin(long path, double sum) {
        assert (numberOfSetNodes(path) == (dist.length - 1));

        int node1 = getBase10NodeIndex(0, path);
        int noden1 = getBase10NodeIndex(dist.length - 2, path);

        double ends = dist[startNode][node1] + dist[noden1][startNode];

        double sum2 = sum + ends;

        //debug
        /*TIntList p = new TIntArrayList();
        readPathIntoBase10(path, p);
        System.out.printf("final: bs=%s (%s) sum=%.2f sum2=%.2f, min=%.2f\n",
            Long.toBinaryString(path),
            Arrays.toString(p.toArray()), sum, sum2, minCost);
         */
        // end debug
        if (sum2 == minCost) {
            minPath.add(path);
        } else if (sum2 < minCost) {
            minCost = sum2;//3, 2, 4, 1, 5
            minPath.clear();
            minPath.add(path);
        }
    }
 
    public double getMinCost() {
        return this.minCost;
    }

    public int getNumberOfMinPaths() {
        return minPath.size();
    }

    public TIntList getMinPath(int idx) {
        long path = minPath.get(idx);
        TIntList out = new TIntArrayList();
        readPathIntoBase10(path, out);
        out.insert(0, startNode);
        out.add(startNode);
        return out;
    }

    /**
     *
     * @param path a memo bit-string
     * @return
     */
    protected int numberOfUnsetNodes(long path) {
        // composed of w-bit bit-strings
        int nn = (dist.length - 1) * w;
        // the first 64-((dist.length-1*w) can be discarded as unused space
        int nd = 64 - nn;
        int n0 = Long.numberOfLeadingZeros(path);
        assert (n0 >= nd);
        int nc = n0 - nd;
        int nUnset = nc / w;
        return nUnset;
    }

    /**
     *
     * @param path a memo bit-string
     * @return
     */
    protected int numberOfSetNodes(long path) {
        int nUnset = numberOfUnsetNodes(path);
        int nSet = (dist.length - 1) - nUnset;
        return nSet;
    }

    public int getMemoLength() {
        return memo.size();
    }
    
    protected void printMemo() {
        TLongDoubleIterator iter = memo.iterator();
        long bitstring;
        double sum;
        TIntList p = new TIntArrayList();
        
        for (int i = 0; i < memo.size(); ++i) {
            iter.advance();
            bitstring = iter.key();
            sum = iter.value();
            
            p.clear();
            readPathIntoBase10(bitstring, p);
            
            System.out.printf("memo: (%s) sum=%.2f min=%.2f\n",
                Arrays.toString(p.toArray()), sum, minCost);
        }
        System.out.flush();
    }
}
