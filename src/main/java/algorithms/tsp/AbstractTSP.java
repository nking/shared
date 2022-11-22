package algorithms.tsp;

import algorithms.PermutationsWithAwait;
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
import java.util.Stack;

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

 * @author nichole
 */
public abstract class AbstractTSP {

    protected double minCost = Double.POSITIVE_INFINITY;
    protected final TLongList minPath = new TLongArrayList();
    protected final int startNode = 0;
    protected final double[][] dist;
    protected final TLongDoubleMap memo;
    protected final double sentinel = Double.POSITIVE_INFINITY;
    protected final BigInteger totalNPerm;
    protected final long totalNSubSet;
    protected final long totalNSubSeq;
    protected final int w; // number of bits a city takes in a path where path is a bitstring of type long

    public AbstractTSP(double[][] dist) {
        if (dist.length < 3) {
            throw new IllegalArgumentException("dist.length must be >= 3");
        }
        
        this.dist = dist;
        int n = dist.length;
        
        BigInteger nPerm = MiscMath0.factorialBigInteger(n); // max for n=13 for limit of array length
        totalNPerm = nPerm.divide(new BigInteger(Integer.toString(n)));

        //TODO: add in the number of permutations for those not in a 3-set, that is,
        //   the 2 node and 1-node permutations
        totalNSubSet = countTotalNumSubSetInvocations(n - 1); // max for n=338 for limit of array length
        totalNSubSeq = countTotalNumSubSeqInvocations(n - 1); 
        
        System.out.printf("nPerm (w/o 1st node)=%s, totalNSubSet=%d  totalNSubSeq=%d\n",
            totalNPerm.toString(), totalNSubSet, totalNSubSeq);
        
        int sz = (int)MiscMath0.computeNDivNMinusK(dist.length-1, 3);
                       
        w = (int)(Math.ceil(Math.log(dist.length)/Math.log(2)));
        
        memo = new TLongDoubleHashMap(sz);
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
    protected static BigInteger count3(int n) {
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
    
    /**
     * roughly counting k-permutations for a dynamic approach where k is
     * increased by a factor of 2 each time and begins with k=3.
     * @param n
     * @return
     */
    protected static BigInteger count2(int n) {
        n = n - 1;
        BigInteger c1;
        BigInteger c2;
        BigInteger c3;
        BigInteger c0 = BigInteger.ZERO;
        int k = 2;
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
     * @return the base10 node number for pathNodeNumber >= 0.
     * if pathNodeNumber is less than 0, -1 is returned.
     */
    protected int getBase10NodeIndex(final long pathNodeNumber, final long path) {
        if (pathNodeNumber < 0) {
            return -1;
        }
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
    protected long createAMemoNodeBitstring(int[] s) {
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
//            assert (si != startNode);
            path2 = setBits(si, path2, nPathNodesSet + i);
        }
        return path2;
    }
    
    /**
     * initialize memo with permutations for all unset path nodes, where the
     * number of unset path nodes < 4.
     */
    protected void initNodePaths() throws InterruptedException {
        assert(memo.isEmpty());
        
        int nNotSet = dist.length - 1;
        
        int i;
        int[] sel = new int[nNotSet];
        for (i = 1; i <= sel.length; ++i) {
            sel[i-1] = i;
        }
        
        TIntList nodes = new TIntArrayList(sel);
        
        long bitstring2 = 0; 
        double sum2 = 0;
        int nNodesSet = 0;
        Stack<StackP> stack = null;
        boolean storeInMemo = true;
        
        createAndStackPermutations(bitstring2, sum2, nNodesSet, 
            nodes, stack, storeInMemo);
    }

    /**
     * initialize memo with permutations for all subsets of 4 path nodes, where the
     * number of unset path nodes is > 4.
     */
    protected void init4NodePaths() throws InterruptedException {  
        initKNodePaths(4);
    }
    
    /**
     * initialize memo with permutations for all subsets of 3 path nodes, where the
     * number of unset path nodes is > 3.
     */
    protected void init3NodePaths() throws InterruptedException {  
        initKNodePaths(3);
    }
    
    /**
     * initialize memo with permutations for all subsets of 3 path nodes, where the
     * number of unset path nodes is > 3.
     * @param k
     */
    protected void initKNodePaths(final int k) throws InterruptedException {        
        assert(memo.isEmpty());
        
        int nNodesSet = 0;
        Stack<StackP> stack = null;
        boolean storeInMemo = true;
        double cost = 0;
        createAndStackSubsetPermutations(0, cost, nNodesSet, k, stack, storeInMemo);
    }
    
    /**
     * 
     * @param bitstring bit-string of ordered path nodes in format for memo key
     * @param sum the cost of the ordered path thus far
     * @param nNodesSet the number of nodes set in bitstring
     * @param k the length of subsets to choose from the unset bits of bitstring
     * @param stack can be null.  if not null, the permuted path and its cost
     * are pushed onto the stack.
     * @param storeInMemo if true, the permuted path and sum are stored in memo
     */
    protected void createAndStackSubsetPermutations(long bitstring, double sum, int nNodesSet, int k,
        Stack<StackP> stack, boolean storeInMemo) throws InterruptedException {
        
        TIntList remaining = new TIntArrayList();
        findUnsetBitsBase10(bitstring, remaining);
        //System.out.println("remaining unset=" + remaining.toString());
   
        int nPerm = (int)MiscMath0.factorial(k);
        
        final int[] sel = new int[k];
        final int[] sel2 = new int[k];
        int s, i;
        int[] selPerm = new int[k];
        
        int lastNode = getBase10NodeIndex(nNodesSet-1, bitstring);
        
        int j, i0, i1;
        long permi, path2;
        double sum2;
        SubsetChooser chooser = new SubsetChooser(remaining.size(), k);
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

            PermutationsWithAwait permutations = new PermutationsWithAwait(sel2);

            //for (i = 0; i < nPerm; ++i) {
            while (permutations.hasNext()) {
                permutations.getNext(selPerm);
                sum2 = sum;
                if (lastNode >= 0) {
                    sum2 += dist[lastNode][selPerm[0]];
                }
                
                permi = createAMemoNodeBitstring(selPerm);
                if (memo.containsKey(permi)) {
                    sum2 += memo.get(permi);
                } else {
                    // add each edge
                    for (j = 1; j < selPerm.length; ++j) {
                        i0 = selPerm[j - 1];
                        i1 = selPerm[j];
                        sum2 += dist[i0][i1];
                    }
                }
                
                path2 = concatenate(bitstring, nNodesSet, selPerm);
                
                if (storeInMemo) {
                    memo.put(path2, sum2);
                }
                if (stack != null) {
                    stack.add(new StackP(path2, sum2, remaining.size() - k));
                }
            }
        }
    }

    protected void compareToMin(long path, double sum) {
//        assert (numberOfSetNodes(path) == (dist.length - 1));

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
            minCost = sum2;
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
        System.out.printf("memo.size()=%d\n", memo.size());
        System.out.flush();
    }
    
    /**
     * 
     * @param bitstring2 the ordered path nodes in format for memo keys
     * @param sum2 the cost of the path of the ordered nodes in bitstring2
     * @param nNodesSet2 the number of nodes set in bitstring2
     * @param nodes nodes to permute
     * @param stack can be null
     * @param storeInMemo
     * @throws InterruptedException 
     */
    protected void createAndStackPermutations(long bitstring2, double sum2, 
        int nNodesSet2, TIntList nodes, Stack<StackP> stack, boolean storeInMemo) 
        throws InterruptedException {
        
        int nNodes = nodes.size();
        //int nBitsSet2 = (dist.length - 1) - nNodesRemaining2;
        //int nBitsSet2 = this.numberOfSetNodes(bitstring2);
                
        // note that this will be -1 if (nNodesSet2-1) < 0
        // it's the last node set into bitstring2
        int lastNode = getBase10NodeIndex(nNodesSet2-1, bitstring2);

        int np = (int)MiscMath0.factorial(nNodes);

        int[] selPerm = new int[nNodes];
            
        PermutationsWithAwait permutations = new PermutationsWithAwait(nodes.toArray());
        
        long permi, path3;
        double sum3;
        int i, j, i0, i1;
        

        for (i = 0; i < np; ++i) {

            permutations.getNext(selPerm);

            sum3 = sum2;
            if (lastNode >= 0) {
                sum3 += dist[lastNode][selPerm[0]];
            }
            
            permi = createAMemoNodeBitstring(selPerm);
            if (memo.containsKey(permi)) {                            
                sum3 += memo.get(permi);
            } else {
                // add each edge
                for (j = 1; j < selPerm.length; ++j) {
                    i0 = selPerm[j - 1];
                    i1 = selPerm[j];
                    sum3 += dist[i0][i1];
                }
            }
            
            path3 = concatenate(bitstring2, nNodesSet2, selPerm);
            
            if (storeInMemo) {
                memo.put(path3, sum3);
            }

            if (stack != null) {
                int nRemaining = (dist.length - 1) - (nNodesSet2 + nodes.size());
                stack.add(new StackP(path3, sum3, nRemaining));
            }
        }
    }
    
    protected static class StackP {
        long bitstring;
        double sum;
        int nNodesRemaining;
        public StackP(long path, double cost, int nRemaining) {
            this.bitstring = path;
            this.sum = cost;
            this.nNodesRemaining = nRemaining;
        }
    }
    
}
