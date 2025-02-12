package algorithms.combPerm;

import algorithms.misc.MiscMath0;
import gnu.trove.set.hash.TIntHashSet;

import java.math.BigInteger;
import java.util.*;

/**
Class to iterate over every combination of sub-sequences within n objects in an 
ordered manner.
<pre>
The number n is the number of distinct items.  
k is the size of the sub-sequence within the set of n numbers.
The number of sequences returned is 
    = n! /(k!*(n − k)!).

The class uses Gosper's hack as demonstrated by Eddie Kohler:
  https://read.seas.harvard.edu/~kohler/class/cs207-s12/lec12.html

This class was implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) Climb With Your Feet
     and using The MIT License (MIT)
     and moved to this project.
     
 Example use: 
      int n=7; int k=3;
      SubsetChooser chooser = new SubsetChooser(n, k);
      long s;
      while (true) {
          s = chooser.getNextSubset64Bitstring();
          if (s == -1) {
              break;
          }
       }
    results in the following 35 numbers
      7 (    111)
     11 (   1011)
     13 (   1101)
     14 (   1110)
     19 (  10011)
     21 (  10101)
     22 (  10110)
     25 (  11001)
     26 (  11010)  
     28 (  11100)  
     35 ( 100011)  
     37 ( 100101)  
     38 ( 100110) 
     41 ( 101001) 
     42 ( 101010) 
     44 ( 101100) 
     49 ( 110001) 
     50 ( 110010) 
     52 ( 110100)
     56 ( 111000)
     67 (1000011)
     69 (1000101)
     70 (1000110)
     73 (1001001)
     74 (1001010)
     76 (1001100)
     81 (1010001)
     82 (1010010)
     84 (1010100)
     88 (1011000)
     97 (1100001)
     98 (1100010)
     100 (1100100)
     104 (1101000)
     112 (1110000)
 </pre>
 * @author nichole
 */
public class SubsetChooser {
    
    private final int n;

    private final int k;

    private long x64;
    
    private final long highBit64;
    
    private BigInteger x;
    
    private final BigInteger highBit;

    private long count = 0;

    private final long np;

    /**
     * constructor with the number of indexes to choose from, n, and the size of
     * the subset, k.
     @param n the number of indexes that the selector will choose from
     @param k the subset size of selected indexes.  the maximum value currently
     accepted is 12.
     @throws ArithmeticException thrown when number of combinations is out of 
     range of type long
     */
    public SubsetChooser(int n, int k) {
        
        /*if (k > 13) {
            throw new IllegalArgumentException(
                "currently, class can only handle k < 13, but changes to accomodate larger could be made");
        }*/
        if (n < 1) {
            throw new IllegalArgumentException("n must be larger than 0");
        }
        if (k < 1) {
            throw new IllegalArgumentException("k must be larger than 0");
        }
        if (k > n) {
            throw new IllegalArgumentException("k must be less than or equal to n");
        }
        
        this.n = n;

        this.k = k;

        count = 1;

        // n!/(k!(n-k)!).  r.t.c. is O(k)
        np = MiscMath0.computeNDivKTimesNMinusK(n, k);

        if (n < 64) {
            
            highBit64 = 1L << n;

            x64 = (1L << k) - 1;

            x = null;
            
            highBit = null;
            
        } else {
            
            // x = BigInteger.ONE;
            // x = x.shiftLeft(k);
            // x = x.subtract(BigInteger.ONE);
            byte[] val = MiscMath0.writeToBigEndianBytes((1L << k) - 1);
            x = new BigInteger(val);
            
            BigInteger hb = BigInteger.ONE;
            highBit = hb.shiftLeft(n);  
            
            highBit64 = Long.MAX_VALUE;
            x64 = Long.MAX_VALUE;
        }
    }

    /**
     * given a constructed array, populates it with the next selected subset
     * of indexes and returns the number of values placed in the subset.
     * Returns a -1 when there are no more subsets to return;
     @param outputIndexes
     @return bitstring for next subset
     */
    public int getNextSubset(int[] outputIndexes) {

        if (outputIndexes == null || outputIndexes.length != k) {
            throw new IllegalArgumentException(
                "outputIndexes cannot be null and has to be size k");
        }
        
        if (n < 64) {
            return getNextSubset64(outputIndexes);
        } else {
            return getNextSubsetBigInteger(outputIndexes);
        }
        
    }
    
    /**
     * get the bitstring for the next subset
     @param outputIndexes
     @return bitstring for next subset
     */
    protected int getNextSubsetBigInteger(int[] outputIndexes) {

        if (outputIndexes == null || outputIndexes.length != k) {
            throw new IllegalArgumentException(
                "outputIndexes cannot be null and has to be size k");
        }
        
        if (count > np) {
            return -1;
        }

        int nValues = selectBigInteger(outputIndexes);

        x = nextSubsetBigInteger(x);

        return nValues;
    }

    /**
     * returns true if there is a next subset, that is, not done
     * with subsets
     * @return true if there is another subset.
     */
    public boolean hasNext() {
        return (count < np);
    }

    /**
     * returns the total number of subsets to expect.
     * @return total number of subsets to expect
     */
    public long getNp() {
        return np;
    }
    
    /**
     * given a constructed array, populates it with the next selected subset
     * of indexes and returns the number of values placed in the subset.
     * Returns a -1 when there are no more subsets to return;
     @param outputIndexes
     @return bitstring of the next subset
     */
    protected int getNextSubset64(int[] outputIndexes) {

        if (outputIndexes == null || outputIndexes.length != k) {
            throw new IllegalArgumentException(
                "outputIndexes cannot be null and has to be size k");
        }

        if (count > np) {
            return -1;
        }

        int nValues = select64(outputIndexes);

        x64 = nextSubset64(x64);

        return nValues;
    }
    
    /**
     * gets bitstring of the next subset
     @return bistring of the next subset
     */
    public long getNextSubset64Bitstring() {

        if (count > np) {
            return -1;
        }
        
        long r = x64;

        x64 = nextSubset64(r);

        return r;
    }

    /**
     * calculates bitstring of next subset
     @param x0 current set bitstring
     @return bitstring of next subset
     */
    private long nextSubset64(long x0) {

        long y = x0 & -x0;  // = the least significant one bit of x0
        long c = x0 + y; // set the next 0 bit that is higher than y

        // c^x0 isolates all the bits from y to c as all 1s.
        // y is a power of 2, so (c ^ x0)/y down shifts (c ^ x0) by LSB(x0).
        // then >>2 down shifts twice more to get the pattern to add to c
        x0 = c + (((c ^ x0) / y) >> 2);

        count++;

        return x0;
    }

    /**
     * calculates bitstring of next subset
     @param x0 current set bitstring
     @return bitstring of next subset
     */
    private BigInteger nextSubsetBigInteger(BigInteger x0) {

        BigInteger y = x0.and(x0.negate()); // = the least significant one bit of x
        BigInteger c = x0.add(y);
        
        //x0 = c + (((c ^ x0) / y) >> 2);
        BigInteger tmp = c.xor(x0).divide(y).shiftRight(2);
        
        x0 = c.add(tmp);

        count++;

        return x0;
    }

    /**
     * populate the array selected with the next subset
     @param selected
     @return the number of set bits in selected
     */
    protected int select64(int[] selected) {

        // interpret the bit string x:  1 is 'selected' and 0 is not

        /*
        String str = Long.toBinaryString(x);
        while (str.length() < n) {
            str = "0" + str;
        }
        System.out.format("%d\t%10s
", x, str);
        */

        int nBits = 0;
        int nOneBits = 0;
        long xp = x64;
        while (xp > 0) {
            if ((xp & 1) == 1) {
                int idx2 = n - 1 - nBits;
                selected[nOneBits] = idx2;
                nOneBits++;
            }
            xp = xp >> 1;
            nBits++;
        }

        return nOneBits;
    }

    /**
     * populate the array selected with the next subset
     @param selected
     @return the number of set bits in selected
     */
    protected int selectBigInteger(int[] selected) {

        // interpret the bit string x:  1 is 'selected' and 0 is not

        /*
        String str = Long.toBinaryString(x);
        while (str.length() < n) {
            str = "0" + str;
        }
        System.out.format("%d\t%10s
", x, str);
        */

        int nBits = 0;
        int nOneBits = 0;
        BigInteger xp = new BigInteger(x.toByteArray());
        while (xp.compareTo(BigInteger.ZERO) > 0) {
            if (xp.and(BigInteger.ONE).equals(BigInteger.ONE)) {
                int idx2 = n - 1 - nBits;
                selected[nOneBits] = idx2;
                nOneBits++;
            }
            xp = xp.shiftRight(1);
            nBits++;
        }

        return nOneBits;
    }

    /**
     * calculate all subsequences of size k from size a.
     * the r.t.c. on the order of O(n!/(k!*(n-k)!)
     * @param a array
     * @param k subsequence size
     * @return the subsequences of size k of a
     */
    public static List<int[]> calcSubSets(int[] a, int k) {
        int n = a.length;
        if (n < 1) {
            throw new IllegalArgumentException("n must be larger than 0");
        }
        if (k < 1) {
            throw new IllegalArgumentException("k must be larger than 0");
        }
        if (k > n) {
            throw new IllegalArgumentException("k must be less than or equal to n");
        }
        // n!/(k!(n-k)!) number of subsequences
        long nck = MiscMath0.computeNDivKTimesNMinusK(n, k);
        if (nck > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("the number of combinations is larger than max length of an array," +
                    "so this algorithm needs to be adjusted to return one element at a time");
        }
        List<int[]> out = new ArrayList<>();

        // memo holds the first index of a subset.
        // when DFS has completed for the first index, it's stored in memo to avoid repeating the work
        Set<Integer> memo = new HashSet<>();

        nIter = 0;
        recurseSet(a, new int[k], 0, 0, out, memo, new HashSet<>());

        return out;
    }

    static int nIter = 0;
    private static void recurseSet(int[] a, int[] s, int sIdx, int iIdx, List<int[]> out, Set<Integer> memo,
                                   Set<Integer> drawn) {
        ++nIter;
        if (sIdx == s.length || iIdx == a.length) {
            if (drawn.size() == s.length) {
                out.add(Arrays.copyOf(s, s.length));
            }
            return;
        }

        if (memo.contains(iIdx)) {
            return;
        }

        // include i:
        drawn.add(iIdx);
        s[sIdx] = a[iIdx];
        recurseSet(a, s, sIdx + 1, iIdx + 1, out, memo, drawn);

        // exclude i
        drawn.remove(iIdx);
        recurseSet(a, s, sIdx, iIdx + 1, out, memo, drawn);

        if (sIdx == 0) {
            // store in memo after returned from DFS for the first index of s
            memo.add(iIdx);
        }

    }
}
