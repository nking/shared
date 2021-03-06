package algorithms;

import java.util.Arrays;

/**
 * 
 * a modifiable holder for a very long bit string of fixed size nbits. 
 * the internal array type is long which means each item is a 64 bit holder
 * (128 bits on the stack of a 64 bit platform).  
 * the number of items in the internal long array is limited to 31 bits
 * (positive portion of signed 32 bit integer), 
 * so that's a capacity of 64 * 31 bits, if the jvm has enough memory for that.
 * 
 * first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.
  
 * @author nichole
 */
public final class VeryLongBitString {
    
    private final long[] bitstring;
    
    protected final long nBits;
    
    protected static final long ITEM_BIT_LENGTH = 64L;
    
    protected static final long N_MAX_BITS = Integer.MAX_VALUE * ITEM_BIT_LENGTH;
    
    protected final long capacityBits;
    
    private long nSetBits = 0;
   
    public VeryLongBitString(long nBits) {
        
        if (nBits > N_MAX_BITS) {
            throw new IllegalArgumentException("cannot hold more than " + N_MAX_BITS + " bits");
        }
        
        int nElements = getRowNumber(nBits) + 1;
        
        bitstring = new long[nElements];
        
        capacityBits = nElements * ITEM_BIT_LENGTH;
        
        this.nBits = nBits;
    }
    
    protected VeryLongBitString(long[] bitstrings, long nBits, long 
        nSetBits) {
        
        if (nBits > N_MAX_BITS) {
            throw new IllegalArgumentException("cannot hold more than " + N_MAX_BITS + " bits");
        }
        
        bitstring = Arrays.copyOf(bitstrings, bitstrings.length);
        
        capacityBits = bitstring.length * ITEM_BIT_LENGTH;
        
        this.nBits = nBits;
        
        this.nSetBits = nSetBits;
    }
    
    public long getCapacity() {
        return capacityBits;
    }
    
    public void setBit(long nthBit) {
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        assert(bitIdx < 64);
        assert(bitIdx >= 0);
        
        // test bit is not set
        if ((bitstring[idx] & (1L << bitIdx)) == 0) {
          
            // set bit
            bitstring[idx] |= (1L << bitIdx);
                        
            nSetBits++;
        }
    }
    
    public void clearBit(long nthBit) {
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        // test bit
        if ((bitstring[idx] & (1L << bitIdx)) != 0) {
        
            // clear bit
            bitstring[idx] &= ~(1L << bitIdx);
                        
            nSetBits--;
        }
    }
    
    public void toggleBit(long nthBit) {
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        // replace toggle so can keep track of nSetBits
        //bitstring[idx] ^= (1L << bitIdx);
        
        if ((bitstring[idx] & (1L << bitIdx)) == 0) {
            
            // set bit
            bitstring[idx] |= (1L << bitIdx);
            
            nSetBits++;
            
        } else {
            
            // clear bit
            bitstring[idx] &= ~(1L << bitIdx);
            
            nSetBits--;
        }
    }
    
    public long getNSetBits() {
        return nSetBits;
    }
    
    public long getInstantiatedBitSize() {
        return nBits;
    }
    
    public boolean isSet(long nthBit) {
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        return ((bitstring[idx] & (1L << bitIdx)) != 0);
    }
    
    public boolean isNotSet(long nthBit) {
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        return ((bitstring[idx] & (1L << bitIdx)) == 0);
    }
         
    protected int getRowNumber(long n) {
        
        int nthElement = (int)(n/ITEM_BIT_LENGTH);
        
        return nthElement;
    }
    
    /**
     * given the arrayIdx, calculate the bit position within that item
     * for the very large bitstring position nthBit.  Note that if the
     * result is out of bounds, -1 is returned and the invoker should
     * handle that.
     * 
     * @param nthBit
     * @param arrayIdx
     * @return 
     */
    int getBitIdx(long nthBit, int arrayIdx) {
        
        int bitIdx = (int)(nthBit - (arrayIdx * ITEM_BIT_LENGTH));
        
        if ((bitIdx < 0) || (bitIdx > (ITEM_BIT_LENGTH - 1))) {
            return -1;
        }
        
        return bitIdx;
    }
    
    public void clearAllBits() {
        
        Arrays.fill(bitstring, 0);
        
        nSetBits = 0;
    }
    
    public VeryLongBitString copy() {
        
        VeryLongBitString c = new VeryLongBitString(bitstring, 
            nBits, nSetBits);
        
        return c;
    }
    
    public void resetAllTo(VeryLongBitString other) {
        if (other.nBits != nBits) {
            throw new IllegalArgumentException("nBits must be the same in both to use this method");
        }
        System.arraycopy(other.bitstring, 0, bitstring, 0, other.bitstring.length);
        
        this.nSetBits = other.nSetBits;
    }
    
    /**
     * get the '01...' bitstring representation of this object
     * @return 
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = (bitstring.length - 1); i > -1 ; i--) {
            sb.append(Long.toBinaryString(bitstring[i]));
        }
        return sb.toString();
    }
    
    protected void recountNSetBits() {
        nSetBits = 0;
        for (int i = 0; i < bitstring.length; ++i) {
            nSetBits += Long.bitCount(bitstring[i]);
        }
    }
    
    /**
     * get a list of the bit numbers that are set.
     * @return 
     */
    public int[] getSetBits() {
        
        int n = 0;
        for (int i = 0; i < bitstring.length; ++i) {
            n += Long.bitCount(bitstring[i]);
        }
        assert(n == nSetBits);
        
        int[] setBits = new int[n];
        int n2 = 0;
        for (int i = 0; i < bitstring.length; ++i) {
            long b = bitstring[i];
            int count = 64 * i;
            while (b != 0) {
                if ((b & 1L) == 1L) {
                    setBits[n2] = count;
                    n2++;
                }
                b >>>= 1L;
                count++;
            }
        }
        
        assert(n == n2);
                
        return setBits;
    }

    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof VeryLongBitString)) {
            return false;
        }
        
        VeryLongBitString other = (VeryLongBitString)obj;
        
        if (nBits != other.nBits) {
            return false;
        }
        
        if (nSetBits != other.nSetBits) {
            return false;
        }
        
        return Arrays.equals(this.bitstring, other.bitstring);
    }    

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 11 * hash + Arrays.hashCode(this.bitstring);
        hash = 11 * hash + (int) (this.nBits ^ (this.nBits >>> 32));
        return hash;
    }
    
    /**
     * approximate the memory (in Bytes) used by this instance and its member variables, and 
     * return that portion on the stack.
     * @return 
     */
    protected long approximateMemoryUsed_Stack() {
       
        String arch = System.getProperty("sun.arch.data.model");
        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;
        int nbits = (is32Bit) ? 32 : 64;
        
        /*
        instance contains members:
            long[], long, long, long
        
        mem Stack:
           3 long primitives (a long is 2 times word size)
           1 array ref
        mem Heap:
           object overhead
           sum of each item in the long array
        */
        
        //TODO: should update this one day soon...
        
        long sumBits = 3*(nbits*2) + 32;
        
        long sumBytes = (sumBits/8);
        if (sumBytes*8 > sumBits) {
            sumBytes++;
        }
        
        return sumBytes;
    }
    
    /**
     * approximate the memory (in Bytes) used by this instance and its member variables, and 
     * return that portion on the heap.
     * @return 
     */
    protected long approximateMemoryUsed_Heap() {
        
        //TODO: update this using ObjectSpaceEstimator
        
        /*
        see comments within approximateMemoryUsed_Stack
        
        mem Heap:
           object overhead
           sum of each item in the long array  (a long is 2*word size)
        */
        
        String arch = System.getProperty("sun.arch.data.model");
        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;
        int nbits = (is32Bit) ? 32 : 64;
        
        long sumBytes = 16;
        if (bitstring != null) {
            sumBytes += bitstring.length * (2 * nbits);
        }
        long padding = (sumBytes % 8);
        sumBytes += padding;
        
        return sumBytes;
    }

    /**
     * perform a bitwise 'AND' on this bitstring and otherBS to find
     * the the intersection of bits in both bitstrings.
     * @param otherBS
     * @return 
     */
    public VeryLongBitString and(VeryLongBitString otherBS) {
                
        VeryLongBitString out;
        
        if (nBits == otherBS.nBits || (nBits > otherBS.nBits)) {
            
            out = copy();
            
            for (int i = 0; i < otherBS.bitstring.length; ++i) {
                long bs = otherBS.bitstring[i];
                out.bitstring[i] &= bs;
            }
            
        } else { //if (nBits < otherBS.nBits) {
            
            out = otherBS.copy();
            
            for (int i = 0; i < bitstring.length; ++i) {
                long bs = bitstring[i];
                out.bitstring[i] &= bs;
            }
        }
        
        out.recountNSetBits();
        
        return out;
    }

    /**
     * return the number of bits different betweeen the bitstring
     * otherBS and this and return that count.
     * Note that otherBS must have the same number of bits.
     * @param otherBS
     * @return 
     */
    public long nBitsDifferent(VeryLongBitString otherBS) {
        
        if (nBits != otherBS.nBits) {
            throw new IllegalArgumentException("otherBS must have same"
                + " number of bits");
        }        
        
        long nDiff = 0;
        
        for (int i = 0; i < otherBS.bitstring.length; ++i) {
            
            long bs = otherBS.bitstring[i];
            long tbs = bitstring[i];
            
            // xor gives the  different bits
            long diff = bs ^ tbs;
            nDiff += Long.bitCount(diff);
        }
        
        return nDiff;
    }
    
    /**
     * perform a bitwise 'or' on this bitstring and otherBS to make
     * a union operation.
     * @param otherBS
     * @return 
     */
    public VeryLongBitString or(VeryLongBitString otherBS) {
                
        VeryLongBitString out;
        
        if (nBits == otherBS.nBits || (nBits > otherBS.nBits)) {
            
            out = copy();
            
            for (int i = 0; i < otherBS.bitstring.length; ++i) {
                long bs = otherBS.bitstring[i];
                out.bitstring[i] |= bs;
            }
            
        } else { //if (nBits < otherBS.nBits) {
            
            out = otherBS.copy();
            
            for (int i = 0; i < bitstring.length; ++i) {
                long bs = bitstring[i];
                out.bitstring[i] |= bs;
            }
        }
        
        out.recountNSetBits();
        
        return out;
    }
    
    /**
     * perform a bitwise 'xor' on this bitstring and otherBS.
     * The result is the bits which are different.
     * @param otherBS
     * @return 
     */
    public VeryLongBitString xor(VeryLongBitString otherBS) {
                
        VeryLongBitString out;
        
        if (nBits == otherBS.nBits || (nBits > otherBS.nBits)) {
            
            out = copy();
            
            for (int i = 0; i < otherBS.bitstring.length; ++i) {
                long bs = otherBS.bitstring[i];
                out.bitstring[i] ^= bs;
            }
            
        } else { //if (nBits < otherBS.nBits) {
            
            out = otherBS.copy();
            
            for (int i = 0; i < bitstring.length; ++i) {
                long bs = bitstring[i];
                out.bitstring[i] ^= bs;
            }
        }
        
        out.recountNSetBits();
        
        return out;
    }

    /**
     * find the bits in this.copy() which are not in otherBS by performing
     * perform a bitwise 'AND' on this bitstring and otherBS to find
     * the intersection bits then clear those bits in the copy of this instance.
     * @param otherBS
     * @return 
     */
    public VeryLongBitString difference(VeryLongBitString otherBS) {
                
        VeryLongBitString out = copy();
        
        if (nBits == otherBS.nBits || (nBits > otherBS.nBits)) {
                        
            for (int i = 0; i < otherBS.bitstring.length; ++i) {
                long bs = otherBS.bitstring[i];
                long intersection = out.bitstring[i] & bs;
                // clear the intersection bits 
                out.bitstring[i] &= ~intersection;
            }
            
        } else { //if (nBits < otherBS.nBits) {
            
            out = otherBS.copy();
            
            for (int i = 0; i < bitstring.length; ++i) {
                long bs = bitstring[i];
                long intersection = out.bitstring[i] & bs;
                // clear the intersection bits 
                out.bitstring[i] &= ~intersection;
            }
        }
        
        out.recountNSetBits();
        
        return out;
    }

    /**
     * where bits2 are set in bits1, unset the bits in bits1.
     * This is the bitwise 'subtract' operation A & ~B.
     * @param bs1
     * @param bs2
     * @return 
     */
    public static VeryLongBitString subtract(VeryLongBitString bs1,
        VeryLongBitString bs2) {
        
        if (bs1.nBits != bs2.nBits) {
            throw new IllegalArgumentException("bs1 and bs2 must be same lengths");
        }
                
        VeryLongBitString out = bs1.copy();
                                
        for (int i = 0; i < bs2.bitstring.length; ++i) {
            out.bitstring[i] &= ~bs2.bitstring[i];
        }
            
        out.recountNSetBits();
        
        return out;
    }

    /**
     * 
     * @param nthBit
     * @return index for next highest bit set, else -1 if there is none.
     */
    public int nextHighestBitSet(int nthBit) {
        
        if (nthBit > nBits) {
            throw new IllegalArgumentException("bit must be < " + nBits);
        }
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        // nthBit is last item on idx array, increment idx
        if (nthBit == (((idx + 1) * ITEM_BIT_LENGTH) - 1)) {
            idx++;
        }
        
        for (int i = idx; i < bitstring.length; ++i) {
            
            long b = bitstring[i];
            
            if (i == idx) {                
                long shift = bitIdx + 1;
                b >>= shift;
                b <<= shift;
            }

            // same as v & -v, but might be using native method:
            long l = Long.lowestOneBit(b);
            
            if (l != 0) { 
                
                int bn = 63 - Long.numberOfLeadingZeros(l);
                
                // bn is bit position within row i
                
                return (int)((i * ITEM_BIT_LENGTH) + bn);                  
            }
        }
        
        return -1;
        
    }
    
    /**
     * returns the index of the lowest bit set ("rightmost), else -1 if no bits are set.
     * @return 
     */
    public int leastBitSet() {
        
        for (int i = 0; i < bitstring.length; ++i) {
            long b = bitstring[i];
                
            // same as v & -v, but might be using native method:
            long l = Long.lowestOneBit(b);
               
            if (l > 0) {
                
                int bn = 63 - Long.numberOfLeadingZeros(l);
                
                // bn is bit position within row i
                
                int bitIdx = (int)((i * ITEM_BIT_LENGTH) + bn);
                
                return bitIdx;
            }
        }
        
        return -1;
    }
    
    /**
     * returns the index of the highest bit set ("leftmost"), else -1 if no bits are set.
     * @return 
     */
    public int highestBitSet() {
        
        for (int i = (bitstring.length - 1); i > -1; --i) {
            long b = bitstring[i];
                
            long l = Long.highestOneBit(b);
               
            if (l > 0) {
             
                int bn = 63 - Long.numberOfLeadingZeros(l);
                
                // bn is bit position within row i
                
                int bitIdx = (int)((i * ITEM_BIT_LENGTH) + bn);
                
                return bitIdx;
            }
        }
        
        return -1;
    }
}
