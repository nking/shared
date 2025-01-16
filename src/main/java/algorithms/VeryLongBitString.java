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
    
    /**
     *
     */
    protected final long nBits;
    
    /**
     * length for signed long positive
     */
    protected static final long ITEM_BIT_LENGTH = 63;
    
    // the length of an array in java is limited to Integer.MAX_VALUE

    /**
     *
     */
    protected static final long N_MAX_BITS = Integer.MAX_VALUE * ITEM_BIT_LENGTH;
    
    /**
     *
     */
    protected final long capacityBits;
    
    private long nSetBits = 0;
   
    /**
     *
     @param nBits the fixed capacity for this bitstring 
     */
    public VeryLongBitString(long nBits) {
        
        if (nBits > N_MAX_BITS) {
            throw new IllegalArgumentException("cannot hold more than " + N_MAX_BITS + " bits");
        }
        if (nBits < 0) {
            throw new IllegalArgumentException("nBits must be non-negative");
        }
        
        int nElements = getRowNumber(nBits) + 1;
        
        bitstring = new long[nElements];
        
        capacityBits = nElements * ITEM_BIT_LENGTH;
        
        this.nBits = nBits;
    }
    
    /**
     *
     @param bitstrings the bitstring to initialize this with
     @param nBits the fixed capacity for this bitstring 
     @param nSetBits the number of set bits in bitstrings.  a convenience to avoid counting the set bits.
     */
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
    
    /**
     *
     @return
     */
    public long getCapacity() {
        return capacityBits;
    }

    public void setAllBits() {
        for (int i = 0; i < bitstring.length - 1; ++i) {
            bitstring[i] = (1L << 63L) - 1;
        }
        // handle last row
        long bitIdx = getBitIdx(nBits, bitstring.length - 1);
        bitstring[bitstring.length - 1] = (1L << bitIdx) - 1;
        nSetBits = nBits;
    }
    /**
     *
     @param nthBit
     */
    public void setBit(long nthBit) {
        if (nthBit < 0) {
            throw new IllegalArgumentException("nthBit must be non-negative");
        }

        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        assert(bitIdx < ITEM_BIT_LENGTH);
        assert(bitIdx >= 0);
        
        // test bit is not set
        if ((bitstring[idx] & (1L << bitIdx)) == 0) {
          
            // set bit
            bitstring[idx] |= (1L << bitIdx);
                        
            nSetBits++;
        }
    }
    
    /**
     *
     @param nthBit
     */
    public void clearBit(long nthBit) {

        if (nthBit < 0) {
            throw new IllegalArgumentException("nthBit must be non-negative");
        }

        int idx = getRowNumber(nthBit);
        assert(idx >= 0 && idx < bitstring.length);
        
        int bitIdx = getBitIdx(nthBit, idx);
        assert(bitIdx >= 0);
        
        // test bit
        if ((bitstring[idx] & (1L << bitIdx)) != 0) {
        
            // clear bit
            bitstring[idx] &= ~(1L << bitIdx);
                        
            nSetBits--;
        }
    }
    
    /**
     *
     @param nthBit
     */
    public void toggleBit(long nthBit) {

        int idx = getRowNumber(nthBit);
        assert(idx >= 0 && idx < bitstring.length);

        int bitIdx = getBitIdx(nthBit, idx);
        assert(bitIdx >= 0);

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
    
    /**
     *
     @return
     */
    public long getNSetBits() {
        return nSetBits;
    }
    
    /**
     *
     @return
     */
    public long getInstantiatedBitSize() {
        return nBits;
    }
    
    /**
     *
     @param nthBit
     @return
     */
    public boolean isSet(long nthBit) {

        if (nthBit < 0) {
            throw new IllegalArgumentException("nthBit must be non-negative");
        }

        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        return ((bitstring[idx] & (1L << bitIdx)) != 0);
    }
    
    /**
     *
     @param nthBit
     @return
     */
    public boolean isNotSet(long nthBit) {
        if (nthBit < 0) {
            throw new IllegalArgumentException("nthBit must be non-negative");
        }

        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        return ((bitstring[idx] & (1L << bitIdx)) == 0);
    }
         
    /**
     *
     @param n
     @return
     */
    protected int getRowNumber(long n) {

        if (n < 0) {
            throw new IllegalArgumentException("n must be non-negative");
        }
        
        int nthElement = (int)(n/ITEM_BIT_LENGTH);
        
        return nthElement;
    }
    
    /**
     * given the arrayIdx, calculate the bit position within that item
     * for the very large bitstring position nthBit.  Note that if the
     * result is out of bounds, -1 is returned and the invoker should
     * handle that.
     * 
     @param nthBit bit  index
     @param arrayIdx array index
     @return  item stored in nthBit bit index of item stored at arrayIdx of internal array
     */
    int getBitIdx(long nthBit, int arrayIdx) {
        
        int bitIdx = (int)(nthBit - (arrayIdx * ITEM_BIT_LENGTH));
        
        if ((bitIdx < 0) || (bitIdx > (ITEM_BIT_LENGTH - 1))) {
            return -1;
        }
        
        return bitIdx;
    }
    
    /**
     *
     */
    public void clearAllBits() {
        
        Arrays.fill(bitstring, 0);
        
        nSetBits = 0;
    }
    
    /**
     *
     @return
     */
    public VeryLongBitString copy() {
        
        VeryLongBitString c = new VeryLongBitString(bitstring, 
            nBits, nSetBits);
        
        return c;
    }
    
    /**
     *
     @param other
     */
    public void resetAllTo(VeryLongBitString other) {
        if (other.nBits != nBits) {
            throw new IllegalArgumentException("nBits must be the same in both to use this method");
        }
        System.arraycopy(other.bitstring, 0, bitstring, 0, other.bitstring.length);
        
        this.nSetBits = other.nSetBits;
    }
    
    /**
     * get the '01...' bitstring representation of this object
     @return binary bitstring representation of this object
     */
    @Override
    public String toString() {
        /*
        bit 0 is at bitstring[0] bit 0
        bit n is at bitstring[bitstring.length-1] bit n % ITEM_LENGTH
        so reverse array order
         */
        StringBuilder sb = new StringBuilder();
        for (int i = (bitstring.length - 1); i > -1 ; i--) {
            sb.append(Long.toBinaryString(bitstring[i]));
        }
        return sb.toString();
    }

    /**
     * appends below the bitstring, 2 arrays to count the bits
     * @return
     */
    public String toStringWithRuler() {
        /*
        bit 0 is at bitstring[0] bit 0
        bit n is at bitstring[bitstring.length-1] bit n % ITEM_LENGTH
        so reverse array order
         */
        StringBuilder sbr1 = new StringBuilder();
        StringBuilder sbr2 = new StringBuilder();
        StringBuilder sb = new StringBuilder();
        long c = nBits - 1;
        for (int i = (bitstring.length - 1); i > -1 ; i--) {
            sb.append(Long.toBinaryString(bitstring[i]));
            for (int j = 0; j < ITEM_BIT_LENGTH; ++j,--c) {
                if (c < 0) break;
                sbr1.append(c % 10);
                if (c%10 == 0) {
                    sbr2.append(c/10);
                } else {
                    sbr2.append(" ");
                }
            }
        }
        sb.append("\n").append(sbr1).append("\n").append(sbr2);
        return sb.toString();
    }
    
    /**
     *
     */
    protected void recountNSetBits() {
        nSetBits = 0;
        for (int i = 0; i < bitstring.length; ++i) {
            nSetBits += Long.bitCount(bitstring[i]);
        }
    }
    
    /**
     * get a list of the bit numbers that are set.
     @return bit number of set bits
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
            int count = (int)ITEM_BIT_LENGTH * i;
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
     @return memory used in Bytes
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
     @return memory used in Bytes
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
     @param otherBS other bitstring
     @return the result of bitwise and of this bitstring and otherBS bitstring
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
     * return the number of bits different between the bit-string
     * otherBS and this and return that count.
     * Note that otherBS must have the same number of bits.
     @param otherBS other bitstring
     @return the nuber of bits that are different in bitwise comparison of otherBS with this bitstring
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
     @param otherBS other bitstring
     @return bitwise or of otherBS with this bitstring
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
     @param otherBS other bitstring
     @return bitwise xor of otherBS with this bitstring
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
     @param otherBS other bitstring
     @return bitwise difference between otherBS and this bitstring
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
     * This is the bitwise 'subtract' operation A 'bitwise and' ~B.
     * It subtracts from A any set bits in B.
     @param bs1 a bitstring
     @param bs2 another bitstring
     @return the difference between otherBS and this bitstring
     */
    public static VeryLongBitString subtract(VeryLongBitString bs1,
        VeryLongBitString bs2) {
        
        if (bs1.nBits != bs2.nBits) {
            throw new IllegalArgumentException("bs1 and bs2 must be same lengths");
        }
                
        VeryLongBitString out = bs1.copy();
                                
        for (int i = 0; i < bs2.bitstring.length; ++i) {
            // set subtraction is A & ~B
            out.bitstring[i] &= ~bs2.bitstring[i];
        }
            
        out.recountNSetBits();
        
        return out;
    }

    /**
     * 
     @param nthBit the bit index
     @return index for next highest bit set, else -1 if there is none.
     */
    public int nextHighestBitSet(int nthBit) {
        
        if (nthBit > nBits) {
            throw new IllegalArgumentException("bit must be < " + nBits);
        }
        if (nthBit < 0) {
            throw new IllegalArgumentException("bit must be >= 0");
        }
        
        int idx = getRowNumber(nthBit);
        
        int bitIdx = getBitIdx(nthBit, idx);
        
        // nthBit is last item on idx array, so start looking at next row, bit 0
        if (bitIdx == ((int)ITEM_BIT_LENGTH - 1)) {
            idx++;
            bitIdx = 0;
        } else {
            ++bitIdx;
        }

        if (bitIdx > 0) {
            // check the same array element for the next higher set bit
            long b = bitstring[idx];
            // mask out bits below bitIdx
            long mask = Long.MAX_VALUE >> (bitIdx - 0);
            mask <<= (bitIdx - 0);
            b &= mask;

            if (b == 0) {
                idx++;
                bitIdx = 0;
            } else {
                long lsbNumber = b & -b;
                int lsb = (int)(Math.log(lsbNumber) / Math.log(2));
                // transform back to bit number
                int nextBit = (idx * (int)ITEM_BIT_LENGTH) + lsb;
                return nextBit;
            }
        }

        for (int i = idx; i < bitstring.length; ++i) {
            
            long b = bitstring[i];

            if (b == 0) continue;

            long lsbNumber = b & -b;
            int lsb = (int)(Math.log(lsbNumber) / Math.log(2));
            // transform back to bit number
            int nextBit = (idx * (int)ITEM_BIT_LENGTH) + lsb;
            return nextBit;
        }
        
        return -1;
    }
    
    /**
     * returns the index of the lowest bit set ("rightmost", LSB), else -1 if no bits are set.
     @return the index of the lowest (rightmost) bit set, else -1 if no bits are set.
     */
    public int leastBitSet() {
        
        for (int i = 0; i < bitstring.length; ++i) {
            long b = bitstring[i];
            if (b == 0) continue;
                
            long l = b & -b;
            int lsb = (int)(Math.log(l) / Math.log(2));
            // transform back to bit number
            int nextBit = (i * (int)ITEM_BIT_LENGTH) + lsb;
            return nextBit;
        }
        
        return -1;
    }
    
    /**
     * returns the index of the highest bit set ("leftmost", MSB), else -1 if no bits are set.
     @return the index of the highest bit set (leftmost) else -1 if no bits are set.
     */
    public int highestBitSet() {
        
        for (int i = (bitstring.length - 1); i > -1; --i) {
            long b = bitstring[i];
            if (b == 0) continue;
                
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
