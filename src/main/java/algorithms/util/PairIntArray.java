package algorithms.util;

import algorithms.Rotate;
import algorithms.sort.MiscSorter;
import java.util.Arrays;
import java.util.Set;

/**
 * class to hold x and y arrays of points
 * 
   first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 * @author nichole
 */
public class PairIntArray {
    
    /**
     *
     */
    protected int[] x = null;
    
    /**
     *
     */
    protected int[] y = null;
    
    /**
     *
     */
    protected int n = 0;
    
    /**
     *
     @param capacity
     */
    public PairIntArray(int capacity) {
        
        x = new int[capacity];
        
        y = new int[capacity];
    }
    
    /**
     *
     */
    public PairIntArray() {
        
        x = new int[10];
        
        y = new int[10];
    }
    
    /**
     *
     @return
     */
    public int getN() {
        return n;
    }
    
    /**
     *
     @param xPoint
     @param yPoint
     */
    public void add(int xPoint, int yPoint) {
        
        expandIfNeeded(n + 1);
        
        x[n] = xPoint;
        y[n] = yPoint;
        
        n++;
    }
    
    /**
     *
     @param xy
     */
    public void add(PairInt xy) {
        
        expandIfNeeded(n + 1);
        
        x[n] = xy.getX();
        y[n] = xy.getY();
        
        n++;
    }
    
    /**
     *
     @param other
     */
    public void addAll(PairIntArray other) {
        
        expandIfNeeded(n + other.getN());        
        
        System.arraycopy(other.getX(), 0, x, n, other.getN());
        System.arraycopy(other.getY(), 0, y, n, other.getN());
            
        n += other.getN();
    }
    
    /**
     *
     @param other
     */
    public void addAll(Set<PairInt> other) {
        
        expandIfNeeded(n + other.size());        
        
        for (PairInt p : other) {
            add(p);
        }
    }
    
    /**
     * remove indexes from idxLo to idxHi, inclusive
     @param idxLo
     @param idxHi 
     */
    public void removeRange(int idxLo, int idxHi) {
        
        if ((idxLo < 0) || (idxLo > (n - 1))) {
            throw new IllegalArgumentException(idxLo + ", idxLo is out of range");
        }
        if ((idxHi < 0) || (idxHi > (n - 1))) {
            throw new IllegalArgumentException(idxHi + ",idxHi is out of range");
        }
        if (idxHi < idxLo) {
            throw new IllegalArgumentException("idxHi has to be >= idxLo");
        }
        
        int nRemove = idxHi - idxLo + 1;

        int moveIdx = idxHi + 1;
        if (moveIdx <= (n - 1)) {
            for (int moveToIdx = idxLo; moveToIdx < (n - nRemove); moveToIdx++) {
                x[moveToIdx] = x[moveIdx];
                y[moveToIdx] = y[moveIdx];
                moveIdx++;
            }
        }
        
        // not necessary, but easier debugging to reset the last nRemove to 0
        for (int i = (n - nRemove); i < n; i++) {
            x[i] = 0;
            y[i] = 0;
        }
        
        n -= nRemove;
    }
    
    /**
     *
     @param numberOfInserts
     */
    public void insertSpaceAtTopOfArrays(int numberOfInserts) {
        
        if (x.length >= (n + numberOfInserts)) {
            
            for (int i = (n - 1); i > -1; i--) {
                x[i + numberOfInserts] = x[i];
                y[i + numberOfInserts] = y[i];
            }
            for (int i = 0; i < numberOfInserts; i++) {
                x[i] = 0;
                y[i] = 0;
            }
            
        } else {
            int[] xx = new int[n + numberOfInserts];
            int[] yy = new int[n + numberOfInserts];
            System.arraycopy(x, 0, xx, numberOfInserts, n);
            System.arraycopy(y, 0, yy, numberOfInserts, n);
            x = xx;
            y = yy;
        }
        
        n += numberOfInserts;
    }
    
    /**
     *
     @param index
     @param xPoint
     @param yPoint
     */
    public void set(int index, int xPoint, int yPoint) {
        
        if (index < 0) {
            throw new IllegalArgumentException("index is out of bounds of arrays");
        }
        
        expandIfNeeded(index + 1);
        
        x[index] = xPoint;
        y[index] = yPoint;
    }
    
    /**
     *
     @param other
     */
    public void swapContents(PairIntArray other) {
        
        int[] swap = x;
        x = other.x;
        other.x = swap;
        
        swap = y;
        y = other.y;
        other.y = swap;
        
        int swap2 = n;
        n = other.n;
        other.n = swap2;
    }
    
    /**
     *
     @param index
     @return
     */
    public int getX(int index) {
        if (index > (n - 1) || (index < 0)) {
            throw new IllegalArgumentException(
            "index is out of range in getX.  n=" + n + 
             " index=" + index);
        }
        return x[index];
    }
    
    /**
     *
     @param index
     @return
     */
    public int getY(int index) {
        if (index > (n - 1) || (index < 0)) {
            throw new IllegalArgumentException(
            "index is out of range in getY.  n=" + n + 
             " index=" + index);
        }
        return y[index];
    }
    
    /**
     *
     @return
     */
    public int[] getX() {
        return x;
    }
    
    /**
     *
     @return
     */
    public int[] getY() {
        return y;
    }
    
    /**
     *
     @param nTotal
     */
    protected void expandIfNeeded(int nTotal) {
        
        if (nTotal > x.length) {
            
            int n2 = x.length + 10;
            
            if (nTotal > n2) {
                n2 = nTotal;
            }
            
            x = Arrays.copyOf(x, n2);
            
            y = Arrays.copyOf(y, n2);
        }
    }
    
    /**
     *
     */
    public void reverse() {
        
        if (n < 2) {
            return;
        }
                
        int end = n >> 1;
        // 0 1 2 3 4
        for (int i = 0; i < end; i++) {
            int idx2 = n - i - 1;
            int swap = x[i];
            x[i] = x[idx2];
            x[idx2] = swap;
            
            swap = y[i];
            y[i] = y[idx2];
            y[idx2] = swap;
        }
    }
    
    /**
     *
     @param offset
     */
    public void rotateLeft(int offset) {
        Rotate r = new Rotate();
        r.rotate2(x, n, offset);
        r.rotate2(y, n, offset);
    }
    
    /**
     *
     @param index
     @param xPoint
     @param yPoint
     */
    public void insert(int index, int xPoint, int yPoint) {
        if (index < 0 || (index > n)) {
            throw new IllegalArgumentException("index is out of bounds of arrays");
        }
        
        expandIfNeeded(n + 1);
                
        // move everything at index thru n-1 to higher index
        for (int i = n; i > index; i--) {
            x[i] = x[i - 1];
            y[i] = y[i - 1];
        }
        
        x[index] = xPoint;
        y[index] = yPoint;
        
        n++;
    }

    /**
     *
     @param insertAtIndex
     @param insert
     */
    public void insertAll(int insertAtIndex, PairIntArray insert) {
        
        if (insertAtIndex < 0) {
            throw new IllegalArgumentException("insertAtIndex must be >= 0");
        } else if (insertAtIndex > n) {
            throw new IllegalArgumentException("insertAtIndex is out of bounds of array");
        }
        
        if (insertAtIndex == n) {
            addAll(insert);
            return;
        }
                
        int nTotal = n + insert.getN();
        
        expandIfNeeded(nTotal);
         
        /*
        copy everything at insertAtIndex           through  (n-1) 
                        to (nTotal - insert.n - 1)    to    (nTotal - 1)
        */
        
        int nMoveLen = (n - insertAtIndex);
        int dest0 = (nTotal - nMoveLen);
        
        System.arraycopy(x, insertAtIndex, x, dest0, nMoveLen);
        System.arraycopy(y, insertAtIndex, y, dest0, nMoveLen);
        
        System.arraycopy(insert.x, 0, x, insertAtIndex, insert.n);
        System.arraycopy(insert.y, 0, y, insertAtIndex, insert.n);
        
        n = nTotal;
        
        /*
        n=5.  insIdx=1,  insert.n=3 --> nTotal=5+3=8
        0 
           <---
        1
        2
        3
        4
        
        0
        1
        2
        3
        4 <-- 1
        5 <-- 2 (n-3)
        6 <-- 3 (n-2)
        7 <-- 4 (n-1)
        */
    }
    
    /**
     * reverse the indexes from 0 to lastSwapIdx, inclusive.
     * 
     @param lastSwapIdx 
     */
    public void reverse0toIdx(int lastSwapIdx) {
        
        if ((lastSwapIdx < 0) || (lastSwapIdx > (n - 1))) {
            throw new IllegalArgumentException("lastSwapIdx is out of bounds");
        }
        
        int nSep = (lastSwapIdx + 1) >> 1;
        
        for (int idx = 0; idx < nSep; ++idx) {
            int idx2 = lastSwapIdx - idx;
            int swapX = x[idx];
            int swapY = y[idx];
            x[idx] = x[idx2];
            y[idx] = y[idx2];
            x[idx2] = swapX;
            y[idx2] = swapY;
        }
    }
    
    /**
     * reverse the values between index firstSwapIdx and the last index in arrays.
     * 
     @param firstSwapIdx 
     */
    public void reverseIdxtoEnd(int firstSwapIdx) {
        
        if ((firstSwapIdx < 0) || (firstSwapIdx > (n - 1))) {
            throw new IllegalArgumentException("firstSwapIdx is out of bounds");
        }
        
        int count = 0;
        int nSep = (n - firstSwapIdx) >> 1;
        for (int idx = firstSwapIdx; idx < (firstSwapIdx + nSep); ++idx) {
            int idx2 = n - count - 1;
            int swapX = x[idx];
            int swapY = y[idx];
            x[idx] = x[idx2];
            y[idx] = y[idx2];
            x[idx2] = swapX;
            y[idx2] = swapY;
            count++;
        }
    }
    
    /**
     * sort by x ascending, then break ties
     * by y ascending.
     */
    public void sortByXAsc() {
        MiscSorter.sortBy1stArgThen2nd(x, y, 0, n - 1);
    }
    
    /**
     * copy inclusively the range startIdx through stopIdx.
     @param startIdx first index copied
     @param stopIdx last index copied
     @return 
     */
    public PairIntArray copyRange(int startIdx, int stopIdx) {
        
        if (startIdx < 0 || startIdx >= stopIdx) {
            throw new IllegalArgumentException("startIdx must be in range and "
                + " < stopIdx");
        }
        if (stopIdx < 0 || stopIdx > (n - 1)) {
            throw new IllegalArgumentException("stopIdx must be in range");
        }
        
        int len = stopIdx - startIdx + 1;
        
        PairIntArray sub = new PairIntArray(len);
        
        sub.insertSpaceAtTopOfArrays(len);
        System.arraycopy(x, startIdx, sub.x, 0, len);
        System.arraycopy(y, startIdx, sub.y, 0, len);
        sub.n = len;
        
        return sub;
    }
    
    /**
     *
     @return
     */
    public PairIntArray copy() {
        
        PairIntArray clone = new PairIntArray(n);
        
        clone.insertSpaceAtTopOfArrays(n);
        System.arraycopy(x, 0, clone.x, 0, n);
        System.arraycopy(y, 0, clone.y, 0, n);
        clone.n = n;
        
        return clone;
    }
    
    /**
     *
     @return
     */
    public PairFloatArray toPairFloatArray() {
        PairFloatArray out = new PairFloatArray();
        for (int i = 0; i < n; i++) {
            out.add(x[i], y[i]);
        }
        return out;
    }
    
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append("[").append(i).append("]")
            .append("(").append(x[i]).append(", ").append(y[i]).append(") ");
        }
        return sb.toString();
    }
 
}
