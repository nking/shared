package algorithms;

import java.util.Arrays;

/**
 * a stack with a single large primitive array internally
 * that can expand dynamically, but will not contract.
 * 
 * This is intended for use when there are a large number of numbers to process
 * and the linked list of objects in a standard stack consumes too much memory.
 * (a java object is 16B overhead + content, so N*16Bytes is a rough approx
 * and can be compared to the memory here which will be
 * capacity * 4Bytes).
 * 
 * @author nichole
 */
public class StackIntLarge {
    
    // circular array
    protected int[] a = null;
    
    protected int idxLast = -1;
    
    public StackIntLarge(int capacity) {
        
        if (capacity < 1) {
            throw new IllegalArgumentException("capacity must be > 0");
        }
        
        a = new int[capacity];
    }
    
    public void push(int value) {
        
        expandIfNeeded();
        
        if (idxLast == -1) {
            idxLast = 0;
            a[idxLast] = value;
            return;
        }
        
        int idxNext = idxLast + 1;
        assert(idxNext < a.length);
        idxLast = idxNext;
        a[idxLast] = value;
    }
    
    public int pop() {
        
        if (idxLast == -1) {
            throw new IllegalStateException("stack is empty");
        }
        
        if (idxLast == 0) {
            idxLast = -1;
            return a[0];
        }

        idxLast--;
        return a[idxLast + 1];        
    }
    
    public int peek() {
        
        if (idxLast == -1) {
            throw new IllegalStateException("stack is empty");
        }
        
        return a[idxLast];
    }
    
    public boolean isEmpty() {
        return (size() == 0);
    }
    
    public int size() {
     
        if (idxLast == -1) {
            return 0;
        }
        
        return idxLast + 1;
    }
        
    private void expandIfNeeded() {
        
        int sz = size();
        
        if (sz != a.length) {
            return;
        }
        
        assert(idxLast > -1);
        
        // expand by 10 percent?
        int nAdd = (int)Math.ceil(sz * 0.1f);
        if (nAdd < 16) {
            nAdd = 16;
        }
                 
        a = Arrays.copyOf(a, sz + nAdd);
    }
    
}
