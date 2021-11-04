package algorithms;

import java.util.Arrays;
import java.util.concurrent.Semaphore;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * permute the given set of numbers in a thread that waits for the getNext()
 * invocation to calculate the next permutation.
     
   The permute code is adapted from 
        from https://en.wikipedia.org/wiki/Heap%27s_algorithm
   who reference:
        Sedgewick, Robert. "a talk on Permutation Generation Algorithms
        http://www.cs.princeton.edu/~rs/talks/perms.pdf
   
   The semaphore model is adapted from Chap 12.1 of "Java Concurrency in Practice"
   by Goetz et al.
 */
public class PermutationsWithAwait {
    
    private Logger log = null;
    private Level LEVEL = Level.FINE;
    
    private final Semaphore availableItem, computationLock;
    
    /**
     * the current permutation
     */
    private final int[] x;
    
    //private final BigInteger nPermutations;
    
    /*
     * TODO: review volatile w.r.t. BigInteger
     * 
     */
    //private volatile BigInteger nCurrent;
       
    /**
     * permute the given set of numbers in a thread that waits for the getNext()
     * invocation to calculate the next permutation.

       The permute code is adapted from 
            from https://en.wikipedia.org/wiki/Heap%27s_algorithm
       who reference:
            Sedgewick, Robert. "a talk on Permutation Generation Algorithms
            http://www.cs.princeton.edu/~rs/talks/perms.pdf

       The semaphore model is adapted from Chap 12.1 of "Java Concurrency in Practice"
       by Goetz et al.
     * @param set
     * @throws java.lang.InterruptedException
     */
    public PermutationsWithAwait(int[] set) throws InterruptedException {
        
        log = Logger.getLogger(this.getClass().getSimpleName());
        
        int n = set.length;
        this.x = new int[n];
        
        this.availableItem = new Semaphore(0);
        this.computationLock = new Semaphore(1);
        
        //nPermutations = MiscMath0.factorialBigInteger(n);
        
        log.log(LEVEL, "BEFORE computationLock.acquire()");
        computationLock.acquire();
        log.log(LEVEL, "   AFTER computationLock.acquire()");
        
        //output(A)
        System.arraycopy(set, 0, x, 0, n);
        //nCurrent = BigInteger.ONE;
        
        //log.log(LEVEL, String.format("*x=%s  nCurr=%s out of nPerm=%s\n", 
        //   Arrays.toString(x), nCurrent.toString(), nPermutations.toString()));
                 
        log.log(LEVEL, "BEFORE availableItem.release()");
        availableItem.release();
        log.log(LEVEL, "   AFTER availableItem.release()");
        
        log.log(LEVEL, "starting Permuter");
        
        Thread thread = new Thread(new Permuter(set));
        thread.start();
    }
    
    /**
     * get the next permutation of the original set
     * @param out an output array that will hold the results.  note that out.length
     * must be the same length as the set used during construction of this object instance.
     * @throws InterruptedException 
     */
    public void getNext(int[] out) throws InterruptedException {
        if (out.length != x.length) {
            throw new IllegalArgumentException("out.length must equal original set.length given to constructor");
        }
        
        //log.log(LEVEL, "getNext BEFORE availableItem.acquire() nc={0}", 
        //    nCurrent.toString());
        
        availableItem.acquire();
        
        System.arraycopy(x, 0, out, 0, out.length);
        
        computationLock.release();
        
        //log.log(LEVEL, "  getNext AFTER computationLock.release() nc={0}", nCurrent.toString());
    }
    
    //TODO: consider using Callable so run can throw an exception
    private class Permuter implements Runnable {
        private final int[] set;
        Permuter(int[] set) {
           this.set = Arrays.copyOf(set, set.length);
        }

        @Override
        public void run() {
            log.log(LEVEL, "Permuter.run()");
            int n = set.length;
            int[] c = new int[n];
            int i = 0;
            int swap;
            while (i < n) {
                if (c[i] < i) {
                    if ((i & 1) != 1) {
                        // i is even number
                        swap = set[0];
                        set[0] = set[i];
                        set[i] = swap;
                    } else {
                        swap = set[c[i]];
                        set[c[i]] = set[i];
                        set[i] = swap;
                    }
                    
                    log.log(LEVEL, "BEFORE computationLock.acquire()");
                    
                    try {
                        computationLock.acquire();
                    } catch (InterruptedException ex) {
                        log.log(Level.SEVERE, "unrecoverable error? " + ex.getMessage());
                        //TODO: change to using Callable so can throw an exception
                        return;
                    }
                    
                    // output permutation to instance member x
                    System.arraycopy(set, 0, x, 0, n);
                    
                    //nCurrent = nCurrent.add(BigInteger.ONE);                    
                    
                    //log.log(LEVEL, String.format("*x=%s  nCurr=%s out of nPerm=%s\n", 
                    //    Arrays.toString(x), nCurrent.toString(), nPermutations.toString()));
                   
                    availableItem.release();
                    
                    //log.log(LEVEL, "   AFTER cavailableItem.release()");
                    
                    //Swap has occurred ending the for-loop. Simulate the increment 
                    //of the for-loop counter
                    c[i] += 1;
                    //Simulate recursive call reaching the base case by bringing the 
                    //pointer to the base case analog in the array
                    i = 0;
                } else {
                    //Calling generate(i+1, A) has ended as the for-loop terminated. 
                    //Reset the state and simulate popping the stack by incrementing the pointer.
                    c[i] = 0;
                    i++;
                }
            }
        }
    }
}
