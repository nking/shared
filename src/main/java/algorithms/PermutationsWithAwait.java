package algorithms;

import java.util.Arrays;
import java.util.concurrent.Semaphore;

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
    
    private final Semaphore availableItem, computationLock;
    
    /**
     * the current permutation
     */
    private final int[] x;
    
    //private final BigInteger nPermutations;
    
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
                
        int n = set.length;
        this.x = new int[n];
        
        this.availableItem = new Semaphore(0);
        this.computationLock = new Semaphore(1);
        
        //nPermutations = MiscMath0.factorialBigInteger(n);
        
        computationLock.acquire();
        
        //output(A)
        System.arraycopy(set, 0, x, 0, n);
        //nCurrent = BigInteger.ONE;
        
        availableItem.release();
                
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
        
        availableItem.acquire();
        
        System.arraycopy(x, 0, out, 0, out.length);
        
        computationLock.release();        
    }
    
    //TODO: consider using Callable so run can throw an exception
    private class Permuter implements Runnable {
        private final int[] set;
        Permuter(int[] set) {
           this.set = Arrays.copyOf(set, set.length);
        }

        @Override
        public void run() {
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
                                        
                    try {
                        computationLock.acquire();
                    } catch (InterruptedException ex) {
                        //TODO: change to using Callable so can throw an exception
                        //or use an exit() here
                        return;
                    }
                    
                    // output permutation to instance member x
                    System.arraycopy(set, 0, x, 0, n);
                    
                    //nCurrent = nCurrent.add(BigInteger.ONE);                    
                    
                    availableItem.release();
                                        
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
