package algorithms;

import java.util.Arrays;
import java.util.concurrent.Semaphore;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * permute the given set of numbers in a thread that waits for the getNext()
 * invocation to calculate the next permutation.
   <pre>
   The permute code is adapted from 
        from https://en.wikipedia.org/wiki/Heap%27s_algorithm
   which further references:
        Sedgewick, Robert. "a talk on Permutation Generation Algorithms
        http://www.cs.princeton.edu/~rs/talks/perms.pdf
   
   The semaphore model is adapted from Chap 12.1 of "Java Concurrency in Practice"
   by Goetz et al.
   </pre>
 */
public class PermutationsWithAwait {
    
    private final Semaphore availableItem, computationLock;
    
    /**
     * the current permutation
     */
    private final int[] x;
    
    private AtomicBoolean finished;
    
    //private final BigInteger nPermutations;
    
    /**
     * permute the given set of numbers in a thread that waits for the getNext()
     * invocation to calculate the next permutation.

       <pre>
       The permute code is adapted from 
            from https://en.wikipedia.org/wiki/Heap%27s_algorithm
       which further references:
            Sedgewick, Robert. "a talk on Permutation Generation Algorithms
            http://www.cs.princeton.edu/~rs/talks/perms.pdf

       The semaphore model is adapted from Chap 12.1 of "Java Concurrency in Practice"
       by Goetz et al.
       </pre>
       
     * @param seq
     * @throws java.lang.InterruptedException
     */
    public PermutationsWithAwait(int[] seq) throws InterruptedException {
                
        int n = seq.length;
        this.x = new int[n];
        
        this.availableItem = new Semaphore(0);
        this.computationLock = new Semaphore(1);
        
        //nPermutations = MiscMath0.factorialBigInteger(n);
        
        computationLock.acquire();
        
        //output(A)
        System.arraycopy(seq, 0, x, 0, n);
        //nCurrent = BigInteger.ONE;
        
        availableItem.release();
        
        finished = new AtomicBoolean(false);
                
        Thread thread = new Thread(new Permuter(seq, finished));
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
        
        if (finished.get()) {
            return;
        }
                
        availableItem.acquire();
        
        System.arraycopy(x, 0, out, 0, out.length);
        
        computationLock.release();        
    }
    
    //TODO: consider using Callable so run can throw an exception
    private class Permuter implements Runnable {
        private final int[] s;
        final AtomicBoolean permDone;
        Permuter(int[] seq, AtomicBoolean permDone) {
           this.s = Arrays.copyOf(seq, seq.length);
           this.permDone = permDone;
        }

        @Override
        public void run() {
            int n = s.length;
            int[] c = new int[n];
            int i = 0;
            int swap;
            while (i < n) {
                if (c[i] < i) {
                    if ((i & 1) != 1) {
                        // i is even number
                        swap = s[0];
                        s[0] = s[i];
                        s[i] = swap;
                    } else {
                        swap = s[c[i]];
                        s[c[i]] = s[i];
                        s[i] = swap;
                    }
                                        
                    try {
                        computationLock.acquire();
                    } catch (InterruptedException ex) {
                        Thread.currentThread().interrupt();
                    }
                    
                    // output permutation to instance member x
                    System.arraycopy(s, 0, x, 0, n);
                    
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
            permDone.set(true);
            System.out.println("finished permutations, exiting runloop");
        }
    }
}
