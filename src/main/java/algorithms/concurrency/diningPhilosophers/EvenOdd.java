package algorithms.concurrency.diningPhilosophers;

import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.ReentrantLock;

/**
 * a toy model close to sequential for timing comparisons.
 * All odd numbered diners eat then all even numbered diners eat.
 */
public class EvenOdd extends AbstractDiningPhilosophers {
    ReentrantLock eatLock = new ReentrantLock();

    // needed for the alternating first last diner when n is an odd number
    int nEvenIter = 0;

    /**
     * constructor for dining philospher's abstract
     * @param n number of diners
     * @param thinkTimeMsec the time to think between trying to eat
     * @param label a name of the implementation useful in debugging
     */
    public EvenOdd(int n, int thinkTimeMsec) {
        super(n, thinkTimeMsec, EvenOdd.class.getSimpleName());
    }

    public void dine(int bitesPer, int timeoutSec) throws InterruptedException {

        initDine(bitesPer);

        Thread thrE = new Thread(groupEat(0, timeoutSec));
        Thread thrO = new Thread(groupEat(1, timeoutSec));
        thrE.start();
        thrO.start();
        thrE.join();
        thrO.join();

        closeDine();
    }//end dine

    protected boolean foodOnPlates(int evenOdd) {
        for (int i = 0; i < n; ++i) {
            if ((i & 1) == evenOdd && bites.containsKey(i)) {
                return true;
            }
        }
        return false;
    }

    protected Runnable groupEat(int evenOdd, int timeOut) {
        return new Runnable() {
            public void run() {
                while (foodOnPlates(evenOdd)) {
                    try {
                        if (eatLock.tryLock(timeOut, TimeUnit.SECONDS)) {
                            for (int i = 0; i < n; ++i) {
                                if ((i & 1) != evenOdd || !bites.containsKey(i)) {
                                    continue;
                                }

                                // if nIsOdd and this is the even eating group, we need to alternate
                                // the eating between first and last diner.
                                //  when nEvenIter is even, 0 eats, else n-1 eats
                                if (nIsOdd && evenOdd == 0) {
                                    if (i == (n - 1) && ((nEvenIter & 1) == 0)) {
                                        continue;
                                    } else if (i == 0 && ((nEvenIter & 1) != 0)) {
                                        continue;
                                    }
                                }
                                if (bites.containsKey(i)) {
                                    if (bites.get(i) == 1) {
                                        bites.remove(i);
                                    } else {
                                        bites.put(i, bites.get(i) - 1);
                                    }
                                    //System.out.printf("%d eats\n", i);
                                }
                            }// end i loop
                        }// end locked eat
                    } catch (InterruptedException e) {
                        System.out.printf("eatLock interrupted\n");
                        Thread.currentThread().interrupt();
                    } finally {
                        if (eatLock.isHeldByCurrentThread()) {
                            eatLock.unlock();
                            if (evenOdd == 0) ++nEvenIter;
                            try {
                                Thread.sleep(thinkTimeMilliSec);
                            } catch (InterruptedException e) {
                                Thread.currentThread().interrupt();
                            }
                        }
                    }
                }
            }// end while food on plate
        };// end run

    }// end groupEat
}
