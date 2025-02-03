package algorithms.concurrency.diningPhilosophers;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.ReentrantLock;

/**
 * class implementing the hierarchical solution.
 * each diner attempts to eat by choosing the lowest order fork
 * from their 2 choices, and if cannot obtain it, puts it down
 * and tries again.
 * if obtains the lower numbered fork, tries to obtain the higher
 * number and if fails, puts both forks down and tries again,
 * else if succeeds, eats, then puts both forks down.
 *
 * there is no interprocess communication between diner's (they're busy
 * thinking is the metaphor).
 */
public class Hierarchical extends AbstractDiningPhilosophers {
    protected ReentrantLock[] forkLocks;

    /**
     * constructor
     * @param n number of diners
     * @param thinkTimeMsec the time to think between trying to eat
     */
    public Hierarchical(int n, int thinkTimeMsec) {
        super(n, thinkTimeMsec, Hierarchical.class.getSimpleName());

        this.forkLocks = new ReentrantLock[n];
        for (int i = 0; i < n; ++i) {
            forkLocks[i] = new ReentrantLock();
        }
    }

    /**
     * constructor
     * @param n number of diners
     * @param thinkTimeMsec the time to think between trying to eat
     * @param label a name of the implementation useful in debugging
     */
    public Hierarchical(int n, int thinkTimeMsec, String label) {
        super(n, thinkTimeMsec, label);

        this.forkLocks = new ReentrantLock[n];
        for (int i = 0; i < n; ++i) {
            forkLocks[i] = new ReentrantLock(true);
        }
    }

    //   0   1   2
    // 2  0    1   2

    /**
     * get lower numbered fork for diner i
     * @param i number of diner
     * @return lower numbered fork
     */
    protected ReentrantLock getLowerFork(int i) {
        if (i == 0) {
            return forkLocks[0];
        }
        return forkLocks[i - 1];
    }

    /**
     * get higher numbered fork for diner 'i'
     * @param i diner number
     * @return the higher numbered fork of diner 'i'
     */
    protected ReentrantLock getHigherFork(int i) {
        if (i == 0) {
            return forkLocks[n - 1];
        }
        return forkLocks[i];
    }

    /**
     * begin dining of all diners
     * @param bitesPer the number of bites to put on each diner's plate
     * @param timeoutSec the timeout for waiting to start eating
     * @throws InterruptedException thrown if a thread is interrupted
     */
    public void dine(int bitesPer, int timeoutSec) throws InterruptedException {

        initDine(bitesPer);

        // start all threads.
        // each tries to obtain low number fork, high number fork, eat, release forks. think
        List<Thread> threadList = new ArrayList<>();
        for (int i = 0; i < n; ++i) {
            Thread thr = new Thread(eat(i, timeoutSec));
            threadList.add(thr);
            thr.start();
        }
        for (int i = 0; i < n; ++i) {
            threadList.get(i).join();
        }

        closeDine();
    }//end dine

    /**
     * a runnable implementing eat for diner 'i'
     * @param i
     * @param timeoutSec
     * @return
     */
    protected Runnable eat(int i, int timeoutSec) {
        return new Runnable() {
            @Override
            public void run() {
                while (bites.containsKey(i)) {
                    try {
                        boolean s = getLowerFork(i).tryLock(timeoutSec, TimeUnit.MILLISECONDS);
                        if (!s) continue;
                        s = getHigherFork(i).tryLock(timeoutSec, TimeUnit.MILLISECONDS);
                        if (!s) continue;
                        if (bites.containsKey(i)) {
                            if (bites.get(i) == 1) {
                                bites.remove(i);
                            } else {
                                bites.put(i, bites.get(i) - 1);
                            }
                            //System.out.printf("%d eats\n", i);
                        }
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                    } finally {
                        if (getLowerFork(i).isHeldByCurrentThread()) {
                            getLowerFork(i).unlock();
                        }
                        if (getHigherFork(i).isHeldByCurrentThread()) {
                            getHigherFork(i).unlock();
                            try {// think after eating
                                Thread.sleep(thinkTimeMilliSec);
                            } catch (InterruptedException e) {
                                Thread.currentThread().interrupt();
                            }
                        }
                    }
                }
            }
        };
    }
}
