package algorithms.concurrency.diningPhilosophers;

import java.util.HashMap;
import java.util.Map;

public abstract class AbstractDiningPhilosophers {
    final int n;
    final boolean nIsOdd;
    long beginDine = Long.MIN_VALUE;
    long endDine = Long.MAX_VALUE;
    final int thinkTimeMilliSec;

    final String impl;

    Map<Integer, Integer> bites = new HashMap<>();

    public AbstractDiningPhilosophers(int n, int thinkTimeMsec, String label) {
        if (n < 1) {
            throw new IllegalArgumentException("n must be positive integer, > 0");
        }
        this.thinkTimeMilliSec = thinkTimeMsec;
        this.n = n;
        this.nIsOdd = (n & 1) != 0;
        this.impl = label;

        int nEven = nIsOdd ? (n / 2) + 1 : n / 2;
        int nOdd = n / 2;
        System.out.printf("%s] init %d diners.  %d odd, %d even\n", impl, n, nOdd, nEven);
    }

    public void initDine(int bitesPer) {
        beginDine = System.nanoTime();

        bites.clear();
        for (int i = 0; i < n; ++i) {
            bites.put(i, bitesPer);
        }
    }

    public void closeDine() {
        endDine = System.nanoTime();
        double time = (endDine - beginDine) * 1E-6;//time in ns * (1E-9 sec/ns) * (1 msec)/1E-3sec
        System.out.printf("%20s] n=%7d, %10.2f msec\n", impl, n, time);
    }

    public abstract void dine(int bites, int timeoutSec) throws InterruptedException;
}
