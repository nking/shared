package algorithms.concurrency.diningPhilosophers;

public interface IDiningPhilosophers {
    /**
     * begin dining
     * @param bites
     * @param timeoutSec
     * @throws InterruptedException
     */
    public void dine(int bitesPer, int timeoutSec) throws InterruptedException;
}
