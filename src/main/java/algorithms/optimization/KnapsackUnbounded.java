package algorithms.optimization;

public class KnapsackUnbounded {

    public static int findMaxValueForCapacity(int[] values, int[] weights, int capacity) {
        int[] memo = new int[capacity + 1];
        int n = values.length;
        int take, i, wIdx;
        for (int wc = 1; wc <= capacity; wc++) { // starting at wc=1 is correct
            for (i = 0; i < n; i++) {
                wIdx = wc - weights[i];
                if (wIdx >= 0) {
                    take = memo[wIdx] + values[i];
                    if (take > memo[wc]) {
                        memo[wc] = take;
                    }
                }
            }
        }
        return memo[capacity];
    }
}
