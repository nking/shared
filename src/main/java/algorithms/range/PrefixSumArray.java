package algorithms.range;

public class PrefixSumArray {
    // length is a.length + 1
    // array[0] = 0.
    // a query for [i0,i1] is array[i1+1] - array[i0]
    // e.g. [1,2] = array[2+1] - array[1]
    protected long[] array;

    public PrefixSumArray(long[] a) {
        array = new long[a.length + 1];
        System.arraycopy(a, 0, array, 1, a.length);
        for (int i = 1; i < array.length; ++i) {
            array[i] += array[i-1];
        }
    }

    public long[] sum(int[][] queries, boolean queriesAre0Based) {
        long[] out = new long[queries.length];
        int i0, i1;
        for (int i = 0; i < queries.length; ++i) {
            int[] q = queries[i];
            i0 = q[0];
            i1 = q[1];
            if (queriesAre0Based) {
                --i0;
                --i1;
            }
            //[i0,i1] is array[i1+1] - array[i0]
            out[i] = array[i1 + 1] - array[i0];
        }
        return out;
    }
}
