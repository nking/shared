package algorithms.misc;

/**
 *
 * @author nichole
 */
public class PseudoRNG {

    //TODO: add tests from
    // https://github.com/nmondal/diehard.c
    // https://web.archive.org/web/20160125103112/http:/stat.fsu.edu/pub/diehard/

    /**
     * a fast pseudo-random number generator suitable for use in testing.
     <pre>
     references: Marsaglia 2003, Journal of Statistical Software, 8(13),
     "XorShift RNGs"
     referenced by Goetz et al. "Java Concurrency in Practice", Listing 12.4.
     </pre>
     @param y
     @return
     */
    public static int xorShift(int y) {
        y ^= (y << 6);
        y ^= (y >> 21);
        y ^= (y << 7);
        return y;
    }

    /**
     * a fast pseudo-random number generator suitable for use in testing.
     <pre>
     references: Marsaglia 2003, Journal of Statistical Software, 8(13),
     "XorShift RNGs"
     referenced by Goetz et al. "Java Concurrency in Practice", Listing 12.4.
     </pre>
     @param y
     @return
     */
    public static long xorShift(long y) {
        y ^= (y << 13);
        y ^= (y >> 7);
        y ^= (y << 17);
        return y;
    }
}
