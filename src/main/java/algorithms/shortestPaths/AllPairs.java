package algorithms.shortestPaths;

public class AllPairs {

    /*
     * goal is to make a runtime complexity efficient use of Dijkstra's with each vertex
     * as a start node, and reuse cached paths and partial paths.
     *
     * the most re-use/cacheing can occur if several of the largest separation pairs,
     * distributed around the convex hull of the data, are solved first.
     * The convex hull can be made in linear runtime complexity.
     *
     * without cacheing, the algorithm runtime complexity is:
     * |V| * O(|V|*log(|V|) + |E|*log(|V|)) for a min heap = min priority queue,
     * or |V| * O(|V|*(C) + |E|*C) for a min heap = yfasttrie with
     * C = log(log(max distance)).
     * (Note: Can only use the yfasttrie on long or integer.  can use Fibonacci Heap on floating point or integer).
     *
     * with cacheing, we store all partial paths and costs of the shortest path.
     * That is, we cache dijkstra shortest path segments path[0] to path[i] and that cost in our distance map.
     * If the graph is an undirected graph, we can also do the reverse in the
     * shortest path, that is, cache path[destination] to path[i] of the reverse
     * path.
     *
     * Using the convex hull to solve all pairs longest paths first may lead to
     * best cacheing, not knowing the point distribution ahead of time,
     * then subsequent all pairs will be using mostly cached values.
     *
     *  e.g. a dense grid of 7 x 7 points:
     *
     * 7  *               *
     * 5    *           *
     * 4       *     *
     * 3          *
     * 2       *     *
     * 1    *           *
     * 0 *                 *
     *   0  1  2  3  4  5  6
     * for undirected graph:  hull is 28 points, left and bottom half only is 14 points.
     *          dense grid: no caching, we have  all srcs=49, all src to all dest is approx (49)/2 ~ 25 dks runs.
     *          dense grid: all combinations of hull start to end for half of hull ~ 14/2 ~ 7 dks runs.
     *              cache all partial paths.  assume that the convex hull of src to dest covers
     *              most of the vertexes that are not a geometric equivalent of courts, culdesacs and dead ends, etc.
     *              then the subsequent iteration over all src and dest will use mostly cached values
     *              and need very few dks iterations.
     *
     * The resulting distance map might be very large for large graphs.
     * A large graph might need to be segmented into smaller pieces and handled in parallel.
     *
     * This algorithm should be useful in building more efficient graph algorithms.
     */
}
