package algorithms.bipartite;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;

import java.util.*;
import java.util.stream.IntStream;

/**
 * find the max value weighted match between vertexes A and vetexes B
 * where A are the row indexes of the weight matrix and B are the column
 * indexes of the weight matrix.
 *
 * The current implementation runtime complexity (r.t.c.)
 * is O(n^4) where |V| = n/2 and |E| = n^2., though
 * it has modified to a r.t.c. of O(n^3) as following Cormen, Leiserson, Rivest and Stein,
 * but the r.t.c. still looks like it could be O(n^4) at worse.
 *
 * Note that maximizing the weight of a matching is the dual of
 * minimizing the sum of feasible vertex labels.
 *
 */
public class Hungarian {

    final protected int nL;
    final protected int nR;
    final protected Map<Integer, Integer> matchL = new HashMap<>();
    final protected Map<Integer, Integer> matchR = new HashMap<>();
    final protected Set<Integer> unmatchedL;
    final protected Set<Integer> unmatchedR;
    final protected double[] labelsL;
    final protected double[] labelsR;
    // note: edges in matchedL are not in adjL
    final protected Map<Integer, Set<Integer>> adjL;
    final protected double[][] weights;
    final protected boolean transposed;

    // for BFS search
    final protected Set<Integer> fL;
    final protected Set<Integer> fR;
    final protected Set<Integer> rSigmaMinR;

    final protected double tol = 1E-10;
    /*
    E_h is adjL and adjR as undirected edges

    E_M,h is:
       edges in adjL that are not in matchL
           edges in matchR
    */

    /**
     *
     * @param weights [n X m] matrix of weights between vertices of row indices and
     *                vertices of column indices.  use 0. if there is no association.
     */
    public Hungarian(double[][] weights) {

        if (weights.length < weights[0].length) {
            transposed = true;
            this.weights = MatrixUtil.transpose(weights);
        } else {
            this.weights = weights;
            transposed = false;
        }

        this.nL = weights.length;
        this.nR = weights[0].length;
        unmatchedL = new HashSet<>();
        unmatchedR = new HashSet<>();
        labelsL = new double[nL];
        labelsR = new double[nR];
        adjL = new HashMap<>();

        this.fL = new HashSet<>();
        this.fR = new HashSet<>();
        this.rSigmaMinR = new HashSet<>();
    }

    protected void init() {
        // init labels
        for (int vL = 0; vL < this.nL; ++vL) {
            this.labelsL[vL] = Arrays.stream(weights[vL]).max().getAsDouble();
        }
        // labelsR is already 0

        IntStream.range(0, nL).forEach(unmatchedL::add);
        IntStream.range(0, nR).forEach(unmatchedR::add);

        // init E_h with the equality subgraph for the L to R edges
        //l.h+r.h=w(l,r)}
        double labelL, labelR;
        for (int vL = 0; vL < this.nL; ++vL) {
            labelL = labelsL[vL];
            for (int vR = 0; vR < this.nR; ++vR) {
                if (weights[vL][vR] == 0.) {
                    continue;
                }
                labelR = labelsR[vR];
                if (Math.abs(weights[vL][vR] - (labelL + labelR)) < tol) {
                    adjL.putIfAbsent(vL, new HashSet<Integer>());
                    adjL.get(vL).add(vR);
                }
            }
        }

        // init G_M,h with greedy bi-partite matching
        greedyBiPartiteMatching();
        assert(!matchL.isEmpty());
        assert(!matchR.isEmpty());
        /*
        E_M,h is:
           edges in adjL that are not in matchL
           edges in matchR
        */
    }

    protected int[][] solve() {
        if (!matchL.isEmpty() || !matchR.isEmpty()) {
            throw new IllegalArgumentException("already solved");
        }

        init();

        //printDebug("init");

        // perfect matching is that all of L or all of R are matched.
        // not (x==m || y==n)
        // x!=m && y!=n
        while ((matchL.size() != this.nL) && (matchR.size() != this.nR)) {

            int[] revP = findAugmentingPath();
            assert(revP != null);
            assert(revP.length >= 2);

            /*M = M (+) P symmetric difference
                = elements that belong to M or P but not to both.
                no edges in BFS forest fL,fR leave G_M,h.
                no edges in M (matchedL, matchedR) leave G_M,h.
                the path may add new edges and unset matches.
            line 10: update G_h and G_M,h
             */

            int left = 0;
            for (int i = 0; i < revP.length-1; i+=2) {
                int v = revP[i];
                int u = revP[i+1];
                if (left == 0) {// v is Right, u is Left
                    if (matchL.containsKey(u)) {
                        // remove current match for u, but store it in adjL to keep it in G_M,h
                        int tmpR = matchL.get(u);
                        adjL.putIfAbsent(u, new HashSet<Integer>());
                        adjL.get(u).add(tmpR);
                        matchR.remove(tmpR);
                        unmatchedR.add(tmpR);
                    }
                    matchL.put(u, v);
                    matchR.put(v, u);// this adds a right link in G_M,h
                    unmatchedR.remove(v);
                    unmatchedL.remove(u);
                    if (adjL.containsKey(u)) {
                        // remove from adjL so it won't be searched in BFS
                        adjL.get(u).remove(v);
                        if (adjL.get(u).isEmpty()) {
                            adjL.remove(u);
                        }
                    }
                } /*else { // v is Left, u is Right
                    // edge already exists
                    assert(adjL.containsKey(v) && adjL.get(v).contains(u));
                }
                left ^= 1;*/
            }

            //printDebug("after update G_M,h for aug path");
        }

        int nRes = Math.min(this.nL, this.nR);
        assert(matchL.size() == nRes);
        assert(matchR.size() == nRes);
        int[][] matched = new int[nRes][];
        int i = 0;
        if (transposed) {
            for (int vL : matchL.keySet()) {
                matched[i++] = new int[]{matchL.get(vL), vL};
            }
        } else {
            for (int vL : matchL.keySet()) {
                matched[i++] = new int[]{vL, matchL.get(vL)};
            }
        }
        return matched;
    }

    protected int[] findAugmentingPath() {
        // forest of nodes.  used to mark "visited" also.
        fL.clear();
        fR.clear();

        Queue<Integer> q = new ArrayDeque<>();
        // 1 if is Left, else 0 if is Right
        Queue<Integer> qL = new ArrayDeque<>();

        Map<Integer, Integer> prevL = new HashMap<>();
        Map<Integer, Integer> prevR = new HashMap<>();

        // init queue:
        for (int vL : unmatchedL) {
            q.add(vL);
            qL.add(1);
            fL.add(vL);
        }

        int nIter = 0;

        int augR = -1;
        // repeat until an augmenting path is found. aug path is when reach an unmatched right vertex
        while (fL.size() < nL || fR.size() < nR) {
            if (q.isEmpty()) {

                augmentLabelsAndGMh();

                fL.clear();
                fR.clear();
                for (int vL : unmatchedL) {
                    q.add(vL);
                    qL.add(1);
                    fL.add(vL);
                }
                ++nIter;
            } // end if empty queue

            augR = bfs(q, qL, prevL, prevR);

            if (augR > -1) {
                break;
            }
        }
        if (augR == -1) {
            throw new IllegalArgumentException("error in algorithm impl or perfect match given the weights is not feasible");
        }

        // construct the augmenting path to return
        int v = augR;
        int left = 0;
        int[] path = new int[2*Math.min(nL,nR)];
        int i = 0;
        path[i] = v;
        int u;
        Map<Integer, Integer> prev;
        while (true) {
            prev = (left == 0) ? prevR : prevL;
            if (!prev.containsKey(v)) {
                break;
            }
            u = prev.get(v);
            left ^= 1;
            path[++i] = u;
            v = u;
        }

        return Arrays.copyOf(path, i+1);
    }

    protected int bfs(Queue<Integer> q, Queue<Integer> qL, Map<Integer, Integer> prevL, Map<Integer, Integer> prevR) {
        //BFS search
        int u = q.poll();
        boolean isL = (qL.poll() == 1);

        int augR = -1;

        Set<Integer> adj = getAdj(u, isL);
        if (adj == null) {
            return augR;
        }
        for (int nei : adj) {
            if (!isL) { // nei is in Left
                if (fL.contains(nei)) {
                    continue;
                }
                prevL.put(nei, u);
                fL.add(nei);
                q.add(nei);
                qL.add(1);
            } else { // nei is in Right
                if (fR.contains(nei)) {
                    continue;
                }
                prevR.put(nei, u);
                fR.add(nei);
                if (unmatchedR.contains(nei)) {//found an unmatched R
                    augR = nei;
                    break;
                }
                q.add(nei);
                qL.add(0);
            }
        } // end loop over neighbors
        return augR;
    }

    protected void augmentLabelsAndGMh() {
        double delta = calcDelta(rSigmaMinR);
        assert(delta > 0.);

        for (int vL : fL) {
            labelsL[vL] -= delta;
        }
        for (int vR : fR) {
            labelsR[vR] += delta;
        }

        //debugFeasibility();

        // line 15: update G_M,h
        // E_M,h is: edges in adjL that are not in matchL
        //           edges in matchR
        // only need to recalc G_M,h for rows in rSigmaMinR
        //TODO: follow up on reducing rSigmaMin to make this O(n) r.t.c.
        for (int vR : rSigmaMinR) {
            for (int vL = 0; vL < this.nL; ++vL) {
                if (matchL.containsKey(vL) && matchL.get(vL) == vR) {
                    continue;
                }
                if (Math.abs((labelsL[vL] + labelsR[vR]) - weights[vL][vR]) < tol) {
                    // it is not a matched pair which was flagged above.
                    adjL.putIfAbsent(vL, new HashSet<Integer>());
                    adjL.get(vL).add(vR);
                } else {
                    if (adjL.containsKey(vL) && adjL.get(vL).contains(vR)) {
                        adjL.get(vL).remove(vR);
                    }
                }
            }
        } // end loop over rSigmaMinR
    }

    private void debugFeasibility() {
        //check that l.h + r.h >= w[l][r]
        for (int vL = 0; vL < nL; ++vL) {
            for (int vR = 0; vR < nR; ++vR) {
                if ((labelsL[vL] + labelsR[vR]) < (weights[vL][vR] - tol)) {
                    printDebug("ERROR");
                    throw new IllegalArgumentException("labels are not feasible for vL=" + vL + ", vR=" + vR);
                }
            }
        }
    }

    private void printDebug(String label) {
        System.out.printf("\n%s\n", label);

        System.out.printf("W=\n%s\n", FormatArray.toString(weights, "%.4f"));

        System.out.printf("L labels = %s\n", FormatArray.toString(labelsL, "%.5f"));
        System.out.printf("R labels = %s\n", FormatArray.toString(labelsR, "%.5f"));

        System.out.printf("matched L: ");
        for (int vL : matchL.keySet()) {
            System.out.printf("%d:%d, ", vL, matchL.get(vL));
        }
        System.out.println();

        System.out.printf("matched R: ");
        for (int vR : matchR.keySet()) {
            System.out.printf("%d:%d, ", vR, matchR.get(vR));
        }
        System.out.println();

        System.out.printf("f_L: ");
        for (int vL : fL) {
            System.out.printf("%d, ", vL);
        }
        System.out.println();

        System.out.printf("f_R: ");
        for (int vR : fR) {
            System.out.printf("%d, ", vR);
        }
        System.out.println();

        System.out.printf("rSigmaMinR: ");
        for (int vR : rSigmaMinR) {
            System.out.printf("%d, ", vR);
        }
        System.out.println();

        System.out.printf("adjL:\n");
        for (int vL : adjL.keySet()) {
            System.out.printf("%d: ", vL);
            for (int vR : adjL.get(vL)) {
                System.out.printf("%d, ", vR);
            }
            System.out.println();
        }

        /*System.out.printf("adjR:\n");
        for (int vR : adjR.keySet()) {
            System.out.printf("%d: ", vR);
            for (int vL : adjR.get(vR)) {
                System.out.printf("%d, ", vL);
            }
            System.out.println();
        }*/

        System.out.printf("unmatched R: ");
        for (int vR : unmatchedR) {
            System.out.printf("%d, ", vR);
        }
        System.out.println();

        System.out.printf("unmatched L: ");
        for (int vL : unmatchedL) {
            System.out.printf("%d, ", vL);
        }
        System.out.println();
    }

    /**
     * calculate delta and store the rows with rSigma == delta when nIterBFSAug>0
     * else stores rows with rSigma <= delta.
     * @param rSigmaMinR
     * @param nIterBFSAug
     * @return
     */
    protected double calcDelta(Set<Integer> rSigmaMinR) {

        /*
        we calc min(delta) from the set R-fR and fL.

        if rSigmaMinR is empty, we add to rSigmaMinR all R with delta <= min
        else we add to rSigmaMinR all (R-fR) with delta == min,
        */

        double min = Double.POSITIVE_INFINITY;
        // calculate min for {set R - set fR}
        for (int vR = 0; vR < nR; ++vR) {
            if (fR.contains(vR)) continue;
            for (int vL : fL) {
                min = Math.min(min, labelsL[vL] + labelsR[vR] - weights[vL][vR]);
            }
        }

        int sigmaMinIdx = -1;
        if (rSigmaMinR.isEmpty()) {
            for (int vR = 0; vR < nR; ++vR) {
                double m = Double.POSITIVE_INFINITY;
                for (int vL : fL) {
                    m = Math.min(m, labelsL[vL] + labelsR[vR] - weights[vL][vR]);
                }
                if (m <= min) {
                    rSigmaMinR.add(vR);
                    sigmaMinIdx = vR;
                }
            }
        } else {
            for (int vR = 0; vR < nR; ++vR) {
                if (fR.contains(vR)) continue;
                double m = Double.POSITIVE_INFINITY;
                for (int vL : fL) {
                    m = Math.min(m, labelsL[vL] + labelsR[vR] - weights[vL][vR]);
                }
                if (m <= min) {
                    rSigmaMinR.add(vR);
                    sigmaMinIdx = vR;
                }
            }
        }

        //System.out.printf("R sigmaMin=%d, delta=%f\n", sigmaMinIdx, min);

        return min;
    }

    // get the neighbors of u
    protected Set<Integer> getAdj(int u, boolean isL) {
        // E_M,h is: edges in adjL that are not in matchL
        //           edges in matchR
        if (isL) {
            return adjL.get(u);
        } else {
            //TODO: this could be improved at invocation to avoid creating another hashset
            if (matchR.containsKey(u)) {
                Set<Integer> adj = new HashSet<>();
                adj.add(matchR.get(u));
                return adj;
            }
            return null;
        }
    }

    protected void greedyBiPartiteMatching() {
        for (int vL = 0; vL < this.nL; ++vL) {
            if (!adjL.containsKey(vL)) {
                continue;
            }
            Set<Integer> rm = new HashSet<>();
            for (int vR : adjL.get(vL)) {
                if (unmatchedR.contains(vR)) {
                    unmatchedR.remove(vR);
                    matchL.put(vL, vR);
                    matchR.put(vR, vL);
                    unmatchedL.remove(vL);

                    // remove matchL from adjL to make ir usable as part of G_M,h
                    rm.add(vR);
                }
            }
            adjL.get(vL).removeAll(rm);
            if (adjL.get(vL).isEmpty()) {
                adjL.remove(vL);
            }
        }
    }
}
