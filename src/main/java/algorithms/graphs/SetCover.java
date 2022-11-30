package algorithms.graphs;

import algorithms.misc.Misc0;
import algorithms.optimization.LinearProgramming;
import algorithms.util.FormatArray;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 Given a set of elements {1,2,...,n} (called the universe) 
 * and a collection S of m sets whose union equals the universe, 
 * the set cover problem is to identify the smallest sub-collection of 
 * S whose union equals the universe. For example, consider the universe 
 * U={1,2,3,4,5} and the collection of sets S={1,2,3},{2,4},{3,4},{4,5}.
 * Clearly the union of S is U. However, we can cover all of the elements 
 * with the following, smaller number of sets: {1,2,3},{4,5}.
 
 Set Cover is one of Karp's 21 NP-complete problems shown to be NP-complete in 1972.
  
 A set cover is the smallest vertex cover of a graph,
 where a vertex cover is a subset of a graph's vertices which represents at 
 * least one vertex from every edge in the full graph
 * 
 * <pre>
 * https://en.m.wikipedia.org/wiki/Set_cover_problem
 * </pre>
 * @author nichole
 */
public class SetCover {
    
    private final Random rand = Misc0.getSecureRandom();
    
    /**
     * machine precision used in evaluating whether a double is different
     * from another double.
     */
    private static double eps = 1e-11;
    
    public SetCover() {
        this(System.nanoTime());
    }
    
    public SetCover(long randomSeed) {
        //seed = 180328550254112L;
        System.out.println("seed=" + randomSeed);
        System.out.flush();
        rand.setSeed(randomSeed);
    }
    
    /**
     * find a minimum weighted set cover using a randomized rounding approximation algorithm
     * of 2*log_e(n) where n is the number of vertexes in the final
     * cover (== nU).
     * TODO: revisit this for element weights.
     * The cost of the sets in the cover is minimized.
     * The problem is NP-complete.
     * @param nU the number of items in U. U is the universe of elements in sets.
     * items in U are the sequential range of integers from 0 to nU-1.
     * @param sets a list of sets for which each set contains integers from the range
     * 0 through weights.length - 1, inclusive.
     * @param weights the weights of each element.
     * @return the list of indexes of sets which comprise the cover, that is the
     * indexes of the minimum subset of sets that together include all numbers 
     * 0 through nU-1, inclusive.
     */
    public TIntSet weightedElementsApprox2LgN(int nU, List<TIntSet> sets, double[] weights) {
        
        double[] sWeights = calcSetWeightsFromElementWeights(nU, sets, weights);

        return weightedSetsApprox2LgN(nU, sets, sWeights);
    }
        
    /**
     * find a minimum weighted set cover using a randomized rounding approximation algorithm
     * of 2*log_e(n) where n is the number of vertexes in the final
     * cover (== nU).
     * The cost of the sets in the cover is minimized.
     * The problem is NP-complete.
     * @param nU the number of items in U. U is the universe of elements in sets.
     * items in U are the sequential range of integers from 0 to nU-1.
     * @param sets a list of sets for which each set contains integers from the range
     * 0 through weights.length - 1, inclusive.
     * @param weights the weights of each set in sets.
     * @return the list of indexes of sets which comprise the cover, that is the
     * indexes of the minimum subset of sets that together include all numbers 
     * 0 through nU-1, inclusive.
     */
    public TIntSet weightedSetsApprox2LgN(int nU, List<TIntSet> sets, double[] weights) {
        /*
        material from lecture slides of Principal lecturer: Dr Thomas Sauerwald
        Advanced Algorithms, University of Cambridge.
        VII. Approximation Algorithms: Randomisation and Rounding
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/materials.html
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/rand.pdf
        
        for the linear program:
            minimize: 
                summation_S_in_Cover( c(S) ) = summation_S_in_sets( c(S)*y(S) )
            subject to:
                summation_S_in_sets : x in S ( y(S) ) >= 1
                y(S) <= 1 for each S in sets ----\
            non-negativity constraints:           \ these 2 rules are derived from y(S) ∈ [0,1]
                y(S) >= 0 for each S in sets ----/
        
        for the weighted set cover w/ LP(X, F, c):
            // where each x in X belongs to at least 1 subset in F
            // F is the list of subsets to choose from when building the cover
            // c is the cost for each set in F
            compute y, an optimal solution to the linear program
            C = empty set
            repeat 2*ln(nU) times
                for each S in F
                    let C = C union S with probabilty y(S)
            return C
                
        misc notes about randomized picking of set S given probability y(S):
            http://theory.stanford.edu/~trevisan/cs261/lecture08.pdf
            probability ≥ 1 − 1/e that at least one subset S covers an element u in U (U=[0:nU-1, incl]).
            For each iteration of the while loop, there is a probability at most 
                1/e that element u in U is not covered by the sets added to C 
                in that iteration.
            The probability that u is not covered after ln |U| + k iterations 
               is then at most = (i/|U|(*e^(-k)).
            So then The probability that, after ln|U|+k iterations, there is
              an element that is not covered, is at most the sum over all u 
              of the probability that u is not covered, which is at most e^(−k)
        NOTE: nU = |U|
        
        can compare this cost O(log n)*OPT_LP  to a greedy weighted set cover H_d*OPT_LP 
        where H_d is harmonic series ~ 0.5+ln(d)
        Augmented Greedy Algorithm of weighted Set Cover: Covered = ∅;
            while Covered ̸= U do
                j ← argmin_k( w_k / |S_k ∩ Uncovered| )
                if i is uncovered and i ∈ S_j, set pi = Covered = Covered ∪ S_j;
                   A = A ∪ {j}.
            end while;
        Output sets in A as cover
        
        */
                
        LinearProgramming.StandardForm standForm = 
            createLinearProgramInStandardFormForWeightedSets(nU, sets, weights);
        LinearProgramming lp = new LinearProgramming();
        LinearProgramming.SlackForm lpSoln = lp.solveUsingSimplexMethod(standForm);
        double[] primalX = lpSoln.calculatePrimalX();
        double[] dualY = lpSoln.calculateDualY();
        System.out.printf("primalX=%s\n", FormatArray.toString(primalX, "%.3f"));
        System.out.printf("dualY=%s\n", FormatArray.toString(dualY, "%.3f"));
        
        return weightedApprox2LgN(nU, sets, weights, primalX);
    }
    
    /**
     * find a minimum weighted set cover using a randomized rounding approximation algorithm
     * of 2*log_e(n) where n is the number of vertexes in the final
     * cover (== nU).
     * The cost of the sets in the cover is minimized.
     * The problem is NP-complete.
     * @param nU the number of items in U. U is the universe of elements in sets.
     * items in U are the sequential range of integers from 0 to nU-1.
     * @param sets a list of sets for which each set contains integers from the range
     * 0 through weights.length - 1, inclusive.
     * @param weights the weights of each set in sets.
     * @param yProbabilities coefficients of the optimal solution to LinearProgram
     * for the given sets and weights.
     * @return the list of indexes of sets which comprise the cover, that is the
     * indexes of the minimum subset of sets that together include all numbers 
     * 0 through nU-1, inclusive.
     */
    private TIntSet weightedApprox2LgN(int nU, List<TIntSet> sets, 
        double[] weights, double[] yProbabilities) {
        
        /*
        material from lecture slides of Principal lecturer: Dr Thomas Sauerwald
        Advanced Algorithms, University of Cambridge.
        VII. Approximation Algorithms: Randomisation and Rounding
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/materials.html
        https://www.cl.cam.ac.uk/teaching/1617/AdvAlgo/rand.pdf
        
        for the linear program:
            minimize: 
                summation_S_in_Cover( c(S) ) = summation_S_in_sets( c(S)*y(S) )
            subject to:
                summation_S_in_sets : x in S ( y(S) ) >= 1
                y(S) <= 1 for each S in sets ----\
            non-negativity constraints:           \ these 2 rules are derived from y(S) ∈ [0,1]
                y(S) >= 0 for each S in sets ----/
        
        for the weighted set cover w/ LP(X, F, c):
            // where each x in X belongs to at least 1 subset in F
            // F is the list of subsets to choose from when building the cover
            // c is the cost for each set in F
            compute y, an optimal solution to the linear program
            C = empty set
            repeat 2*ln(nU) times
                for each S in F
                    let C = C union S with probabilty y(S)
            return C
                
        misc notes about randomized picking of set S given probability y(S):
            http://theory.stanford.edu/~trevisan/cs261/lecture08.pdf
            probability ≥ 1 − 1/e that at least one subset S covers an element u in U (U=[0:nU-1, incl]).
            For each iteration of the while loop, there is a probability at most 
                1/e that element u in U is not covered by the sets added to C 
                in that iteration.
            The probability that u is not covered after ln |U| + k iterations 
               is then at most = (i/|U|(*e^(-k)).
            So then The probability that, after ln|U|+k iterations, there is
              an element that is not covered, is at most the sum over all u 
              of the probability that u is not covered, which is at most e^(−k)
        NOTE: nU = |U|
        
        can compare this cost O(log n)*OPT_LP  to a greedy weighted set cover H_d*OPT_LP 
        where H_d is harmonic series ~ 0.5+ln(d)
        Augmented Greedy Algorithm of weighted Set Cover: Covered = ∅;
            while Covered ̸= U do
                j ← argmin_k( w_k / |S_k ∩ Uncovered| )
                if i is uncovered and i ∈ S_j, set pi = Covered = Covered ∪ S_j;
                   A = A ∪ {j}.
            end while;
        Output sets in A as cover
        
        */
               
        // TODO: consider visiting the sets ordered by highest probability, 
        //      though this is increasingly not the probability of the remaining 
        //      members summed.
          
        int i;
        TIntSet indexes = new TIntHashSet();
        for (i = 0; i < sets.size(); ++i) {
            indexes.add(i);
        }
        TIntIterator iter;
        TIntSet rm = new TIntHashSet();
        double r;
        
        TIntSet c = new TIntHashSet();
        int nMax = (int)Math.ceil(2*Math.log(nU));
        int idx;
        for (i = 0; i < nMax; ++i) {
            iter = indexes.iterator();
            while (iter.hasNext()) {
                idx = iter.next();
                r = rand.nextDouble();
                //System.out.printf("%d) r=%.3f, yP[%d]=%.3f\n", i, r, idx, yProbabilities[idx]);
                if (r <= yProbabilities[idx]) {
                    c.add(idx);
                    rm.add(idx);
                }
            }
            iter = rm.iterator();
            while (iter.hasNext()) {
                idx = iter.next();
                indexes.remove(idx);
            }
            if (c.size() == nU) {
                break;
            }
        }
        return c;
    }
    
    /**
     * find a set cover that is O(log-n)-approx, that is no more than O(log-n) times
     * as large as the optimal set cover
     * (e.g., for n=100, this would be up to 6.64 times as large as the optimal solution).
     * This is a greedy approach.
     * The runtime complexity is polynomial.
     * <pre>
     * The algorithm implements pseudocode from
     *  from https://www.ics.uci.edu/~goodrich/teach/graph/notes/Approximation.pdf
     * and Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms" chap 35.3.
     * </pre>
     * @param sets
     * @return the list of indexes of sets which comprise the cover.
     */
    public TIntList approxLgN(List<TIntSet> sets) {
        
        // make a copy of the sets to edit it
        TIntObjectMap<TIntSet> setsMap = copy(sets);
        TIntObjectIterator<TIntSet> iter;
        
        TIntList c = new TIntArrayList();
        
        TIntSet si, siMaxN = null;        
        int maxN, i, n, maxNIdx = -1, idx;
        TIntList rm = new TIntArrayList();
        while (!setsMap.isEmpty()) {
            maxN = Integer.MIN_VALUE;
            iter = setsMap.iterator();
            for (i = 0; i < setsMap.size(); ++i) {
                iter.advance();
                si = iter.value();
                idx = iter.key();
                n = si.size();
                //System.out.printf("   idx=%d (%s) n=%d\n", idx, si.toString(), n);
                if (n > 0 && n > maxN) {
                    maxN = n;
                    maxNIdx = idx;
                }
            }
            if (maxN == Integer.MIN_VALUE) {
                break;
            }
            siMaxN = setsMap.remove(maxNIdx);
            //System.out.printf("max idx=%d (%s)\n", maxNIdx, siMaxN.toString());
            c.add(maxNIdx);
            
            // remove siMaxN from each set in s.  and store the empty sets as indexes to remove after use of iterator
            iter = setsMap.iterator();
            for (i = 0; i < setsMap.size(); ++i) {
                iter.advance();
                si = iter.value();
                si.removeAll(siMaxN);
                if (si.isEmpty()) {
                    rm.add(iter.key());
                }
            }
            for (i = 0; i < rm.size(); ++i) {
                setsMap.remove(rm.get(i));
            }
            rm.clear();
        }
        
        return c;
    }

    protected TIntObjectMap<TIntSet> copy(List<TIntSet> s) {
        TIntObjectMap<TIntSet> c = new TIntObjectHashMap<TIntSet>(s.size());
        int i;
        TIntSet cSet, si;
        TIntIterator iter;
        for (i = 0; i < s.size(); ++i) {
            cSet = new TIntHashSet();
            si = s.get(i);
            iter = si.iterator();
            while (iter.hasNext()){
                cSet.add(iter.next());
            }
            c.put(i, cSet);
        }
        return c;
    }

    protected static LinearProgramming.StandardForm 
        createLinearProgramInStandardFormForWeightedSets(
        int nX, List<TIntSet> sets, double[] weights) {
        
        /*
        For the primal:
        the 'a' matrix will have rows that represent each x 0 through nX-1
        while the columns will represent whether the set within sets contains
        that integer x.
        
        e.g.  List<TIntSet> sets = [ {0, 2, 3}, {1, 2} ];
              nX = 4
              would create:
              
              a = | x=0 row | = | 1  0 | because 0 is in sets[0] but not in sets[1]
                  | x=1 row |   | 0  1 |  ...    1 ...
                  | x=2 row |   | 1  1 |  ...    2 ...
                  | x=3 row |   | 1  0 |  ...    1 ...
                  | 2nd set of constraints |
        
        minimize: 
                summation_S_in_Cover( c(S) ) = summation_S_in_sets( c(S)*y(S) )
            subject to:
                summation_S_in_sets : x in S ( y(S) ) >= 1
                y(S) <= 1 for each S in sets ----\
            non-negativity constraints:           \ these 2 rules are derived from y(S) ∈ [0,1]
                y(S) >= 0 for each S in sets ----/

        a length is nU + nS.  a[0] length is nS.
        
        solves for y.         
        */
        
        int nS = sets.size();
        
        double[] c = Arrays.copyOf(weights, weights.length);
        double[] b = new double[nX + nS];
        Arrays.fill(b, 1);
        
        double[][] a = new double[nX + nS][nS];
        int i, j;
        for (i = 0; i < nX; ++i) {
            a[i] = new double[nS];
            for (j = 0; j < sets.size(); ++j) {
                if (sets.get(j).contains(i)) {
                    a[i][j] = 1;
                }
            }
        }
        int i2;
        for (i = 0, i2=nX; i < nS; ++i, ++i2) {
            a[i2] = new double[nS];
            a[i2][i] = 1;
        }
                
        boolean isMaximization = false;
        int[] constraintComparisons = new int[nX + nS];
        Arrays.fill(constraintComparisons, 0, nX, 1);
        Arrays.fill(constraintComparisons, nX, nX+nS, -1);
        
        boolean[] nonnegativityConstraints = new boolean[nS];
        Arrays.fill(nonnegativityConstraints, true);
        
        LinearProgramming.StandardForm standForm = LinearProgramming
            .convertLinearProgramToStandardForm(isMaximization, a, b, c, 
            constraintComparisons, nonnegativityConstraints);
        
        //System.out.printf("Linear Program in standard form=\n%s\n", standForm.toString());

        return standForm;
    }

    public static double[] calcSetWeightsFromElementWeights(int nU,
        List<TIntSet> sets, double[] weights) {
        double[] weightsSets = new double[sets.size()];
        int i, idx;
        TIntSet set;
        TIntIterator iter;
        double sum = 0;
        for (i = 0; i < sets.size(); ++i) {
            set = sets.get(i);
            iter = set.iterator();
            sum = 0;
            while (iter.hasNext()) {
                idx = iter.next();
                sum += weights[idx];
            }
            weightsSets[i] = sum;
        }
        return weightsSets;
    }
    
}
