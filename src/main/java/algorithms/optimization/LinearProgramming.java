package algorithms.optimization;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.Misc0;
import algorithms.optimization.LinearProgramming.SlackForm.STATE;
import algorithms.sort.MiscSorter;
import algorithms.util.FormatArray;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;

/**
 * A program to find the optimal values of x if possible given
 * an objective and constraints in Standard Form.
  * 
 
 * The Simplex Method is used and will return infeasible or unbounded
 * if an optimal Basic feasible solution is not possible.
 * The runtime complexity is exponential at worse, but often performs as
 * a polynomial runtime complexity.
 * (pivot operations can be expensive for large problems, having combinatorial
 * complexity C(m+n, m) at worse.
 * 
 * Alternative algorithms are ellipsoid, and interior-point.
 * Interior-point is a class of algorithms that solve linear and nonlinear 
 * convex optimization problems in provably polynomial time.
 * https://en.m.wikipedia.org/wiki/Interior-point_method
 * The interior-point algorithms can use linear algebra more effecitvely
 * for large and sparse problems.
 * 
 * when all variables must be integers, the problem is more specifically, Linear Integer
       Programming problem and it is NP-hard, no known polynomial time algorithm.
       
 * <pre>
 * The algorithm implements Chap 29 of Cormen et al. "Introduction to Algorithms"
 * 
  Some definitions:
   
   Standard Form:
       maximization of a linear function subject to linear *inequalities*

   Slack Form:
       maximization of a linear function subject to linear *equalities*

algebraic solution: write the linear program in slack form which express the equalities
        in terms of "basic variables" and "nonbasic variables".
        involves making a basic variable nonbasic and a nonbasic variable basic using
        "pivot" operations.
        
Simplex Method:
      Notation:  variables x1, x2, ... xn
                 a realization of the variables: xHat0, xHat1, ...xHatn
      `
      Standard Form:
          real numbers: c1...cn; b1,...bm; and aij for a=1:m and j=1:n.

                        Find numbers x1,...xn:
        objective ->    maximize summation_j=1:n(cj*xj)
        constraints ->   subject to: summation_j=1:n(aij*xj) .leq. bi for i=1:m
        constraints ->   xj .geq. 0 for j=1:n
                         the later is a nonnegativity constraint
                         
              OR expressed more compactly:
          A = (aij) =  mXn matrix
          b = (bi) an m-dimensional vector
          c = (cj) an n-dimensional vector
          x = (xj) an n-dimensional vector
        objective ->    maximize c^T*x
        constraints ->   subject to: A*x .leq. b
                         x .leq. 0

       A feasible solution is the numbers xHat as x1,...xn that satisfy all constraints.

       An unfeasible solution fails to satisfy all constraints.

       The optimal solution is the feasible solution which maximizes the objective.

       An unbounded solution has feasible solution, but does not have finite optimal objective.
      
 * </pre>
 * <pre>
 * More regarding the Simplex Method runtime complexity:
 * From wikipedia: http://en.wikipedia.org/wiki/Simplex_algorithm
 * 
 * '...Analyzing and quantifying the observation that the simplex algorithm is 
 * efficient in practice, even though it has exponential worst-case complexity, 
 * has led to the development of other measures of complexity. The simplex 
 * algorithm has polyxnomial-time average-case complexity under various 
 * probability distributions, with the precise average-case performance of the 
 * simplex algorithm depending on the choice of a probability distribution for 
 * the random matrices.   Another approach to studying "typical phenomena" uses 
 * Baire category theory from general topology, and to show that (topologically) 
 * "most" matrices can be solved by the simplex algorithm in a polynomial number 
 * of steps. Another method to analyze the performance of the simplex algorithm 
 * studies the behavior of worst-case scenarios under small perturbation – are 
 * worst-case scenarios stable under a small change (in the sense of structural 
 * stability), or do they become tractable?   Formally, this method uses random 
 * problems to which is added a Gaussian random vector 
 * ("smoothed complexity")....'
 * </pre>
 * 
 * NOTE: some Simplex algorithms use a matrix composed of the objective function 
 * (c, v) and the constraints (a, b) in a format called the "tableau". 
 * The tableau format varies in some implementations.
 <pre>
     |  1  (-1*objective function coefficients)         0         |
     |  0  (constraint 1 coefficients)      (constraint 1 bound)  |
     |  0   ... through all constraints          ...              |
     
   
 notes from Linear Programming, CSE 6331 Algorithms Steve Lai:
 
 Another form frequently used for Standard Form tableau:
        |  x1         x2         ...        xn  |
  ------|---------------------------------------|------------
    y1  |  a11        a12      ...      a1n     |  .leq. b1
    y2  |  a21        a22      ...      a2n     |  .leq. b2
    ... |  ...        ...      ...      ...     |  ...
    ym  |  am1        am2      ...      amn     |  .leq. bm
  ------|---------------------------------------|-------------
        |.geq. c1    .geq. c2          .geq. cn |
</pre>
   The written text of the equations is "Dictionary" form.
 <pre>
 The original linear program before transformations is called "primal".
 
 Regarding duality:
 
 A Primal Standard Form LP:
     max: c^T*x
     subject to: A*x .leq. b and x .geq. 0
 A Primal Minimization LP:
     min: y^T*b
     subject to: y^T*A .geq. c^T and y .geq. 0    

Example:
   A Primal Standard Form LP  |  Dual Minimization LP
      max: x1 + x2            |    min: 4*y1 + 12*y2 + y3
      subject to:             |    subject to:
        x1 + 2*x2 .leq. 4     |      y1 + 4*y2 - y3 .geq. 1
      4*x1 + 2*x2 .leq. 12    |     2*y1+ 2*y2 + y3 .geq. 1
       -x1 +   x2 .leq. 1     |
      x1, x2 .geq. 0          |    y1, y2, y3 .geq. 0
 
        |  x1          x2 |
  ------|----------------------------
    y1  |  1           2  | .leq. 4
    y2  |  4           2  | .leq. 12
  ------|-----------------|----------
        |.geq. 1  .geq. 1 |
 </pre>
 The number of basic solutions is at most C(m+n, m).
 <pre>
 Basic feasible solutions are extreme points
     A set C is said to be convex if for all points x1, x2 in C 
         we have lambda*x1 + (1 - lambda)*x2 in C for all 0 .leq. lambda .leq. 1.
     A point x in a convex set C is called an extreme point if there exist no 
         x1, x2 in C, x1 != x2 , and 0 .lt. lambda .lt. 1,
         such that x = lambda*x1 + (1 - lambda)*x2.
     The feasible region of an LP is a convex set.
     Basic feasible solutions of an LP (in the equality form)
         are extreme points of the feasible region.
     If the LP is feasible bounded, then at least one
         optimum solution occurs at an extreme point.
 </pre>
 * @author nichole
 */
public class LinearProgramming {

    private final Random rand = Misc0.getSecureRandom();
    
    /**
     * machine precision used in evaluating whether a double is different
     * from another double.
     */
    private static double eps = 1e-11;
    
    public LinearProgramming() {
        this(System.nanoTime());
    }
    
    public LinearProgramming(long randomSeed) {
        //seed = 180328550254112L;
        System.out.println("seed=" + randomSeed);
        System.out.flush();
        rand.setSeed(randomSeed);
    }
    
    /**
        given a linear program in Standard Form if the problem is unfeasible, 
        this method returns a slackform with state set to unfeasible, else 
        if the problem is unbounded the method returns a slackform with state 
        set to unbounded, else
        returns a slack form for which a basic solution is feasible and optimal.
     <pre>
       The method is implemented from pseudocode in Section 29.3 of Cormen et al.
       
       Some definitions:
           A feasible solution is the numbers xHat as x1,...xn that satisfy all constraints.

           An unfeasible solution fails to satisfy all constraints.

           The optimal solution is the feasible solution which maximizes the objective.

           An unbounded solution has feasible solution, but does not have finite optimal objective.
      
     </pre>
     * @param standForm
     * @return 
     */
    public SlackForm solveUsingSimplexMethod(StandardForm standForm) {
      
        SlackForm slackForm = initializeSimplex(standForm);
        
        //System.out.printf("initialized:\n%s\n", slackForm.toString());
        
        if (!slackForm.state.equals(STATE.FEASIBLE) && !slackForm.state.equals(STATE.OPTIMAL)) {
            System.out.println("exiting: " + slackForm.state.name());
            return slackForm;
        }
                
        //lines 2-11 of the SIMPLEX method of Cormen et al. Chap 29.:
        
        // the j's in N w/ c_j > 0
        TIntSet positiveCIndexes = findPositiveCIndexes(slackForm);
        
        double[] delta = new double[slackForm.b.length];
        
        //choose eIdx from the j's in N w/ c_j > 0
        int eIdx = chooseEnteringIndex(slackForm, positiveCIndexes);
        
        int lIdx = -1, i;
        double minDelta; // lIdx is chosen with this
        
        // iterate until all c's are positive
        
        //choose index lIdx in bIndices that minimizes delta[i]
        while (eIdx > -1) {
                 
            minDelta = Double.POSITIVE_INFINITY; // lIdx is chosen with this
            
            /*
                for each index i in B
                    if a[i][eIdx] > 0
                        delta[i] = b[i]/a[i][eIdx]
                    else 
                        delta[i] = inf.
                choose index lIdx in bIndices that minimizes delta[i]
            */
            
            for (i = 0; i < delta.length; ++i) {
                if (slackForm.a[i][eIdx] > 0) {
                    delta[i] = slackForm.b[i]/slackForm.a[i][eIdx];
                    if (delta[i] < minDelta) {
                        minDelta = delta[i];
                        lIdx = i;
                    }
                } else {
                    delta[i] = Double.POSITIVE_INFINITY;
                }
            }
            
            //System.out.printf("eIdx=%d lIdx=%d minDelta=%.3f\n", eIdx, lIdx, minDelta);
            
            if (Double.isInfinite(minDelta)) {
                System.out.println("exit: UNBOUNDED");
                slackForm.state = STATE.UNBOUNDED;
                return slackForm;
            }
            assert(lIdx != -1);
           
            slackForm = pivot(slackForm, lIdx, eIdx);
            
            slackForm.computeBasicSolution();
            double z = slackForm.evaluateObjective();
            
            //System.out.printf("after pivot of eIdx=%d lIdx=%d, z=%.3f:\n%s\n", 
            //    eIdx, lIdx, z, slackForm.toString());
                        
            positiveCIndexes = findPositiveCIndexes(slackForm);
                        
            eIdx = chooseEnteringIndex(slackForm, positiveCIndexes);            
            //System.out.printf("remaining cIndexes=\n%s\n", positiveCIndexes.toString());
        }
        
        // lines 12-16
        /*
        for i=1:n
            if i is in B
                xHat_i = b_i
            else
                xhat_i = i
        return xHat_1:n
        */
        slackForm.computeBasicSolution();
        slackForm.state = STATE.OPTIMAL;
        
        return slackForm;
    }
    
    /**
        given a linear program in Standard Form if the problem is unfeasible, 
        this method returns a slackform with state set to unfeasible, , else 
        if the problem is unbounded the method returns a slackform with state 
        set to unbounded, else
        returns a slack form for which a basic solution is feasible.
     * <pre>
     * The method is implemented from pseudocode in Section 29.3 of Cormen et al.
     * </pre>
     * @param standForm
     * @return 
     */
    protected SlackForm initializeSimplex(StandardForm standForm) {
        
        SlackForm _slackForm = convertConstraints(standForm);
        //System.out.printf("init: convert to slackForm =\n%s\n", _slackForm.toString());
                
        int m = _slackForm.b.length;
        int n = _slackForm.nIndices.length;
        
        int i;
        int lIdx = findMinIndex(_slackForm.b);
                
        //is the initial basic solution feasible?
        if (lIdx > -1 && _slackForm.b[lIdx] >= 0) {
            int[] nIndices = new int[n];
            int[] bIndices = new int[m];
            //nIndices = [1,2,...n]. NOTE: using 0-based indexes instead
            for (i = 1; i < n; ++i) {
                nIndices[i] = i;
            }
            //bIndices = [n+1,n+2,...n+m]. NOTE: using 0-based indexes instead
            for (i = 0; i < m; ++i) {
                bIndices[i] = n+i;
            }
            
            SlackForm slackForm = new SlackForm(nIndices, bIndices, 
                _slackForm.a, _slackForm.b, _slackForm.c, _slackForm.v);
            
            slackForm.state = STATE.FEASIBLE;
            
            return slackForm;
        }
        
        // arrive here if any of b[i] were negative in the slack form
                       
        /*
        L is in Standard Form.
        Slack Form L_aux is an alteration of L for n+1 variables", following eqns (29.109) - (29.111)
        
        maximize: -x0
        subject to: summation_j=1:n(a_i_j*x_j) - x0 <= b_i for i=1:m
                    x_j >= 0 for j=1:n     
        
        Form L_aux by adding -x0 to the LHS of each equation
            and setting the objective function to -x0
        */
        
        //let (N, B, A, b, c, v) be the resulting Slack Form for L_aux
        SlackForm slackFormAux = createAuxiliarySlackForm(standForm);
        
        //System.out.printf("init: L_aux=\n%s\n", slackFormAux.toString());
        
        //System.out.printf("init: pivot eIdx=%d (=x%d), lIdx=%d (=x%d)\n",
        //    0, slackFormAux.nIndices[0], lIdx, slackFormAux.bIndices[lIdx]);
        
        //L_aux has n+1 nonbasic variables and m basic variables
        //(N, B, A, b, c, v) = pivot((N, B, A, b, c, v, lIdx, 0)
        slackFormAux = pivot(slackFormAux, lIdx, 0);
        
        //the basic solution is now feasible for L_aux
        double[] xBasicSoln = slackFormAux.computeBasicSolution();
        
        //System.out.printf("init: L_aux after pivot=\n%s\n",
        //    slackFormAux.toString());
        
        boolean isFeasible = slackFormAux.isFeasible();
        //System.out.printf("init: L_aux basicSoln=%s\n", FormatArray.toString(xBasicSoln, "%.3f"));
        //System.out.println("   isFeasible=" + isFeasible);
        slackFormAux.state = STATE.FEASIBLE;
         
        /*        
        iterate the while loop of lines 2-11 of the SIMPLEX method
            until an optimal solution to L_aux is found
        
        OPTIMAL FEASIBLE BASIC SOLN: the objective function has only negative 
            signs in front of the coefficients
            (the feasiblity is maintained in the nonnegative constraints and 
            its a basic solution).
        Example: z = 28 - x3/6 - x5/6 - 2x6/3
           then c = (c3, c5, c6)^T = ( -1/6  -1/6 -2/3 )^T
           for optimal
        
        lines 2-11 of the SIMPLEX method of Cormen et al. Chap 29.:
         
        while some index j in nIndices has c_j > 0
            do choose an index eIdx in nIndices for which c_eIdx > 0
            for each index i in bIndices
                if a[i][eIdx] > 0
                    delta[i] = b[i]/a[i][eIdx]
                else 
                    delta[i] = inf.
            choose index lIdx in bIndices that minimizes delta[i]
            if delta[lIdx] = inf {
                return "unbounded"
            }
            slackForm2 = pivot(N, B, A, b, c, v, lIdx, eIdx);
        
        */
        
        TIntSet positiveCIndexes = findPositiveCIndexes(slackFormAux);
        
        //System.out.printf("posIndexes of C=%s\n", Arrays.toString(positiveCIndexes.toArray()));
        
        double[] delta = new double[slackFormAux.bIndices.length];
        int eIdx = chooseEnteringIndex(slackFormAux, positiveCIndexes);
        
        double minDelta; // lIdx is chosen with this
        
        //choose index lIdx in bIndices that minimizes delta[i]
        while (eIdx > -1) {
                        
            minDelta = Double.POSITIVE_INFINITY; // lIdx is chosen with this
            
            for (i = 0; i < delta.length; ++i) {
                //bIdx = slackForm2.bIndices[i]; used with x
                if (slackFormAux.a[i][eIdx] > 0) {
                    delta[i] = slackFormAux.b[i]/slackFormAux.a[i][eIdx];
                    if (delta[i] < minDelta) {
                        minDelta = delta[i];
                        lIdx = i;
                    }
                } else {
                    delta[i] = Double.POSITIVE_INFINITY;
                }
            }
            
            //System.out.printf("init: pivot eIdx=%d (=x%d), lIdx=%d (=x%d)\n",
            //    eIdx, slackFormAux.nIndices[eIdx], lIdx, slackFormAux.bIndices[lIdx]);
            
            if (Double.isInfinite(minDelta)) {
                slackFormAux.state = STATE.UNBOUNDED;
                return slackFormAux;
            }
            
            slackFormAux = pivot(slackFormAux, lIdx, eIdx);
            
            //System.out.printf("init: after pivot=\n%s\n", slackFormAux.toString());
                        
            positiveCIndexes = findPositiveCIndexes(slackFormAux);
            
            eIdx = chooseEnteringIndex(slackFormAux, positiveCIndexes);
            
            //System.out.printf("positiveIndexes=\n%s\n", positiveCIndexes.toString());
        }
        
        xBasicSoln = slackFormAux.computeBasicSolution();
        if (Math.abs(xBasicSoln[0]) > eps) {
            slackFormAux.state = STATE.UNFEASIBLE;
            return slackFormAux;
        }
                
        slackFormAux.state = STATE.OPTIMAL;
        
        // remove x0 and restore original objective function
        SlackForm optimal = truncateAuxiliarySlackForm(slackFormAux,
            standForm);
        double[] x = optimal.computeBasicSolution();
        //System.out.printf("init: after truncation\n%s\n", optimal.toString());
        
        //System.out.printf("init: optimal isFeasible=%b\n%s\n", 
        //    optimal.isFeasible(), optimal.toString());
        
        optimal.state = STATE.OPTIMAL;
        
        //System.out.printf("init: slack form =\n%s\n", slackFormAux.toString());
        
        return optimal;
    }

    /**
        given a slack form and the indices for a leaving and entering variable,
        the pivot returns another slack form.
     Its the  geometrical operation of moving from a basic feasible solution to 
     * an adjacent basic feasible solution.
     * https://en.m.wikipedia.org/wiki/Simplex_algorithm
     * <pre>
     * The method is implemented from pseudocode in Section 29.3of Cormen et al.
     * </pre>
     * @param slackForm
     * @param lIdx the leaving index of x. the "leaving variable" is the basic 
     * variable in the constraint, on the lhs (so it's an index in bIndices).
     * @param eIdx the exiting index of x.  the "entering variable" is the 
     * non-basic variable in the constraint, on the rhs (so it's an index in nIndices).
     * @return 
     */
    protected SlackForm pivot(SlackForm slackForm, int lIdx /*m*/, int eIdx /*n*/) {
                        
        int m = slackForm.bIndices.length;
        int n = slackForm.nIndices.length;
                        
        double[] bHat = new double[slackForm.b.length]; /*m*/
        
        // bIndices.length = |B| = m and nIndices.length = |N| = n
        
        // NOTE: for the entering indexes in bHat and aHat, 
        // need to use the locations of the leaving indexes to store the
        // new information.  
        // bHat[eIdx] needs to be stored in bHat[lIdx] as eIdx doesn't exist in bHat yet.
        // Similarly, aHat[eIdx][*] is stored in aHat[lIdx][*]
        // and aHat[*][lIdx] is stored in aHat[*][eIdx].
        //the remaining entries continue to have same indexes.

        double eVarCoeff = slackForm.a[lIdx][eIdx];
        
        //compute the coefficients of the equation of the new basic variable x_entering.
        bHat[lIdx] = slackForm.b[lIdx]/eVarCoeff;
        
        double[][] aHat = MatrixUtil.zeros(slackForm.a.length, slackForm.a[0].length);
        
        int j;
        
        for (j = 0; j < slackForm.nIndices.length; ++j) {
            if (j == eIdx) {
                continue;
            }
            aHat[lIdx][j] = slackForm.a[lIdx][j]/eVarCoeff;
        }
        aHat[lIdx][eIdx] = 1./eVarCoeff; 
        
        // compute the coefficients of the remaining constraints
        int i;
        
        for (i = 0; i < slackForm.bIndices.length; ++i) {
            if (i == lIdx) {
                continue;
            }            
            bHat[i] = slackForm.b[i] - slackForm.a[i][eIdx] * bHat[lIdx];
            for (j = 0; j < slackForm.nIndices.length; ++j) {
                if (j == eIdx) {
                    continue;
                }
                aHat[i][j] = slackForm.a[i][j] - 
                    (slackForm.a[i][eIdx]*(slackForm.a[lIdx][j]/slackForm.a[lIdx][eIdx]));
            }
            aHat[i][eIdx] = - slackForm.a[i][eIdx]/slackForm.a[lIdx][eIdx];
        }
       
        //compute the objective function
        double vHat = slackForm.v + slackForm.c[eIdx] * bHat[lIdx];
        double[] cHat = new double[slackForm.c.length]; /*n*/
        for (j = 0; j < slackForm.nIndices.length; ++j) {
            if (j == eIdx) {
                continue;
            }
            cHat[j] = slackForm.c[j] - 
                (slackForm.c[eIdx]*( slackForm.a[lIdx][j]/slackForm.a[lIdx][eIdx]));
        }
        cHat[eIdx] = - slackForm.c[eIdx]/slackForm.a[lIdx][eIdx]; 
              
        //compute new sets of basic and nonbasic variables
        int[] nHatIndices = new int[slackForm.nIndices.length];
        for (j = 0; j < slackForm.nIndices.length; ++j) {
            if (j == eIdx) {
                nHatIndices[j] = slackForm.bIndices[lIdx];
            } else {
                nHatIndices[j] = slackForm.nIndices[j];
            }
        }
        int[] bHatIndices = new int[slackForm.bIndices.length];
        for (i = 0; i < slackForm.bIndices.length; ++i) {
            if (i == lIdx) {
                bHatIndices[i] = slackForm.nIndices[eIdx];
            } else {
                bHatIndices[i] = slackForm.bIndices[i];
            }            
        }
                
        //re-order RHS by sorted nHatIndices
        sortRHS(nHatIndices, cHat, aHat);
        
        //re-order LHS by sorted bHatIndices
        sortLHS(bHatIndices, bHat, aHat);
        
        SlackForm out = new SlackForm(nHatIndices, bHatIndices, aHat, bHat, cHat, vHat);
                
        return out;
    }

    protected void sortRHS(int[] nHatIndices, double[] cHat, double[][] aHat) {
        
        int[] idxs = new int[nHatIndices.length];
        int i;
        for (i = 1; i < nHatIndices.length; ++i) {
            idxs[i] = i;
        }
        MiscSorter.sortBy1stArg(nHatIndices, idxs);
        
        double[] c = Arrays.copyOf(cHat, cHat.length);
        double[][] a = MatrixUtil.copy(aHat);
        
        int idx, row;
        for (i = 0; i < nHatIndices.length; ++i) {
            idx = idxs[i];
            cHat[i] = c[idx];
            // columns of a get re-ordered
            for (row = 0; row < aHat.length; ++row) {
                aHat[row][i] = a[row][idx];
            }
        }
    }

    protected void sortLHS(int[] bHatIndices, double[] bHat, double[][] aHat) {
        int[] idxs = new int[bHatIndices.length];
        int i;
        for (i = 1; i < bHatIndices.length; ++i) {
            idxs[i] = i;
        }
        MiscSorter.sortBy1stArg(bHatIndices, idxs);
        
        double[] b = Arrays.copyOf(bHat, bHat.length);
        double[][] a = MatrixUtil.copy(aHat);
        
        int idx;
        for (i = 0; i < bHatIndices.length; ++i) {
            idx = idxs[i];
            bHat[i] = b[idx];
            System.arraycopy(a[idx], 0, aHat[i], 0, a[idx].length);
        }        
    }

    /**
     * finds the smallest value in b and returns its index;  if b.length is 0
     * or if b contains only Double.POSITIVE_INFINITY, -1 is returned.
     * @param b an array of double numbers.
     * @return returns the index of the smallest value in b, but if b.length is 0
     * or if b contains only Double.POSITIVE_INFINITY, -1 is returned.
     */
    protected int findMinIndex(double[] b) {
        int minIdx = -1;
        double min = Double.POSITIVE_INFINITY;
        int i;
        for (i = 0; i < b.length; ++i) {
            if (b[i] < min) {
                min = b[i];
                minIdx = i;
            }
        }
        return minIdx;
    }
    
    /**
     * Given a linear program L whose objective is minimization or maximization.
     * and which has constraints that are .leq., .eq., or .geq. the 
     * constants in b, convert L into standard form.
     * Standard form has a maximization objective, nonnegative constraints
     * on all x_j's, and constraints which are all .leq. b_i's.
     * The program implements material from Section 29.1 of Cormen et al.
     * "Introduction to Algorithms"
     <pre>
       Example Linear Program for minimization from Cormen et al., Chap 29.
             minimize:
               -2*x1 + 3*x2
             subject to constraints:
                  x1 +   x2  .eq. 7
                  x1 - 2*x2 .leq. 4
             nonnegativity constraints:
                  x1        .geq. 0
       Converted to a Standard Form:
              maximize:
                2*x1 - 3*x2 + 3*x3
              subject to constraints: 
                  x1 +   x2 - x3 .leq. 7
                 -x1 -   x2 + x3 .leq. -7
                  x1 - 2*x2 + 2*x3 .leq. 4
              nonnegativity constraints:
                  x1, x2, and x3 .geq. 0
     </pre>
     * @param a constraint coefficients
     * @param b the constants in each constraint.
     * @param c objective coefficients
     * @param constraintComparisons an index of size c.length containing indicators
     * for whether the constraint is .leq., .eq., or .geq..
     * -1 is used for .leq., 0 for .eq. and +1 for .geq.
     * <pre>
     * e.g. constraints:
     *    x1 +   x2  .eq. 7
     *    x1 - 2*x2 .leq. 4
     * would have constraintComparisons = [0, -1]
     * </pre>
     * @param isMaximization true if the linear program goal is to maximize the objective,
     * else false if the goal is to minimize the objective.
     * @param nonnegativityConstraints indicates whether x_j has a non-negativity
     * constraint
     * @return 
     */
    public static StandardForm convertLinearProgramToStandardForm(
        boolean isMaximization,
        double[][] a, double[] b, double[] c,
        int[] constraintComparisons, boolean[] nonnegativityConstraints) {
                
        int m = a.length;
        int n = c.length;
        if (constraintComparisons.length != m) {
            throw new IllegalArgumentException("constraintComparisons.length "
            + " must be equal to a.length");
        }
        if (a[0].length != n) {
            throw new IllegalArgumentException("c.length "
            + " must be equal to a[0].length");
        }
        if (a.length < 1) {
            throw new IllegalArgumentException("a.length must be > 0 ");
        }
        if (c.length < 1) {
            throw new IllegalArgumentException("c.length must be > 0 ");
        }
        if (nonnegativityConstraints.length != n) {
            throw new IllegalArgumentException("nonnegativityConstraints.length "
            + " must be equal to c.length");
        }
                
        // handle the minimization to maximization
        c = Arrays.copyOf(c, c.length);
        if (!isMaximization) {
            MatrixUtil.multiply(c, -1.);
        }
        
        // rewrite datastructures as arrays for now to make book-keeping easier 
        //    at the expense of potentially inefficient array expansion.
        // TODO: loop through nonnegativityConstraints to count the 
        //       number of needed columns and loop through constraintComparisons
        //       to count the number of needed rows,
        //       then replace these with arrays.
        List<TDoubleList> a2 = copy(a);
        TDoubleList b2 = copy(b);
        TDoubleList c2 = copy(c);
        TIntList ac2 = copy(constraintComparisons);
        
        // handle the each variable missing a non-negative constraint
        //    by replacing it with 2 non-negative variables
        // === NOTE this section only expands along j, no new rows added to a ====
        int j, i;
        int nn = 0;// number of missing non-negative constraints
        for (j = 0; j < nonnegativityConstraints.length; ++j) {
            if (!nonnegativityConstraints[j]) {
                //replace the x variable by 2 non-negative variables
                
                //handle c by appending a new variable at the end of the constraint with - coefficient of current
                c2.add(-c[j]);
                
                //handle the expansion in each constraint (row of a)
                for (i = 0; i < a.length; ++i) {
                    if (Math.abs(a[i][j]) > eps) {// x_j is not 0 for this constraint, so replace it with 2
                        a2.get(i).add(-a[i][j]);
                    } else {
                        // add a 0 for the new variable
                        a2.get(i).add(0);
                    }
                }
                nn++;
            }
        }
        assert(a2.size() == a.length);
        assert(c2.size() == (c.length + nn));
        assert(ac2.size() == constraintComparisons.length);
        assert(b2.size() == b.length);
        assert(a2.get(0).size() == a[0].length + nn);
           
        //System.out.printf("a0=\n%s\n", FormatArray.toString(a,  "%.3f"));
        //System.out.printf("a2=\n%s\n", print(a2));
        //System.out.printf("b0=\n%s\n", FormatArray.toString(b,  "%.3f"));

        // handle the constraints: .eq. and .geq.
        int m2 = 0;
        TDoubleList aRow, aRow2;
        for (i = 0; i < a.length; ++i) {
            switch (constraintComparisons[i]) {
                case -1: {
                    // an .leq. constraint is in standard form
                    break;
                }
                case 1: {
                    //convert .geq. to .leq. by multiplying the row of a by -1 and b[i] by -1
                    aRow = a2.get(i);
                    for (j = 0; j < aRow.size(); ++j) {
                        aRow.replace(j, -1.*aRow.get(j));
                    }
                    b2.replace(i, -1.*b2.get(i));
                    break;
                }
                default: {
                    //convert equalities to 2 .leq. constraints
                    // the first conversion uses the same coefficients in a2[i] and b[i] so no need to change.
                    // the second conversion uses -1 times the coefficients in a2[i] and b[i] and gets inserted at end of both
                    aRow = a2.get(i);
                    aRow2 = new TDoubleArrayList(aRow.size());
                    for (j = 0; j < aRow.size(); ++j) {
                        aRow2.add(-1.*aRow.get(j));
                    }
                    a2.add(aRow2);
                    b2.add(-1.*b2.get(i));
                    break;
                }
            }
        }
        
        //System.out.printf("a2=\n%s\n", print(a2));
        //System.out.printf("b=\n%s\n", FormatArray.toString(b2.toArray(),  "%.3f"));

        double[][] a3 = copy(a2);
         
        StandardForm standForm = new StandardForm(a3, b2.toArray(), c2.toArray(), 0);
        return standForm;
    }
   
    /**
     * given a standard form containing .leq. inequalities expressed by 
     * standForm.a and standForm.b, convert the problem to a SlackForm
     * of equality constraints and their slack variables.
     * NOTE that the slack form 'a' matrix will have sign conventions
     * that are the same as the signs of standForm.a.
     * The slack form constraints are slack variable = b_i - summation_j=1:n(a_i_j * x_j),
     * and so the negative sign read in the written slack form is not present
     * in the matrix 'a'.
     * @param standForm a linear program with an objective of maximization,
       subject to constraints that are inequality constraints of the form .leq.
       and non-negativity constraints on x.
       
     * @return 
     */
    public static SlackForm convertConstraints(StandardForm standForm) {
        /*
          Standard Form:
              real numbers: c1...cn; b1,...bm; and aij for a=1:m and j=1:n.

                            Find numbers x1,...xn:
            objective ->    maximize summation_j=1:n(cj*xj)
            constraints ->   subject to: summation_j=1:n(aij*xj) .leq. bi for i=1:m
            non-negativity constraints ->   xj .geq. 0 for j=1:n

            OR expressed more compactly:
              A = (aij) =  mXn matrix
              b = (bi) an m-dimensional vector
              c = (cj) an n-dimensional vector
              x = (xj) an n-dimensional vector
            objective ->    maximize c^T*x
            constraints ->   subject to: A*x .leq. b
                             x .leq. 0
        */
        int n = standForm.c.length;
        int m = standForm.b.length;
        
        double[][] aHat = MatrixUtil.copy(standForm.a);
        assert(aHat.length == m);
        assert(aHat[0].length == n);
        
        double[] bHat = Arrays.copyOf(standForm.b, standForm.b.length);
        double[] cHat = Arrays.copyOf(standForm.c, standForm.c.length);
        assert(cHat.length == n);
        
        // removing redundant constraints
        int i, i2, j;
        double[] aRow, aRow2;
        boolean same;
        for (i = aHat.length - 1; i >= 0; i--) {
            aRow = aHat[i];
            for (i2 = i - 1; i2 >= 0; i2--) {
                aRow2 = aHat[i2];
                // if b's are same, ignoring signs, compare the rest of the constraint
                if (Math.abs(bHat[i] - bHat[i2]) < eps) {
                    same = true;
                    for (j = 0; j < aHat[i].length; ++j) {
                        if (Math.abs(aRow[j] - aRow2[j]) > eps) {
                            same = false;
                            break;
                        }
                    }
                } else if (Math.abs(bHat[i] + bHat[i2]) < eps) {
                    same = true;
                    for (j = 0; j < aHat[i].length; ++j) {
                        if (Math.abs(aRow[j] + aRow2[j]) > eps) {
                            same = false;
                            break;
                        }
                    }
                } else {
                    same = false;
                }
                if (same) {
                    bHat = removeElement(bHat, i);
                    aHat = removeRow(aHat, i);
                    break;
                }
            } 
        }
        
        m = bHat.length;
                   
        // writing in 0-based indexes, the x index of non-basic and basic variables
        int[] nHatIndices = new int[n]; 
        int[] bHatIndices = new int[m];
        
        for (i = 0; i < n; ++i) {
            nHatIndices[i] = i;
        }
        for (i = 0; i < m; ++i) {
            bHatIndices[i] = n + i;
        }
        double vHat = standForm.v;
        
        SlackForm slackForm = new SlackForm(nHatIndices, bHatIndices, aHat, bHat, cHat, vHat);
        
        return slackForm;
    }

    /*
        L is in standard Form.
        L_aux adds an artificial variable x0, following eqns (29.109) - (29.111)
    <pre>    
        maximize: -x0  (which is the same as minimize x0)
        subject to: summation_j=1:n(a_i_j*x_j) - x0 <= b_i for i=1:m
                    x_j >= 0 for j=1:n  
    
        If all of b_i are non-negative, then the minimum x0 for a feasible solution is x0=0.
        If any of b_i are negative (and all of the constraint's 'a' coefficients are
        non-negative), then the minimum x0 for a feasible solution is the smallest b_i:
            e.g. 
             a11*x1 + a12*x2 - x0 <= -10
             a21*x1          - x0 <= -5
                      a22*x2 - x0 <= 2
              if x1,x2=0, then x0=10 which is non-negative
              but if any of the coefficients of 'a' are negative in the
                constraint with a negative b, then the non-basic variables in
                that constraint (the other x's) can be > 0 and x0 can be 0.
    
    Sect 29.5 of Cormen et al., method INITIALIZE-SIMPLEX(A, b,c)
    line 4:   form L_AUX by adding -x0 to the LHS of each equation and setting the
              objective function to -x0.
    </pre>    
    */
    protected SlackForm createAuxiliarySlackForm(StandardForm standForm) {
        
        //System.out.printf("createAuxiliarySlackForm\n");
        
        SlackForm slackForm = convertConstraints(standForm);
        
        //System.out.printf("SlackForm=%S\n", slackForm.toString());
        
        int m = slackForm.b.length;
        int n = slackForm.c.length;
                
        int i, j;
        
        /*
        double minBForPosA = Double.POSITIVE_INFINITY;
        
        for (i = 0; i < slackForm.b.length; ++i) {
            if (slackForm.b[i] < 0) {
                boolean allPosA = true;
                for (j = 0; j < slackForm.a[i].length; ++j) {
                    if (slackForm.a[i][j] < 0) {
                        allPosA = false;
                        break;
                    }
                }
                if (allPosA) {
                    if (slackForm.b[i] < minBForPosA) {
                        minBForPosA = slackForm.b[i];
                    }
                }
            }
        }
        
        double x0 = 0;
        if (minBForPosA < 0) {
            x0 = -minBForPosA;
        }
        */
                                
        //also see end of Section 5.6, pg 57 and pg 58 of Matousek & Gartner "Undegstanding and Using Linear Programming"
        
        //double[] xHat = new double[m + n + 1];
        //System.arraycopy(xBasicSoln, 0, xHat, 1, m + n);
               
        int[] nHatIndices = new int[n + 1];
        for (i = 1; i < nHatIndices.length; ++i) {
            nHatIndices[i] = slackForm.nIndices[i - 1] + 1;
        }
        int[] bHatIndices = Arrays.copyOf(slackForm.bIndices, m);
        for (i = 0; i < m; ++i) {
            bHatIndices[i]++;
        }
        
        double[] cHat = new double[n + 1];
        cHat[0] = -1;
                
        double[] bHat = Arrays.copyOf(slackForm.b, m);
        
        double[][] aHat = new double[m][];//mX(n+1)
        for (i = 0; i < m; ++i) {
            aHat[i] = new double[n+1];
            System.arraycopy(slackForm.a[i], 0, aHat[i], 1, n);
            aHat[i][0] = -1;
        }
        
        //double[] xBasicSoln = slackForm.computeBasicSolution();
        //double eval = slackForm.evaluateObjective();
        double vHat = 0;//eval;
        
        SlackForm slackForm2 = new SlackForm(nHatIndices, bHatIndices, 
            aHat, bHat, cHat, vHat);
         
        return slackForm2;
    }
    
    // remove x0 and restore original objective function
    private SlackForm truncateAuxiliarySlackForm(SlackForm slackFormAux,
        StandardForm origForm) {
        
        int m = slackFormAux.bIndices.length;
        int n = slackFormAux.nIndices.length;
        
        int i, i2 = 0;
        int[] nHatIndices = new int[slackFormAux.nIndices.length - 1];
        for (i = 0; i < n; ++i) {
            if (slackFormAux.nIndices[i] != 0) {
                nHatIndices[i2] = slackFormAux.nIndices[i] - 1;
                i2++;
            }
        }
        int[] bHatIndices = new int[slackFormAux.bIndices.length];
        for (i = 0; i < m; ++i) {
            bHatIndices[i] = slackFormAux.bIndices[i] - 1;
        }
        
        double[] bHat = Arrays.copyOf(slackFormAux.b, slackFormAux.b.length);
        
        double[][] aHat = new double[m][];//mX(n-1)
        int j;
        for (i = 0; i < m; ++i) {
            aHat[i] = new double[n-1];
            i2 = 0;
            for (j = 0; j < n; ++j) {
                if (slackFormAux.nIndices[j] != 0) {
                    aHat[i][i2] = slackFormAux.a[i][j];
                    i2++;
                }
            }
        }
        
        /*
        line 11: create a slack form for L  for which basic soln is feasible.
        to do so, delete x0 terms from constraints and restore the orig
        obj func for L.
        the orig obj func may contain both basic and nonbasic vars.
        therefore in the obj func we replace each basic var by the RHS of its assoc constraint.
        
        e.g.
            L_aux
                x2 =  4/5 -  x0/5 +  x1/5  +  x4/5 
                x3 = 14/5 + 4x0/5 - 9x1/5  +  x4/5 
            truncated L_aux
                x2 =  4/5 +  x1/5  +  x4/5   NOTE: x2 is x[1[ and x3 is x[2] etc for 0-based indexes
                x3 = 14/5 - 9x1/5  +  x4/5
                nIndices=[0, 3] <-- 0-based indexes
                bIndices=[1, 2] <-- 0-based indexes
        
            orig obj func = 2x1 - x2
                          = 2x1 - (4/5 + x1/5 + x4/5)
                          = -4/5 + 9x1/5 -x4/5
        */
        TIntSet nHatIdxs = new TIntHashSet(nHatIndices);
        // orig c indexes are 0 thru origC.length-1
        TIntIntMap bHatIdxMap = new TIntIntHashMap();
        int idx;
        for (i = 0; i < bHatIndices.length; ++i) {
            idx = bHatIndices[i];
            bHatIdxMap.put(idx, i);
        }
        
        TIntIntMap nHatIdxMap = new TIntIntHashMap();
        for (i = 0; i < nHatIndices.length; ++i) {
            idx = nHatIndices[i];
            nHatIdxMap.put(idx, i);
        }
        
        double cCoeff;
        int jIdx, jj;
        double vHat = origForm.v;
        double tmp;
        double[] cHat = new double[origForm.c.length];
        double[] aHatRow;
        double bHatI;
        for (i = 0; i < origForm.c.length; ++i) {
            if (nHatIdxs.contains(i)) { //i=1
                cHat[i] = origForm.c[i];
            }
        }
        // any orig c indexes not in nHatIdxs need to be rewritten using the basic constraint
        for (i = 0; i < origForm.c.length; ++i) {
            if (!nHatIdxs.contains(i)) {
                cCoeff = origForm.c[i];
                //find the row for index in bindices
                assert(bHatIdxMap.containsKey(i));
                j = bHatIdxMap.get(i);
                
                aHatRow = aHat[j]; 
                bHatI = bHat[j]; 
                
                vHat += (cCoeff * bHatI);
                
                for (jj = 0; jj < aHatRow.length; ++jj) {
                    jIdx = nHatIndices[jj];
                    assert(nHatIdxMap.containsKey(jIdx));
                    jIdx = nHatIdxMap.get(jIdx);
                    tmp = cCoeff * -aHatRow[jj];
                    cHat[jIdx] += tmp;
                }
            }
        }        
        
        SlackForm slackForm = new SlackForm(nHatIndices, bHatIndices, 
            aHat, bHat, cHat, vHat);
        
        //System.out.printf("truncated=\n%s\n", slackForm.toString());
        
        return slackForm;
    }

    protected int chooseEnteringIndex(SlackForm slackForm, TIntSet positiveCIndexes) {
        
        /*
        In the selection of c_j (the entering variable), 
           can use Bland's rule w/ caveat that it is slow:
               choose the entering (nonbasic) variable to be the smallest index
               among the eligible (that is, the positive remaining c coefficients)
               qnd if there is more than one leaving variable w/ same value,
               choose the one w/ smallest index.
           or can randomly select c_j out of the positive remaining coefficients
        */
        
        // choose eIdx randomly from nonNegativeC
        int eIdx = -1;
        int nC = positiveCIndexes.size();
        switch (nC) {
            case 0:
                return eIdx;
            case 1:
                return positiveCIndexes.iterator().next();
            default:
                eIdx = rand.nextInt(nC);
                return positiveCIndexes.toArray()[eIdx];
        }
    }

    /**
     * find the indexes eligible to be entering variable indexes.  They are
     * found as the indexes of the c coefficients which are still positive.
     * @param slackForm
     * @return 
     */
    protected TIntSet findPositiveCIndexes(SlackForm slackForm) {
        TIntSet cs = new TIntHashSet();
        int n = slackForm.c.length;
        int i;
        for (i = 0; i < n; ++i) {
            if (slackForm.c[i] > 0) {
                cs.add(i);
            }
        }
        return cs;
    }
    
    private static List<TDoubleList> copy(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        List<TDoubleList> a2 = new ArrayList<TDoubleList>(m);
        int i;
        for (i = 0; i < m; ++i) {
            a2.add(new TDoubleArrayList(Arrays.copyOf(a[i], a[i].length)));
            assert(a2.get(i).size() == n);
        }
        return a2;
    }
    
    private static double[][] copy(List<TDoubleList> a) {
        int m = a.size();
        int n = a.get(0).size();
        double[][] a2 = new double[m][];
        int i;
        for (i = 0; i < m; ++i) {
            a2[i] = a.get(i).toArray();
            assert(a2[i].length == n);
        }
        return a2;
    }
    
    private static double[] removeElement(double[] b, int idx) {
        double[] out = new double[b.length - 1];
        int i, i2 = 0;
        for (i = 0; i < b.length; ++i) {
            if (i != idx) {
                out[i2] = b[i];
                i2++;
            }
        }
        return out;
    }

    private static double[][] removeRow(double[][] a, int rowIdx) {
        double[][] out = new double[a.length - 1][];
        int i, i2 = 0;
        for (i = 0; i < a.length; ++i) {
            if (i != rowIdx) {
                out[i2] = Arrays.copyOf(a[i], a[i].length);
                i2++;
            }
        }
        return out;
    }
    
    private static String print(List<TDoubleList> a) {
        StringBuilder sb = new StringBuilder();
        int m = a.size();
        int n = a.get(0).size();
        int i;
        for (i = 0; i < m; ++i) {
            sb.append(FormatArray.toString(a.get(i).toArray(), "%.3f")).append("\n");
        }
        return sb.toString();
    }
        
    private static TDoubleList copy(double[] b) {
        TDoubleList b2 = new TDoubleArrayList(Arrays.copyOf(b, b.length));
        return b2;
    }
    
    private static TIntList copy(int[] cc) {
        TIntList cc2 = new TIntArrayList(Arrays.copyOf(cc, cc.length));
        return cc2;
    }
    
    /**
     * maximization of a linear function subject to linear *inequalities*
     */
    public static class StandardForm extends FormTuple {
    
        public StandardForm(double[][] a, double[] b, double[] c, double v) {
            super(a, b, c, v);
        }
        
        /**
         * when x has been calculated, evaluate the objective.
         * see eqn (29.42) of Cormen et al.
         * <pre> z = v + summation_j_in_nIndices(c_j*x_j) </pre>
         * @return 
         */
        @Override
        public double evaluateObjective() {
            if (x == null) {
                throw new IllegalStateException("x has not been calculated to evaluate");
            }
            double sum = v;
            int i;
            for (i = 0; i < c.length; ++i) {
                // in standard form, x is only the non-basic variables, no basic variables,
                //   so can use the index i.
                sum += c[i]*x[i];
            }
            return sum;
        }
    
        public boolean isFeasible() {
            /*
            satifies constraints
                A*x .leq. b = summation_j=1:n(aij*xj) .leq. bi for i=1:m
                xj .geq. 0 for j=1:n
                  the later is a nonnegativity constraint
            */
            if (x == null) {
                throw new IllegalStateException("x cannot be null.  calculate a basic solution first");
            }
            int m = b.length;
            int n = c.length;
            int i, j;
            double sum;
            for (i = 0; i < m; ++i) {
                sum = 0;
                for (j = 0; j < n; ++j) {
                    // in standard form, x is only the non-basic variables, no basic variables,
                //   so can use the index i.
                    sum += (a[i][j]*x[i]);
                }
                if (sum > b[i]) {
                    return false;
                }
            }
            for (i = 0; i < x.length; ++i) {
                if (x[i] < 0) {
                    return false;
                }
            }
            return true;
        }
        
        /**
         * using the present forms, sets all rhs x's to 0 to solve for the lhs x's.
         * The x vector is returned.
         * @return 
         */
        @Override
        public double[] computeBasicSolution() {
            //following Sect 29.3 of Cormen et al. "Introduction to Algorithms"
            double[] xt = new double[a.length + b.length];
            int i;
            for (i = 0; i < b.length; ++i) {
                //x_i = b_i - summation_j_in_nIndices(a_i_j*x_j) for i in bIndices
                // setting rhs x to 0
                // in standard form, x is only the non-basic variables, no basic variables,
                //   so can use the index i.
                xt[i] = b[i];
            }
            this.x = xt;
            return Arrays.copyOf(x, x.length);
        }
        
    }
    
    /**
     * maximization of a linear function subject to linear *equalities*
     */
    public static class SlackForm extends FormTuple {
        
        /**
         * denotes the indices of the nonbasic variables (rhs vars).
         * nIndices.length = n.
         */
        public int[] nIndices;
        
        /**
         * denotes the indices of the basic variables (lhs vars).
         * bIndices.length = m.
         */
        public final int[] bIndices;
        
        public static enum STATE {
            UNBOUNDED, UNFEASIBLE, FEASIBLE, OPTIMAL;
        }
        
        protected STATE state;
        
        /**
         * 
         * @param nIndices denotes the indices of the nonbasic variables (rhs vars).
         * nIndices.length = n.
         * @param bIndices denotes the indices of the basic variables (lhs vars).
         * bIndices.length = m.
         * @param a
         * @param b
         * @param c
         * @param v 
         */
        public SlackForm(int[] nIndices, int[] bIndices, double[][] a, 
            double[] b, double[] c, double v) {
            super(a, b, c, v);
            
            int n = c.length;
            int m = b.length;
            if (nIndices.length != n) {
                throw new IllegalArgumentException("nIndices.length must equal c.length");
            }
            if (bIndices.length != m) {
                throw new IllegalArgumentException("bIndices.length must equal b.length");
            }
            this.nIndices = Arrays.copyOf(nIndices, nIndices.length);
            this.bIndices = Arrays.copyOf(bIndices, bIndices.length);
        }
        
        /**
         * 
         * @param nIndices denotes the indices of the nonbasic variables (rhs vars).
         * nIndices.length = n.
         * @param bIndices denotes the indices of the basic variables (lhs vars).
         * bIndices.length = m.
         * @param a
         * @param b
         * @param c
         * @param v
         * @param state 
         */
        public SlackForm(int[] nIndices, int[] bIndices, double[][] a, 
            double[] b, double[] c, double v, STATE state) {
            this(nIndices, bIndices, a, b, c, v);
            this.state = state;
        }
    
        public boolean isFeasible() {
            /*
            satifies constraints
                A*x .leq. b = summation_j=1:n(aij*xj) .leq. bi for i=1:m
                xj .geq. 0 for j=1:n
                  the later is a nonnegativity constraint
            */
            if (x == null) {
                throw new IllegalStateException("x cannot be null.  calculate a basic solution first");
            }
            int m = b.length;
            int n = c.length;
            int i, j, idxJ;
            double sum;
            for (i = 0; i < m; ++i) {
                sum = 0;
                for (j = 0; j < n; ++j) {
                    idxJ = nIndices[j];
                    sum += (a[i][j]*x[idxJ]);
                }
                if (sum > b[i]) {
                    return false;
                }
            }
            for (i = 0; i < x.length; ++i) {
                if (x[i] < 0) {
                    return false;
                }
            }
            return true;
        }
        
        /**
         * when x has been calculated, evaluate the objective.
         * see eqn (29.42) of Cormen et al.
         * <pre> z = v + summation_j_in_nIndices(c_j*x_j) </pre>
         * @return 
         */
        @Override
        public double evaluateObjective() {
            if (x == null) {
                throw new IllegalStateException("x has not been calculated to evaluate");
            }
            double sum = v;
            int i, idx;
            for (i = 0; i < nIndices.length; ++i) {
                idx = nIndices[i];
                sum += c[i]*x[idx];
            }
            return sum;
        }
        
        /**
         * using the present forms, sets all rhs x's to 0 to solve for the lhs x's.
         * The x vector is returned.
         * @return 
         */
        @Override
        public double[] computeBasicSolution() {
            /*
            following Sect 29.3 of Cormen et al. "Introduction to Algorithms"

            Slack Form example:
                 maximize 3*x1 + x2 + 2*x3
                 subject to:
                      x4 = 30 -   x1 -   x2 - 3*x3
                      x5 = 24 - 2*x1 - 2*x2 - 5*x3
                      x6 = 36 - 4*x1 -   x2 - 2*x3
                      x1, x2, x3, x4, x5, x6 .geq. 0
            */  
            double[] xt = new double[nIndices.length + bIndices.length];
            int ii, i;
            for (ii = 0; ii < bIndices.length; ++ii) {
                i = bIndices[ii];
                //x_i = b_i - summation_j_in_nIndices(a_i_j*x_j) for i in bIndices
                // setting rhs x to 0
                xt[i] = b[ii];
            }
            this.x = xt;
            return Arrays.copyOf(x, x.length);
        }
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder(super.toString());
            sb.append("nIndices=").append(Arrays.toString(nIndices)).append("\n");
            sb.append("bIndices=").append(Arrays.toString(bIndices)).append("\n");
            if (state != null) {
                sb.append("STATE=").append(state.name()).append("\n");
            }
            return sb.toString();
        }
    
    }
    
    /**
        example Standard Form:
        
        Find numbers x1,...xn:
        objective ->    maximize 
                           c^T*x = summation_j=1:n(cj*xj)
        constraints ->  subject to: 
                           A*x .leq. b = summation_j=1:n(aij*xj) .leq. bi for i=1:m
        constraints ->   xj .geq. 0 for j=1:n
                         the later is a nonnegativity constraint
        
        where c1...cn; b1,...bm; and aij are real numbers.
    */
    public abstract static class FormTuple {
        
        /**
         * mXn matrix
         */
        public double[][] a;
        
        /**
         * an m-dimensional vector
         */
        public double[] b;
        
        /**
         * an n-dimensional vector
         */
        public double[] c;
       
        /**
         * an n-dimensional vector
         */
        public double[] x = null;
        
        /**
         * an optional term v is sometimes present in the objective
         */
        public double v = 0;
        
        /**         
         * @param a mXn matrix of constraint coefficients.   Careful with the signs
         * of the values for the SlackForm. see eqn (29.43) of Cormen et al.  
         * <pre> x_i = b_i - summation_j(a_i_j*x_j) for i in bIndices. </pre>
         * @param b m-dimensional vector of constraint inequalities
         * @param c n-dimensional vector of objective coefficients.
         * @param v an optional term sometimes present in the objective.  can be 0 if no v is used.
         */
        public FormTuple(double[][] a, double[] b, double[] c, double v) {
            
            int n = c.length;
            int m = b.length;
            if (a.length != m || a[0].length != n) {
                throw new IllegalArgumentException("a must be of dimensions [b.length][c.length]");
            }
            this.a = MatrixUtil.copy(a);
            this.b = Arrays.copyOf(b, b.length);
            this.c = Arrays.copyOf(c, c.length);
            this.v = v;
            
        }
        
        /**
         * convert the c coefficients from minimization coefficients to
         * maximization coefficients by multiplying them by -1.
         */
        public void convertMinimizationToMaximiation() {
            MatrixUtil.multiply(c, -1.);
        }
        
        public abstract double evaluateObjective();
        
        /**
         * using the present forms, sets all rhs x's to 0 to solve for the lhs x's.
         * The x vector is returned.
         * @return 
         */
        public abstract double[] computeBasicSolution();
        
        
        @Override
        public String toString() {

            StringBuilder sb = new StringBuilder();
            
            sb.append(String.format("v=%.3f\n", v));
            sb.append("c=").append(FormatArray.toString(c, "%.3f")).append("\n");
            
            sb.append("b=").append(FormatArray.toString(b, "%.3f")).append("\n");
            sb.append(String.format("a=\n%s", FormatArray.toString(a, "%.3f")));
            if (x != null) {
                sb.append("x=").append(FormatArray.toString(x, "%.3f")).append("\n");
            }
            return sb.toString();
        }
    }
}
