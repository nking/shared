package algorithms.optimization;

import algorithms.matrix.MatrixUtil;
import algorithms.sort.MiscSorter;
import algorithms.util.FormatArray;
import java.util.Arrays;

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
        
Simplext Method:
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
 * @author nichole
 */
public class LinearProgramming {
    
    /**
        given a linear program in Standard Form if the problem is unfeasible, 
        this method returns null, else it
        returns a slack form for which a basic solution is feasible.
     * <pre>
     * The method is implemented from pseudocode in Section 29.3 of Cormen et al.
     * </pre>
     * @param standForm
     * @return 
     */
    protected SlackForm initialzeSimplex(StandardForm standForm) {
        
        int m = standForm.bIndices.length;
        int n = standForm.nIndices.length;
        
        int i;
        int lIdx = findMinIndex(standForm.b);
        
        //is the inital basic soltuion feasible?
        if (lIdx > -1 && standForm.b[lIdx] >= 0) {
            int[] nIndices = new int[n];
            int[] bIndices = new int[m];
            //nIndices = [1,2,...n]. NOTE: using 0-based indexes instead
            for (i = 1; i < n; ++i) {
                nIndices[i] = i;
            }
            //bIndices = [n+1,n+2,...n+m]. NOTE: using 0-based indexes instead
            for (i = 1; i < m; ++i) {
                bIndices[i] = n+i;
            }
            
            SlackForm slackForm = new SlackForm(nIndices, bIndices, 
                standForm.a, standForm.b, standForm.c, standForm.v);
            
            return slackForm;
        }
               
        /*
        L is in standard Form.
        L_aux is an alternation for n+1 variables", following eqns (29.109) - (29.111)
        
        maximize: -x0
        subject to: summation_j=1:n(a_i_j*x_j) - x0 <= b_i for i=1:m
                    x_j >= 0 for j=1:n     
        
        set all RHS x's to 0
        
        when z==v, the summation_j_in_nIndices(c_j*x_j) is 0
        
        */
        
        double[] xHat = new double[m + n + 1];
        System.arraycopy(xHat, 0, standForm.computeBasicSolution(), 1, m + n);
        double eval = standForm.evaluateObjective();
                
        /*
        Form L_aux by adding -x0 to the LHS of each equation
            and setting the objective function to -x0
        
        let (N, B, A, b, c, v) be the resulting Slack Form for L_aux
        
        //L_aux has n+1 nonbasic variables and m basic variables
        (N, B, A, b, c, v) = pivot((N, B, A, b, c, v, ell, 0)
        
        //the basic solution is now feasible for L_aux
        iterate the while loop of lines 7-11 of the SIMPLEX method
            until an optimal solution to L_aux is found
        
        if the basic solution sets xHat0 to 0
            then return the final slack form with x0 removed and the original
               objective function restored.
        else
            return "unfeasible"        
        */
        
        throw new UnsupportedOperationException("not yet implemented");
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
     * nonbasic variable in the constraint, on the rhs (so it's an index in nIndices).
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

        //compute the coefficients of the equation of the new basic variable x_entering.
        // bHat[eIdx] =  b[lIdx]/a[lIdx][eIdx]
        bHat[lIdx] = slackForm.b[lIdx]/slackForm.a[lIdx][eIdx];
        
        double[][] aHat = MatrixUtil.zeros(slackForm.a.length, slackForm.a[0].length);
        
        int j, jX;
        for (j = 0; j < slackForm.nIndices.length; ++j) {
            if (j == eIdx) {
                continue;
            }
            jX = slackForm.nIndices[j];
            //aHat[eIdx][j] = slackForm.a[lIdx][j]/slackForm.a[lIdx][eIdx];
            aHat[lIdx][j] = slackForm.a[lIdx][j]/slackForm.a[lIdx][eIdx];
        }
        //aHat[eIdx][lIdx] = 1./slackForm.a[lIdx][eIdx];
        aHat[lIdx][eIdx] = 1./slackForm.a[lIdx][eIdx];
        
        // compute the coefficients of the remaining constraints
        int i, iX;
        for (i = 0; i < slackForm.bIndices.length; ++i) {
            if (i == lIdx) {
                continue;
            }
            iX = slackForm.bIndices[i];
            
            //bHat[i] = slackForm.b[i] - slackForm.a[i][eIdx] * bHat[eIdx];
            bHat[i] = slackForm.b[i] - slackForm.a[i][eIdx] * bHat[lIdx];
            for (j = 0; j < slackForm.nIndices.length; ++j) {
                if (j == eIdx) {
                    continue;
                }
                jX = slackForm.nIndices[j];
                       
                //aHat[i][j] = slackForm.a[i][j] - slackForm.a[i][eIdx]*aHat[eIdx][j];
                aHat[i][j] = slackForm.a[i][j] - slackForm.a[i][eIdx]*aHat[lIdx][j];
            }
            //aHat[i][lIdx] = - slackForm.a[i][eIdx]*aHat[eIdx][lIdx];
            aHat[i][eIdx] = - slackForm.a[i][eIdx]*aHat[lIdx][eIdx];
        }
       
        //compute the objective function
        //double vHat = slackForm.v + slackForm.c[eIdx] * bHat[eIdx];
        double vHat = slackForm.v + slackForm.c[eIdx] * bHat[lIdx];
        double[] cHat = new double[slackForm.c.length]; /*n*/
        for (j = 0; j < slackForm.nIndices.length; ++j) {
            if (j == eIdx) {
                continue;
            }
            jX = slackForm.nIndices[j];
            //cHat[j] = slackForm.c[j] - slackForm.c[eIdx] * aHat[eIdx][j];
            cHat[j] = slackForm.c[j] - slackForm.c[eIdx] * aHat[lIdx][j];
        }
        //cHat[lIdx] = -slackForm.c[eIdx] * aHat[eIdx][lIdx];
        cHat[eIdx] = -slackForm.c[eIdx] * aHat[lIdx][eIdx];
                
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
            for (row = 0; row < aHat[i].length; ++row) {
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
     * maximization of a linear function subject to linear *inequalities*
     */
    public static class StandardForm extends FormTuple {
    
        public StandardForm(int[] nIndices, int[] bIndices, double[][] a, 
            double[] b, double[] c, double v) {
            super(nIndices, bIndices, a, b, c, v);
        }
    
    }
    
    /**
     * maximization of a linear function subject to linear *equalities*
     */
    public static class SlackForm extends FormTuple {
    
        public SlackForm(int[] nIndices, int[] bIndices, double[][] a, 
            double[] b, double[] c, double v) {
            super(nIndices, bIndices, a, b, c, v);
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
    public static class FormTuple {
        
        /**
         * mXn matrix
         */
        double[][] a;
        
        /**
         * an m-dimensional vector
         */
        double[] b;
        
        /**
         * an n-dimensional vector
         */
        double[] c;
        
        /**
         * an n-dimensional vector
         */
        double[] x = null;
        
        /**
         * an optional term v is sometimes present in the objective
         */
        double v = 0;
        
        /**
         * denotes the indices of the nonbasic variables (rhs vars).
         * nIndices.length = n.
         */
        final int[] nIndices;
        
        /**
         * denotes the indices of the basic variables (lhs vars).
         * bIndices.length = m.
         */
        final int[] bIndices;
        
        /**
         
         * @param nIndices denotes the indices of the nonbasic variables (rhs vars).
         * nIndices.length = n.
         * @param bIndices denotes the indices of the basic variables (lhs vars).
         * bIndices.length = m.
         * @param a mXn matrix of constraint coefficients.   Careful with the signs
         * of the values. see eqn (29.43) of Cormen et al.  
         * <pre> x_i = b_i - summation_j_in_nIndices(a_i_j*x_j) for i in bIndices. </pre>
         * @param b m-dimensional vector of constraint inequalities
         * @param c n-dimensional vector of objective coefficients.
         * @param v an optional term sometimes present in the objective.  can be 0 if no v is used.
         */
        public FormTuple(int[] nIndices, int[] bIndices, double[][] a,
            double[] b, double[] c, double v) {
            
            this.nIndices = Arrays.copyOf(nIndices, nIndices.length);
            this.bIndices = Arrays.copyOf(bIndices, bIndices.length);
            int n = nIndices.length;
            int m = bIndices.length;
            if (a.length != m || a[0].length != n) {
                throw new IllegalArgumentException("a must be of dimensions [bIndicies.length][nIndices.length]");
            }
            if (b.length != m) {
                throw new IllegalArgumentException("b must be length bIndicies.length");
            }
            if (c.length != n) {
                throw new IllegalArgumentException("c must be length nIndicies.length");
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
        
        /**
         * when x has been calculated, evaluate the objective.
         * see eqn (29.42) of Cormen et al.
         * <pre> z = v + summation_j_in_nIndices(c_j*x_j) </pre>
         * @return 
         */
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

            StringBuilder sb = new StringBuilder();
            
            sb.append(String.format("v=%.3f\n", v));
            sb.append("c=").append(FormatArray.toString(c, "%.3f")).append("\n");
            
            sb.append("b=").append(FormatArray.toString(b, "%.3f")).append("\n");
            sb.append(String.format("a=\n%s\n", FormatArray.toString(a, "%.3f")));
            sb.append("bIndices=").append(Arrays.toString(bIndices)).append("\n");
            sb.append("nIndices=").append(Arrays.toString(nIndices)).append("\n");
            if (x != null) {
                sb.append("x=").append(FormatArray.toString(x, "%.3f")).append("\n");
            }
            return sb.toString();
        }
    }
}
