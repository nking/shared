/*
Class PolynomialRootSolver is adapted from numpy.roots source
code available at
https://github.com/numpy/numpy/blob/v1.18.1/numpy/lib/polynomial.py#L168-L251
which is licensed under BSD-3 Clause "New" or "Revised" License.
https://github.com/numpy/numpy/blob/v1.18.1/LICENSE.txt

Copyright (c) 2005-2019, NumPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

    * Neither the name of the NumPy Developers nor the names of any
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package algorithms.misc;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.util.Arrays;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * calculates polynomial roots by forming "the companion matrix"
 * and using matrix eigenvalue decomposition to find the eigenvalues as the roots.
 * The code is adapted from numpy.roots() source code.
 * 
 * @author nichole
 */
public class PolynomialRootSolver {
    
    /**
     * calculates polynomial roots by forming "the companion matrix"
     * and using matrix eigenvalue decomposition to find the eigenvalues as the roots.
     * The code is adapted from numpy.roots() source code.
 
     * @param coeffs coefficients of a polynomial given in the order of decreasing 
     * exponential, e.g. expecting [4, 3, 2, 1] for 4*x^3 + 3*x^2 + 2*x + 1 = 0.
     * 
     * @return
     * @throws NotConvergedException 
     */
    public static Complex[] roots(double[] coeffs) throws NotConvergedException {
           
        // http://web.mit.edu/18.06/www/Spring17/Eigenvalue-Polynomials.pdf
        // form the companion matrix as the characteristic polynomial
        // https://en.wikipedia.org/wiki/Companion_matrix
        // Also: The Vandermonde determinant was sometimes called the discriminant, 
        // although, presently, the discriminant of a polynomial is the square 
        // of the Vandermonde determinant of the roots of the polynomial.
        
        
        double eps = 1.e-5;
        
        int[] non_zero = nonzero(coeffs, eps);
        
        //System.out.println("non_zero=" + Arrays.toString(non_zero));
        
        if (non_zero.length == 0) {
            return new Complex[]{};
        }
        
        int trailing_zeros = coeffs.length - non_zero[non_zero.length - 1] - 1;
        
        coeffs = Arrays.copyOfRange(coeffs, non_zero[0], non_zero[non_zero.length - 1] + 1);
        
        int n = coeffs.length;
        
        Complex[] roots;
                
        if (n > 1) {
            // create a matrix of size n-1 all zeros. place ones along
            // the diagonal below the main diagonal
            double[][] m = MatrixUtil.zeros(n-1, n-1);
            for (int i = 1; i < m.length; ++i) {
                m[i][i-1] = 1;
            }
            // divide -p by first element and insert all into first row of m except first entry of p
            MatrixUtil.multiply(coeffs, -1./coeffs[0]);
            System.arraycopy(coeffs, 1, m[0], 0, coeffs.length-1);
            
            //System.out.printf("A=\n%s\n", FormatArray.toString(m, "%.3e"));
            
            EVD evd = EVD.factorize(new DenseMatrix(m));
            
            double[] vR = evd.getRealEigenvalues();
            double[] vI = evd.getImaginaryEigenvalues();
            
            roots = new Complex[vR.length];
            
            for (int i = 0; i < vR.length; ++i) {
                roots[i] = new Complex(vR[i], vI[i]);
            }
            
        } else {
            roots = new Complex[0];
        }
        
        if (trailing_zeros > 0) {
            Complex[] r2 = new Complex[roots.length + trailing_zeros];
            for (int i = 0; i < roots.length; ++i) {
                r2[i] = roots[i];
            }
            for (int i = roots.length; i < (roots.length + trailing_zeros); ++i) {
                r2[i] = new Complex(0, 0);
            }
            roots = r2;
        }
        
        return roots;
    }
    
    public static int[] nonzero(double[] p, double eps) {
        TIntList idx = new TIntArrayList();
        for (int i = 0; i < p.length; ++i) {
            if (Math.abs(p[i]) > eps) {
                idx.add(i);
            }
        }
        return idx.toArray();
    }
    
    public static double[] realRoots(double[] p) throws NotConvergedException {
        
        Complex[] roots = roots(p);
        
        if (roots == null || roots.length == 0) {
            return new double[]{};
        }
        
        double tol = 1e-5;
        
        // only returning real roots;
        TDoubleList realRoots = new TDoubleArrayList();
        int n = 0;
        for (int i = 0; i < roots.length; ++i) {
            if (Math.abs(roots[i].im()) < tol) {
                realRoots.add(roots[i].re());
            }
        }

        return realRoots.toArray();
    }
}
