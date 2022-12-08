package algorithms.dimensionReduction;

import algorithms.imageProcessing.SummedAreaTable0;
import algorithms.matrix.MatrixUtil;
import algorithms.matrix.MatrixUtil.SVDProducts;
import algorithms.statistics.CDFRandomSelect;
import algorithms.misc.Misc0;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import gnu.trove.map.TIntIntMap;
import java.util.Arrays;
import java.util.Random;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 * a CUR-Decomposition is an approximation to the SVD, but the result is that 
 * all large matrices of the decomposition are sparse in contrast to the same
 * in SVD.
 * 
 * This class implements CUR-decomposition following the book "Mining of Massive 
 * Datasets" by Jure Leskovec, Anand Rajaraman, Jeff Ullman.
 * http://www.mmds.org/
 * 
 * "CUR matrix decompositions for improved data analysis"
 * Runtime is O(m*n) Drineas et al. 2008
 * (https://www.pnas.org/content/pnas/106/3/697.full.pdf)
 * "CUR decompositions are low-rank matrix decompositions that are explicitly 
 * expressed in terms of a small number of actual columns and/or actual rows of 
 * the data matrix...
 * We present an algorithm that preferentially chooses columns and rows that 
 * exhibit high “statistical leverage” and, thus, in a very precise statistical 
 * sense, exert a disproportionately large “influence” on the best low-rank fit 
 * of the data matrix. By selecting columns and rows in this manner, we obtain 
 * improved relative-error and constant-factor approximation guarantees in 
 * worst-case analysis, as opposed to the much coarser additive-error guarantees 
 * of prior work. In addition, since the construction involves computing 
 * quantities with a natural and widely studied statistical interpretation, 
 * we can leverage ideas from diagnostic regression analysis to employ these 
 * matrix decompositions for exploratory data analysis."
 * 
 * 
 * NOTE: consider also implementing in future:
   http://www.cs.cmu.edu/~christos/PUBLICATIONS/sdm07-lsm.pdf
   "Less is More: Compact Matrix Decomposition for Large Sparse Graphs"
   2007, Sun, Xie, Zhang, and Faloutsos
   Proceedings of the 2007 SIAM International Conference on Data Mining
   ...the Compact Matrix Decomposition (CMD), to compute sparse low rank 
   * approximations. CMD dramatically reduces both the computation cost and 
   * the space requirements over existing decomposition methods (SVD, CUR). 
   * Using CMD as the key building block, we further propose procedures to 
   * efficiently construct and analyze dynamic graphs from real-time 
   * application data. ...
   * Also, kernel PCA, nonmetric multidimensional scaling (NMDS), and a fast
   * version of isomap.
   * 
   * 
 * @author nichole
 */
public class CURDecomposition {
    
    /**
     * a CUR-Decomposition is an approximation to the SVD, but the result is that 
     * all large matrices of the decomposition are sparse in contrast to the same
     * in SVD.
     * 
     *  Runtime is O(m*n) Drineas et al.
     * 
     * NOTE: Let C = CUR decomposition(A).  The eigenvectors of C are equal
     * to the eigenvectors of A and both are equal to the eigenvectors
     * C^T*C and A^T*A.
     * Let lambda_CC = eigenvalues of C^T*C.
     * The determinant( C^T*C - lambda * I) = 0.
     * The trace of (C^T*C) = trace(lambda * I).
     * 
     * NOTE: could improve the speed and space used for very large matrices
     * by using FlexCompColMatrix and FlexCompRowMatrix for C, U, R, and 
     * result matrices.
     * Also, for very large matrices, consider implementing in CDFRandomSelect.java
     * the integer transformation and storage in YFastTrie.
     * 
     @param a is an mxn matrix.
     @param k the rank to approximate
     @return the cur decomposition C, U, and R.  c plays the role of svd's U
     * and R plays he rol of SVD's v^T.
     * @throws no.uib.cipr.matrix.NotConvergedException thrown if the MTJ SVD method cannot converge
     * for a matrix.
     */
    public static CUR calculateDecomposition(double[][] a, int k) throws NotConvergedException {

        CDFs cdfs = _calculateCDFs(a, k);
       
        SelectedFromA r = _calculateR(a, cdfs.rowsSelected, cdfs.pdfs.rowPDF);
        
        SelectedFromA c = _calculateR(MatrixUtil.transpose(a), cdfs.colsSelected, cdfs.pdfs.colPDF);
        c.r = MatrixUtil.transpose(c.r);
        
        double[][] u = _calculateU(a, r.indexesUnique, c.indexesUnique);
        
        CUR cur = new CUR();
        cur._rowsSelected = r.indexesUnique;
        cur._colsSelected = c.indexesUnique;
        cur.c = c.r;
        cur.r = r.r;
        cur.u = u;
        cur.result = MatrixUtil.multiply(cur.getC(), cur.getU());
        cur.result = MatrixUtil.multiply(cur.getResult(), cur.getR());
        
        return cur;        
    }
    
    static CDFs _calculateCDFs(double[][] a, int k) {
     
        PDFs pdfs = _calculatePDFs(a);
        
        // randomly select k columns and k rows
        
        double[] rowCDF = MiscMath0.cumulativeSum(pdfs.rowPDF);
        
        double[] colCDF = MiscMath0.cumulativeSum(pdfs.colPDF);
        
        // normalize so that the last bin is "1".
        double norm = rowCDF[rowCDF.length - 1];
        for (int i = 0; i < rowCDF.length; ++i) {
            rowCDF[i] /= norm;
        }
        
        norm = colCDF[colCDF.length - 1];
        for (int i = 0; i < colCDF.length; ++i) {
            colCDF[i] /= norm;
        }
                        
        int[] k0, k1;
        
        boolean useBinarySearch = true;
        
        /*
        In the future, may implement more methods than binary search in
        CDFRandomSelect.
        
        (1) could provide an overloaded method for N_CDF > 32767 
            that uses a YFastTrie, but would have to convert the CDF's real values
            to integers of a comparetively small range such as 32767 for the 
            trie and the same for the conversion of the random number [0,1].
        (2) a faster query time may exist for a hash structure like an LSH?
           or any other constant time query for Nearest Neighbor or 
           Approx Nearest Neighbor?
        (3) the lecture slides from http://www.mmds.org
            suggest an algorithm they call ColumnSelect which has a runtime of
            O(k *lg(k) / eps^2 )
        */
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        System.out.println("seed=" + seed);
        System.out.flush();
        rand.setSeed(seed);

        //if (useBinarySearch) {
            k0 = CDFRandomSelect.chooseKFromBinarySearch(rowCDF, k, rand);
            k1 = CDFRandomSelect.chooseKFromBinarySearch(colCDF, k, rand);
        /*} else {
            // transform CDFs to integers
            k0 = CDFRandomSelect.chooseKFromIntegerTransformAndTrie(rowCDF, k, rand);
            k1 = CDFRandomSelect.chooseKFromIntegerTransformAndTrie(colCDF, k, rand);
        }*/

        CDFs cds = new CDFs();
        cds.rowsSelected = k0;
        cds.colsSelected = k1;
        cds.rowCDF = rowCDF;
        cds.colCDF = colCDF;
        cds.pdfs = pdfs;

        return cds;
    }
    
    /**
     * calculate the col and row Frobenius norm discrete probabilities from a
     * (a.k.a. calculating Marginal probability mass function of the contingency table).
     * runtime is O(N) where N = mxn.
     * NOTE: the Frobenius norm is the square root of the sum of the squares of 
     * all elements of a matrix (2.2-4 of Golub and van Loan).
     @param a an mxn matrix.
     @return  the column and row PDFs of matrix a where a row PDF is the
     * Frobenius norm of the column divided by the Frobenius norm of the
     * matrix, and the row PDF is similar but calculated for rows instead of
     * columns.
     */
    static PDFs _calculatePDFs(double[][] a) {
     
        // copy a, square each item, create a summed area table from that.
        // create column sums as a vector and normalize it (= discrete pdf for col)
        // create row sums as a vector and normalize it (= discrete pdf for row)
        
        a = MatrixUtil.copy(a);
        int i, j;
        for (i = 0; i < a.length; ++i) {
            for (j = 0; j < a[i].length; ++j) {
                a[i][j] *= a[i][j];
            }
        }
        
        SummedAreaTable0 sat0 = new SummedAreaTable0();
        double[][] s = sat0.createAbsoluteSummedAreaTable(a);
        
        double fN = s[s.length-1][s[0].length-1];
        
        double[] r = new double[s.length];
        double[] c = new double[s[0].length];
        
        /*
        imgS 
        startX - coordinate for x start of window 
        stopX - coordinate for x stop of window 
        startY - coordinate for y start of window 
        stopY - coordinate for y stop of window 
        output - one dimensional array of size 2 in which the sum of the window 
             will be returned and the number of pixels in the window. 
             double[]{sum, nPixels}
        */
        double[] sumAndN = new double[2];
        
        int startCol = 0;
        int stopCol = s[0].length - 1;
        // create row vectors:
        for (i = 0; i < s.length; ++i) {
            
            // s[x:x][y:y]
            sat0.extractWindowFromSummedAreaTable(s, i, i, startCol, stopCol,
                sumAndN);
            
            r[i] = sumAndN[0]/fN;
        }
        
        int startRow = 0;
        int stopRow = s.length - 1;
        // create col vectors:
        for (j = 0; j < s[0].length; ++j) {
            
            sat0.extractWindowFromSummedAreaTable(s, startRow, stopRow, j, j, 
                sumAndN);
            
            c[j] = sumAndN[0]/fN;
        }
        
        PDFs pdfs = new PDFs();
        pdfs.colPDF = c;
        pdfs.rowPDF = r;

        return pdfs;        
    }

    /**
      calculate R: 
      R is an k x n matrix composed of a randomly chosen set of k rows of M.
      each row of R gets normalized by sqrt( k * rowPDF_{selected_row}).
      if a row of A is present more than once in the selection, it is 
      included only the first time and that included column is then multiplied
      by the square root of the number of duplicates.
      @param a the m x n matrix, of m samples of n variates or parameters.
      @param selectedRowIdxs a vector of indexes of rows of matrix a, which have
      been randomly selected from the row marginal probabilities.
      @param rowPDF the row marginal probabilities.
      @return matrix of selected rows from A normalized and an array of the
      unique selected indexes (corrections for duplicated indexes are performed
      in this method).
     */
    static SelectedFromA _calculateR(double[][] a, int[] selectedRowIdxs, 
        double[] rowPDF) {
        
        if (a.length != rowPDF.length) {
            throw new IllegalArgumentException("a.length must equal rowPDF.length");
        }
        
        // make frequency map of selected indexes and remove the key after a
        //   column is visited (to help correct for multiplicity).
        TIntIntMap freq = MiscMath0.makeFrequencyMap(selectedRowIdxs);
        
        int kCorr = freq.size();
        
        int k = selectedRowIdxs.length;
        int m = a.length;
        int n = a[0].length;
        
        int i, j, d, ii;
        int rIdx = 0;
        double factor;
        
        SelectedFromA sa = new SelectedFromA();
        sa.r = new double[kCorr][n];
        sa.indexesUnique = new int[kCorr];
        
        for (i = 0; i < selectedRowIdxs.length; ++i) {
            
            j = selectedRowIdxs[i];
            assert(j >= 0 && j < m);
            
            if (!freq.containsKey(j)) {
                continue;
            }
            
            d = freq.get(j);
            
            factor = 1./Math.sqrt(k * rowPDF[j]);
            
            if (d > 1) {
                factor *= Math.sqrt(d);
            }
            
            //extract row j from a and multiply each item by factor
            //and store it as c[cIdx][*];
            sa.r[rIdx] = Arrays.copyOf(a[j], a[j].length);
            for (ii = 0; ii < sa.r[rIdx].length; ++ii) {
                sa.r[rIdx][ii] *= factor;
            }
            
            sa.indexesUnique[rIdx] = j;
            
            freq.remove(j);
                    
            rIdx++;
        }
        
        return sa;
    }

    static double[][] _calculateU(double[][] a, int[] selectedRIdxs, 
        int[] selectedCIdxs) throws NotConvergedException {
        
        int nr = selectedRIdxs.length;
        int nc = selectedCIdxs.length;
        
        int i, j, ar, ac;
        double[][] w = new double[nr][nc];
        for (i = 0; i < nr; ++i) {
            w[i] = new double[nc];
            ar = selectedRIdxs[i];
            for (j = 0; j < nc; ++j) {
                ac = selectedCIdxs[j];
                w[i][j] = a[ar][ac];
            }
        }
        
        SVD svd = SVD.factorize(new DenseMatrix(w));
        DenseMatrix X = svd.getU();
        DenseMatrix YT = svd.getVt();
        double[] s = svd.getS();
        
        double eps = 1.e-15;
        
        int rank = 0;
        for (i = 0; i < s.length; ++i) {
            if (Math.abs(s[i]) > eps) {
                rank++;
            }
        }
        
        if (rank == 0) {
            throw new IllegalArgumentException("factorization cannot continue as"
            + " the number of selected indexes is too small");
        }
        
        // this will be transposed
        double[][] sInvSq = MatrixUtil.zeros(nr, nc);
        int len = Math.min(nr, nc);
        for (i = 0; i < len; ++i) {
            if (s[i] > eps) {
                sInvSq[i][i] = 1./(s[i]*s[i]);
            }
        }
        sInvSq = MatrixUtil.transpose(sInvSq);
        
        // X is nr x nr
        // S is nr x nc ==> (Σ+)^2 = nc x nr
        // Y is nc x nc
        
        //U = Y * (Σ+)^2 * X^T
        //    [  nc x nr ] * [nr x nr] 
        //  = [nc][nr]
                
        double[][] y = MatrixUtil.convertToRowMajor(YT);
        y = MatrixUtil.transpose(y);
        double[][] xt = MatrixUtil.convertToRowMajor(X);
        xt = MatrixUtil.transpose(xt);
        double[][] u = MatrixUtil.multiply(y, sInvSq);
        u = MatrixUtil.multiply(u, xt);
        
       return u;
    }
    
    /**
     *
     */
    public static class CUR {
        private double[][] c;
        private double[][] u;
        private double[][] r;
        private double[][] result;
        int[] _colsSelected;
        int[] _rowsSelected;
        
        private SVDProducts svd = null;

        /**
         @return the c
         */
        public double[][] getC() {
            return c;
        }

        /**
         @return the u
         */
        public double[][] getU() {
            return u;
        }

        /**
         @return the r
         */
        public double[][] getR() {
            return r;
        }

        /**
         @return the result
         */
        public double[][] getResult() {
            return result;
        }
        
        /**
         *
         @return
         */
        public SVDProducts getApproximateSVD() {
            if (svd != null) {
                return svd;
            }
            SVDProducts svd2 = new SVDProducts();
            svd2.u = MatrixUtil.copy(c);
            svd2.vT = MatrixUtil.copy(r);
            
            // normalize U by columns and V^T by rows
            MatrixUtil.normalizeColumnsL2(svd2.u);
            MatrixUtil.normalizeRowsL2(svd2.vT);
            
            /*
            a = result which is mXn
            u is mXm
            vT is nXn
            s is mXn
            
            A*V = S*U
            S = U^T*A*V
            */
            int m = svd2.u.length;
            int n = svd2.vT.length;
            
            svd2.sigma = MatrixUtil.multiply(MatrixUtil.transpose(svd2.u), result);
            svd2.sigma = MatrixUtil.multiply(svd2.sigma, MatrixUtil.transpose(svd2.vT));

            int sn = Math.min(svd2.sigma.length, svd2.sigma[0].length);
            svd2.s = new double[sn];
            for (int i = 0; i < sn; ++i) {
                svd2.s[i] = svd2.sigma[i][i];// sigma is [m X n]
            }
            
            svd = svd2;
            
            return svd;
        }
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("c=\n");
            if (c != null) {
                sb.append(FormatArray.toString(c, "%.4e"));
            } 
            sb.append("\nu=\n");
            if (u != null) {
                sb.append(FormatArray.toString(u, "%.4e"));
            }
            sb.append("\nr=\n");
            if (r != null) {
                sb.append(FormatArray.toString(r, "%.4e"));
            }
            sb.append("\ncur=\n");
            if (result != null) {
                sb.append(FormatArray.toString(result, "%.4e"));
            }
            return sb.toString();
        }
    }
    
    /**
     *
     */
    public static class PDFs {
        double[] colPDF;
        double[] rowPDF;
    }
    
    /**
     *
     */
    public static class CDFs {
        int[] colsSelected;
        int[] rowsSelected;
        PDFs pdfs;
        double[] colCDF;
        double[] rowCDF;
    }
    
    /**
     *
     */
    public static class SelectedFromA {
        int[] indexesUnique;
        double[][] r;
    }
}
