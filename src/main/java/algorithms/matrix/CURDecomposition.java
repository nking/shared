package algorithms.matrix;

import algorithms.imageProcessing.SummedAreaTable0;
import algorithms.misc.CDFRandomSelect;
import algorithms.misc.Misc0;
import algorithms.misc.MiscMath0;
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
 * Runtime is O(m*n) Drineas et al.
 * (https://www.pnas.org/content/pnas/106/3/697.full.pdf)
 * 
 * NOTE: consider implementing in future:
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
     * NOTE: could improve the speed and space used for very large matrices
     * by using FlexCompColMatrix and FlexCompRowMatrix for C, U, R, and 
     * result matrices.
     * Also, for very large matrices, consider implementing in CDFRandomSelect.java
     * the integer transformation and storage in YFAstTrie.
     * 
     * @param a is an mxn matrix.
     * @param k
     * @return 
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
        cur.result = MatrixUtil.multiply(cur.c, cur.u);
        cur.result = MatrixUtil.multiply(cur.result, cur.r);
        
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
            to integers of a competetively small range such as 32767 for the 
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

        if (useBinarySearch) {
            k0 = CDFRandomSelect.chooseKFromBinarySearch(rowCDF, k, rand);
            k1 = CDFRandomSelect.chooseKFromBinarySearch(colCDF, k, rand);
        } else {
            // transform CDFs to integers
            k0 = CDFRandomSelect.chooseKFromIntegerTransformAndTrie(rowCDF, k, rand);
            k1 = CDFRandomSelect.chooseKFromIntegerTransformAndTrie(colCDF, k, rand);
        }

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
     * @param a an mxn matrix.
     * @return  the column and row PDFs of matrix a where a row PDF is the
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
        double[][] sInvSq = new double[nr][nc];
        for (i = 0; i < nr; ++i) {
            sInvSq[i] = new double[nc];
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
                
        double[][] y = Matrices.getArray(YT);
        y = MatrixUtil.transpose(y);
        double[][] xt = Matrices.getArray(X);
        xt = MatrixUtil.transpose(xt);
        double[][] u = MatrixUtil.multiply(y, sInvSq);
        u = MatrixUtil.multiply(u, xt);
        
       return u;
    }
    
    public static class CUR {
        double[][] c;
        double[][] u;
        double[][] r;
        double[][] result;
        int[] _colsSelected;
        int[] _rowsSelected;
    }
    
    public static class PDFs {
        double[] colPDF;
        double[] rowPDF;
    }
    
    public static class CDFs {
        int[] colsSelected;
        int[] rowsSelected;
        PDFs pdfs;
        double[] colCDF;
        double[] rowCDF;
    }
    
    public static class SelectedFromA {
        int[] indexesUnique;
        double[][] r;
    }
}
