package algorithms.matrix;

import algorithms.imageProcessing.SummedAreaTable0;
import algorithms.misc.CDFRandomSelect;
import algorithms.misc.Misc0;
import algorithms.misc.MiscMath0;
import java.util.Arrays;
import java.util.Random;

/**
 * a CUR Decomposition is an approximation to the SVD, but the result is that 
 * all large matrices of the decomposition are sparse in contrast to the same
 * in SVD.
 * 
 * This class implements CUR decomposition following the book "Mining of Massive 
 * Datasets" by Jure Leskovec, Anand Rajaraman, Jeff Ullman.
 * 
 * 
 * NOTE: in contrast to SVD and CUR, consider implementing in future:
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
    
    public static CUR calculateDecomposition(double[][] a, double k) {

        /*
        *Note this is a randomized algorithm, same column can be sampled more than once and that 
         corrections for such redundancies are made at end of algorithm.
          M: a mxn matrix
          r: Pick a target number of “concepts” to be used in the decomposition.
          CUR-decomposition of M
            C is mxr matrix composed of a randomly chosen set of r columns of M.
            R is rxn matrix composed of a randomly chosen set of r rows of M.
            U is rxr matrix constructed from C and R as follows:
       1. Let W be the r×r matrix that is the intersection of the chosen
          columns of C and the chosen rows of R.
          That is, the element in row i and column j of W
          is the element of M whose
            column is the jth column of C and whose row is the ith row of R.

          NOTE: the choosing of random columns and rows has to be biased by
            importance of the row or column, and that is gauged by
            a column based frobenius norm and a row based frobenius norm
          Algorithm design notes:
            - can create matrix colRowF in O(N) by copying M, squaring each item,
                then creating a summed area table from that.
                any rectangle extracted from the summed area table only takes 4 steps.
            - then can create the probability vectors for columns and rows from the colRowF
            - then can create cdfs for those 2 vectors, normalized to peak at sum (last bin) = 1.
            - then the random selection of [0,1] is chosen and one must find that as a
              probability in a bin of the col or row probability CDF vector (and hence the col or row from that).
              -- the data structure used to store the CDFs have a fast nearest neighbor query for
                 manhattan distance, one-dimension.
                 * a balanced binary search tree has an O(lg(N_CDF)) query
                 * could provide an overloaded method for N_CDF > 32767 
                    that uses YFastTrie, but would have to convert the CDFs real values
                    to integers of a competetively small range such as 32767 for the trie and
                    the same for the conversion of the random number [0,1].
                 * a faster query time may exist for a hash structure like an LSH?
                   or any other constant time query for Nearest Neighbor or Approx Nearest Neighbor?
              -- the lecture slides from the data mining book suggest an algorithm they call 
                 ColumnSelect which has a runtime of O(k *lg(k) / eps^2 )

        2. Compute the SVD of W; say W = XΣY^T.
        3. Compute Σ+, the Moore-Penrose pseudoinverse of the diagonal matrix Σ.
           That is, if the ith diagonal element of Σ is σ 6= 0, then replace it by
           1/σ. But if the ith element is 0, leave it as 0.
        4. Let U = Y (Σ+)^2 X^T

         *Constructing C and R:
           select r columns for the matrix C randomly from the columns of M, with probability q_j
        Having selected each of the columns of M, we scale each column by dividing
        its elements by the square root of the expected number of times this column
        would be picked. That is, we divide the elements of the jth column of M, if it
        is selected, by sqrt(r*q_j) . The scaled column of M becomes a column of C (NOTE: the operation 
        is not retained in M, that is, when row operations begin, they operate on the same M from 
        Fig 11.12, and not on an M altered by column selection.)
            Rows of M are selected for R in the analogous way. For each row of R we
        select from the rows of M, choosing row i with probability p_i.  scale each chosen row by dividing
        by sqrt(r*p_i) if it is the ith row of M that was chosen.

        *Construct U:
        
        *Handle Duplicate columns or rows:
        One can combine k rows of R that are each the same row of the matrix M into a single row of R, 
        thus leaving R with fewer rows. 
        Likewise, k columns of C that each come from the same column of M can be combined into one column of C. 
        However, for either rows or columns, the remaining combined vector should have each of its elements 
        multiplied by sqrt(k). 
        When we merge some rows and/or columns, it is possible that R has fewer rows than C has columns, 
        or vice versa. As a consequence, W will not be a square matrix, however, the pseudoinverse still 
        works though there will be entries with 0 in the diagonal.
        */
        throw new UnsupportedOperationException("implementation is not yet finished.");        
    }
    
    public static CDFs _calculateCDFs(double[][] a, int k) {
     
        PDFs pdfs = _calculatePDFs(a);
        
        double[] rowCDF = MiscMath0.cumulativeSum(pdfs.rowPDF);
        
        double[] colCDF = MiscMath0.cumulativeSum(pdfs.colPDF);
        
        int[] k0, k1;
        
        boolean useBinarySearch = true;
        
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
     * calculate the col and row Frobenius norm discrete probabilities from a.
     * runtime is O(N) where N = mxn.
     * @param a an mxn matrix.
     * @return  the column and row PDFs of matrix a where a row PDF is the
     * Frobenius norm of the column divided by the Frobenius norm of the
     * matrix, and the row PDF is similar but calculated for rows instead of
     * columns.
     */
    public static PDFs _calculatePDFs(double[][] a) {
     
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

    public static class CUR {
        double[] c;
        double[] u;
        double[] r;
        double[] result;
        double[] _colsSelected;
        double[] _rowsSelected;
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
}
