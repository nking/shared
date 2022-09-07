package algorithms.matrix;

import java.util.Arrays;

/**
 * a class to hold a matrix that is partitioned into blocks of the same size.
 * @author nichole
 */
public class BlockMatrixIsometric {
    /**
     * the matrix
     */
    private final double[][] a;
    /**
     * the blocks' sizes in the row dimension
     */
    private final int bSize0; 
    
    /**
     * the blocks' sizes in the column dimension
     */
    private final int bSize1;
    
    /**
     * constructor.  The matrix a is not copied in.   This is a 
     * pass-by-reference constructor.
     * @param a
     * @param bSize0 block size along the first dimension (rows)
     * @param bSize1 block size along the second dimension (columns)
     */
    public BlockMatrixIsometric(double[][] a, int bSize0, int bSize1) {
        int n = a.length;
        int m = a[0].length;
        if ((n % bSize0) != 0) {
            throw new IllegalArgumentException("a.length must be evenly divided"
                    + " by bSize0");
        }
        if ((m % bSize1) != 0) {
            throw new IllegalArgumentException("a[0].length must be evenly divided"
                    + " by bSize1");
        }
        this.a = a;
        this.bSize0 = bSize0;
        this.bSize1 = bSize1;
    }
    
    /**
     * set all values to 0
     */
    public void reset() {
        fill(0);
    }
    
    void fill(double value) {
        int i;
        for (i = 0; i < a.length; ++i) {
            Arrays.fill(a[i], value);
        }
    }
    
    /**
     * set values in internal matrix a for block number (blockNumber0, blockNumber1) 
     * to the given block b.  b must be the same size as the block size
     * of this instance.
     * @param b a given block of values to replace a block of this instance with.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void setBlock(double[][] b, int blockNumber0, int blockNumber1) {
        if (b.length != bSize0 || b[0].length != bSize1) {
            throw new IllegalArgumentException("b size must equal the block size"
                    + " within this instance");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        for (int i = 0; i < this.bSize0; ++i) {
            System.arraycopy(b[i], 0, a[start0 + i], start1, bSize1);
        }
    }
    
    /**
     * set values in internal matrix a for block number (blockNumber0, blockNumber1) 
     * to the given row block b which is a single dimension array.  
     * b.length must equal this.bSize1, and this.bSize0 must equal 1.
     * @param b a given block of values in a single dimension array to replace a 
     * block of this instance with.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void setRowBlock(double[] b, int blockNumber0, int blockNumber1) {
        if (bSize0 != 1) {
            throw new IllegalArgumentException("block b is a row array "
                    + " and so this instance, when constructed, should "
                    + " have used bSize0=1");
        }
        if (b.length != bSize1) {
            throw new IllegalArgumentException("this.bSize1 must equal b.length");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        
        //a[start0][start1] = b
        System.arraycopy(b, 0, a[start0], start1, b.length);
    }
    
    /**
     * set values in internal matrix a for block number (blockNumber0, blockNumber1) 
     * to the given column block b which is a single dimension array.  
     * b.length must equal this.bSize0, and this.bSize1 must equal 1.
     * @param b a given block of values in a single dimension array to replace a 
     * block of this instance with.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void setColumnBlock(double[] b, int blockNumber0, int blockNumber1) {
        if (bSize1 != 1) {
            throw new IllegalArgumentException("block b is a column array "
                    + " and so this instance, when constructed, should "
                    + " have used bSize1=1");
        }
        if (b.length != bSize0) {
            throw new IllegalArgumentException("this.bSize0 must equal b.length");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
                        
        for (int i = 0; i < this.bSize0; ++i) {
            a[start0 + i][start1] = b[i];
        }
    }
    
    /**
     * get values in internal matrix a for block number (blockNumber0, blockNumber1) 
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     * @return contents of (blockNumber0, blockNumber1) of matrix a in this instance.
     */
    public double[][] getBlock(int blockNumber0, int blockNumber1) {
        double[][] b = MatrixUtil.zeros(bSize0, bSize1);
        getBlock(b, blockNumber0, blockNumber1);
        return b;
    }
    
    /**
     * get values in internal matrix a for block number (blockNumber0, blockNumber1) 
     * and set them into the given block b.
     * @param b the output block for the contents of (blockNumber0, blockNumber1) 
     * of matrix a in this instance.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void getBlock(double[][] b, int blockNumber0, int blockNumber1) {
        if (b.length != bSize0 || b[0].length != bSize1) {
            throw new IllegalArgumentException("b size must equal the block size"
                    + " within this instance");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        for (int i = 0; i < this.bSize0; ++i) {
            System.arraycopy(a[start0 + i], start1, b[i], 0, bSize1);
        }
    }
    
    /**
     * get values in internal matrix a for block number (blockNumber0, blockNumber1) 
     * and set them into the given row block b.
     * @param b the output row block for the contents of (blockNumber0, blockNumber1) 
     * of matrix a in this instance.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void getRowBlock(double[] b, int blockNumber0, int blockNumber1) {
        if (bSize0 != 1) {
            throw new IllegalArgumentException("expecting this.bSize0=1 for row blocks."
            + " this.bSize0 is set at construction.");
        }
        if (b.length != bSize1) {
            throw new IllegalArgumentException("b.length must equal this.bSize1");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        
        System.arraycopy(a[start0], start1, b, 0, bSize1);
    }
    
    /**
     * get values in internal matrix a for block number (blockNumber0, blockNumber1) 
     * and set them into the given column block b.
     * @param b the output column block for the contents of (blockNumber0, blockNumber1) 
     * of matrix a in this instance.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void getColumnBlock(double[] b, int blockNumber0, int blockNumber1) {
        if (bSize1 != 1) {
            throw new IllegalArgumentException("expecting this.bSize1=1 for column blocks."
            + " this.bSize1 is set at construction.");
        }
        if (b.length != bSize0) {
            throw new IllegalArgumentException("b.length must equal this.bSize0");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        
        for (int i = 0; i < this.bSize0; ++i) {
            b[i] = a[start0 + i][start1];
        }
    }
    
    /**
     * add block b to internal matrix, pointwise, for (blockNumber0, blockNumber1).
     * b must be the same size as the block size of this instance.
     * @param b a given block of values to add to a block of this instance.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void addToBlock(double[][] b, int blockNumber0, int blockNumber1) {
        if (b.length != bSize0 || b[0].length != bSize1) {
            throw new IllegalArgumentException("b size must equal the block size"
                    + " within this instance");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        int i, j;
        
        for (i = 0; i < this.bSize0; ++i) {
            for (j = 0; j < this.bSize1; ++j) {
                a[start0 + i][start1 + j] += b[i][j];
            }
        }
    }
    
    /**
     * add column block b to internal matrix, pointwise, for (blockNumber0, blockNumber1).
     * b must be the same size as the block size of this instance, including
     * this.bSize1 must equal 1 at construction.
     * @param b a given column block of values to add to a block of this instance.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void addToColumnBlock(double[] b, int blockNumber0, int blockNumber1) {
        if (bSize1 != 1) {
            throw new IllegalArgumentException("this.bSize1 must equal 1 for column blocks."
            + "  That is set at construction.");
        }
        if (b.length != bSize0) {
            throw new IllegalArgumentException("b.length must equal this.bSize0");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        int i;
        for (i = 0; i < b.length; ++i) {
            a[start0 + i][start1] += b[i];
        }
    }
    
    /**
     * add row block b to internal matrix, pointwise, for (blockNumber0, blockNumber1).
     * b must be the same size as the block size of this instance, including
     * this.bSize0 must equal 1 at construction.
     * @param b a given row block of values to add to a block of this instance.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void addToRowBlock(double[] b, int blockNumber0, int blockNumber1) {
        if (bSize0 != 1) {
            throw new IllegalArgumentException("this.bSize0 must equal 1 for row blocks."
            + "  That is set at construction.");
        }
        if (b.length != bSize1) {
            throw new IllegalArgumentException("b.length must equal this.bSize1");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        int i;
        for (i = 0; i < b.length; ++i) {
            a[start0][start1 + i] += b[i];
        }
    }
    
    /**
     * subtract block b from internal matrix, pointwise, for (blockNumber0, blockNumber1).
     * b must be the same size as the block size of this instance.
     * @param b a given block of values to replace a block of this instance with.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void subtractFromBlock(double[][] b, int blockNumber0, int blockNumber1) {
        if (b.length != bSize0 || b[0].length != bSize1) {
            throw new IllegalArgumentException("b size must equal the block size"
                    + " within this instance");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        int i, j;
        
        for (i = 0; i < this.bSize0; ++i) {
            for (j = 0; j < this.bSize1; ++j) {
                a[start0 + i][start1 + j] -= b[i][j];
            }
        }
    }
    
    /**
     * subtract column block b from internal matrix, pointwise, for (blockNumber0, blockNumber1).
     * b must be the same size as the block size of this instance, including
     * this.bSize1 must equal 1 at construction.
     * @param b a given column block of values to subtract from a block of this instance.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void subtractFromColumnBlock(double[] b, int blockNumber0, int blockNumber1) {
        if (bSize1 != 1) {
            throw new IllegalArgumentException("this.bSize1 must equal 1 for column blocks."
            + "  That is set at construction.");
        }
        if (b.length != bSize0) {
            throw new IllegalArgumentException("b.length must equal this.bSize0");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        int i;
        for (i = 0; i < b.length; ++i) {
            a[start0 + i][start1] -= b[i];
        }
    }
    
    /**
     * subtract row block b from internal matrix, pointwise, for (blockNumber0, blockNumber1).
     * b must be the same size as the block size of this instance, including
     * this.bSize0 must equal 1 at construction.
     * @param b a given row block of values to subtract from a block of this instance.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void subtractFromRowBlock(double[] b, int blockNumber0, int blockNumber1) {
        if (bSize0 != 1) {
            throw new IllegalArgumentException("this.bSize0 must equal 1 for row blocks."
            + "  That is set at construction.");
        }
        if (b.length != bSize1) {
            throw new IllegalArgumentException("b.length must equal this.bSize1");
        }
        if ((blockNumber0 * bSize0) >= a.length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= a[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        int i;
        for (i = 0; i < b.length; ++i) {
            a[start0][start1 + i] -= b[i];
        }
    }
    
    public BlockMatrixIsometric copy() {
        BlockMatrixIsometric c = new BlockMatrixIsometric(a, bSize0, bSize1);
        c.set(this);
        return c;
    }
    
    /**
     * set internal contents to equal those of the given b where b has the
     * same dimensions and block size.
     * @param b 
     */
    public void set(BlockMatrixIsometric b) {
        if (b.a.length != a.length) {
            throw new IllegalArgumentException("");
        }
        if (b.a[0].length != a[0].length) {
            throw new IllegalArgumentException("");
        }
        if (b.bSize0 != bSize0) {
            throw new IllegalArgumentException("");
        }
        if (b.bSize1 != bSize1) {
            throw new IllegalArgumentException("");
        }
        int i;
        for (i = 0; i < a.length; ++i) {
            System.arraycopy(b.a[i], 0, a[i], 0, a[i].length);
        }
    }
    
    /**
     * get the entire matrix of this instance.  Note that this is a reference
     * to the instance and is not a copy of it, so modifications of the
     * returned matrix affect this instance too.
     * @return the a
     */
    public double[][] getA() {
        return a;
    }

    /**
     * get the size of dimension 0 of the blocks in this instance
     * @return the bSize0
     */
    public int getBlockSize0() {
        return bSize0;
    }

    /**
     * get the size of dimension 1 of the blocks in this instance
     * @return the bSize1
     */
    public int getBlockSize1() {
        return bSize1;
    }
    
}
