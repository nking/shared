package algorithms.graphs;

import algorithms.matrix.MatrixUtil;
import algorithms.util.SimpleLinkedListNode;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVectorSub;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.SVD;
import no.uib.cipr.matrix.sparse.ArpackSym;
import no.uib.cipr.matrix.sparse.LinkedSparseMatrix;

import java.util.Map;

/**
 *
 * @author nichole
 */
public class Laplacian {

    //TODO: include normalizations and pre-weighted edges

    /**
     * calculate L = D - A where D is a diagonal matrix.
     * For L_i_j, D_i_i is the number of edges into node i and A_i_j is the adjacency of nodes i and j (== 1 if adjacent, else 0).
     @param g
     @return laplacian for D holding in-degrees.
     */
    public static double[][] createInDegreeLaplacian(TIntObjectMap<TIntSet> g) {
        // to get edges into a vertex:
        TIntObjectMap<TIntSet> rev = MatrixUtil.createReverseMap(g);
        return createOutDegreeLaplacian(rev);
    }

    /**
     *
     @param g
     @return
     */
    public static LinkedSparseMatrix createInDegreeLaplacianSparse(TIntObjectMap<TIntSet> g) {
        // to get edges into a vertex:
        TIntObjectMap<TIntSet> rev = MatrixUtil.createReverseMap(g);
        return createOutDegreeLaplacianSparse(rev);
    }

    /**
     *
     @param g
     @return
     */
    public static LinkedSparseMatrix createOutDegreeLaplacianSparse(TIntObjectMap<TIntSet> g) {

        int[] minMax = GraphUtil.minAndMaxVertexNumbers(g);
        int n = minMax[1] + 1;
        LinkedSparseMatrix lS = new LinkedSparseMatrix(n, n);

        int u, v;
        TIntSet vSet;
        TIntObjectIterator<TIntSet> iter = g.iterator();
        TIntIterator iter2;
        while (iter.hasNext()) {
            iter.advance();
            u = iter.key();
            iter2 = iter.value().iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                // add to the in-degree of u
                lS.add(u, u, 1);
                // OR lS.set(u, u, lS.get(u, u) + 1);
                // set u to v -1
                lS.set(u, v, -1);
            }
        }

        return lS;
    }

    /**
     *
     @param g
     @return
     */
    public static double[][] createOutDegreeLaplacian(TIntObjectMap<TIntSet> g) {

        int[] minMax = GraphUtil.minAndMaxVertexNumbers(g);
        int n = minMax[1] + 1;

        double[][] lM = MatrixUtil.zeros(n, n);
        int u, v;
        TIntSet vSet;
        TIntObjectIterator<TIntSet> iter = g.iterator();
        TIntIterator iter2;
        while (iter.hasNext()) {
            iter.advance();
            u = iter.key();
            iter2 = iter.value().iterator();
            while (iter2.hasNext()) {
                v = iter2.next();
                // add to the in-degree of u
                lM[u][u] += 1;
                // set u to v -1
                lM[u][v] = -1;
            }
        }
        return lM;
    }

    /**
     *
     @param g
     @return
     */
    public static double[][] createInDegreeLaplacian(SimpleLinkedListNode[] g) {
        TIntObjectMap<TIntSet> g2 = GraphUtil.convertGraph(g);
        return createInDegreeLaplacian(g2);
    }

    /**
     * calculate the 2nd smallest eigenvector of the Laplacian of undirected graph g.
     * it approximates the smallest cut in the graph.
     @param g
     @return
     */
    public static double[] calculateFieldlerVector(TIntObjectMap<TIntSet> g) {
        LinkedSparseMatrix lS = Laplacian.createInDegreeLaplacianSparse(g);
        //System.out.printf("L=%s\n", lS.toString());

        int nEig = 2;
        ArpackSym.Ritz ritz = ArpackSym.Ritz.SA;
        // calculate the 2 smallest eigenvectors of the laplacian
        Map<Double, DenseVectorSub> eigenVectors = MatrixUtil.sparseEigen(lS, nEig, ritz);

        double maxEig = Double.NEGATIVE_INFINITY;
        for (Map.Entry<Double, DenseVectorSub> eigen : eigenVectors.entrySet()) {
            double eig = eigen.getKey();
            if (eig > maxEig) {
                maxEig = eig;
            }
        }
        double[] eigenVector = Matrices.getArray(eigenVectors.get(maxEig));
        return eigenVector;
    }
}
