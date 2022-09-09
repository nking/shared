package algorithms.sort;

/**
 *
 * @author nichole
 */
public class InsertionSort {
    
    /**
     * sorts in place with worse case runtime of O(n^2), best case of O(n) and
     * average case O(n^2).
     * 
     * following pseudocode in Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms" and
     * wikipedia.
     * 
     * @param a 
     */
    public static void sort(int[] a) {
        
        int i = 1;
        int j;
        int n = a.length;
        int x;
        
        while (i < n) {
            x = a[i];
            j = i - 1;
            while (j >= 0 && a[j] > x) {
                a[j+1] = a[j];
                j--;
            }
            a[j+1] = x;
            i++;
        }
    }
                            
}
