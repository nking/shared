package algorithms.tsp;

/**
 <pre>
 A completely dynamic solution requires a large amount of memory (see statiement below).
 Here is an outline of one, based upon what I learned from making the hybrid
 dynamic and brute force class:
 
 subset chooser for n, k w/ k=3
   permutations for k!
   store each k=3 permutation in memo
 subset chooser for n, k w/ k=6
   permutations for k!
   path should be composed of complete subparts + a connecting edge:
       p = p3_0 + p3_1 + edge where p3_0 and p3_1 are in memo from k=3 permutations
   store each k=6 permutation in memo
 subset chooser for n, k w/ k=12
   permutations for k!
   path should be composed of complete subparts + a connecting edge:
       p = p6_0 + p6_1 + edge where p6_0 and p6_1 are in memo from k=6 permutations
   store each k=12 permutation in memo
...

For the dynamic approach just outlined:
    n    c             compare to n!
    8    210           40320
   14    1.24e6        8.7e10
   29    1.46e16       8.8e31
  731    5.58e538      7.7e1778
 where n is the number of nodes (=dist.length),
 and c is the number of elements one needs to store in a memo datastructure.
 
</pre>

see AbstractTSP for details about the keys used for memo.

 * @author nichole
 */
public class TSPDynamic extends AbstractTSP {
    
    public TSPDynamic(double[][] dist) {
        super(dist);
    }
}
