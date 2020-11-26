package algorithms.graphs;

import algorithms.util.SimpleLinkedListNode;

/**
Newman 2006, "Modularity and community structure in networks"
https://arxiv.org/pdf/physics/0602124.pdf

modularity is a graph quality function over possible divisions of a network.
...
The modularity is, up to a multiplicative constant, the number of edges
 * falling within groups minus the expected number in an equivalent network with
 * edges placed at random. (A precise mathematical formulation is given below.)
 * The modularity can be either positive or negative, with positive values
 * indicating the possible presence of community structure. Thus, one can search
 * for community structure precisely by looking for the divisions of a network
 * that have positive, and preferably large, values of the modularity
 
 NOTE: a paper with tables summarizing different modularity algorithms for
     different definitions of communities as disjoint or sharing members is
     Table 2 on page 5 of
     https://arxiv.org/pdf/1708.00977.pdf
      
     
 <pre>
 Newman 2006:
    defines modularity using eignevectors of a characteritic matrix
    modularity Q = 
      

 Girvan-Newman 2002, 2004:
    modularity Q = summation_{over i}( eii - ai^2 ) = Tr e - ||e^2||
       where ||e^2|| is the sum of the square of elements of matrix e.
       where matrix e is a k x k matrix with each 
           element eij being the fraction of network edges connecting community i to j.
           matrix e is upper right triangle or eij is split in half between 
               eij and eji for symmetry (either way ensures i:j gets counted only once).
       
       where Tr e = summation_{over i}( eii )
           Tr e gives fraction of network edges in same community connecting.
           this should be large for good community divisions, but is not
           a good indicator of the quality of the division.
       
       where ai = summation_{over j}( eij ) for columns
           and represents the fraction of edges that connect to community i.

       where bj = summation_{over i}( eij ) for rows
           and represents the fraction of edges that connect to community i.

   in a network in which edges fall between vertices without regard for the 
   communities they belong to, we would have
       eij = ai*aj.
   
   If the number of within-community edges is no better than random, we will get
       Q = 0.
     value approaching the maximum,
       Q = 1, indicate strong community structure [50].
     values tend to be 0.3 to 0.7 and higher values are rare.
     
Duch and Arenas 2005 "Community detection in complex networks using extremal 
optimization"
    Q =
     
  
 Clauset, Newman, Moore, 2004
    Q = 


Nadakudiki and Newman 2012, "Graph spectra and the detectability of community 
structure in networks"
https://arxiv.org/pdf/1205.1813.pdf

 </pre>
 J. Duch and A. Arenas, Community detection in complex
networks using extremal optimization. Phys. Rev. E 72,
027104 (2005).

A. Clauset, M. E. J. Newman, and C. Moore, Finding
community structure in very large networks. Phys. Rev.
E 70, 066111 (2004).


NOTE: Newman generates random graphs for testing too:
   https://arxiv.org/pdf/cond-mat/0209450.pdf
   https://arxiv.org/pdf/cond-mat/0210146.pdf
   
 * @author nichole
 */
public class Modularity {
    
    public static class ModularityResults {
        float q;
    }
    
    public ModularityResults girvanNewman2002(SimpleLinkedListNode[] adjacencyList) {
        
        // NOTE: communities are disjoint sets for Girvan-Newman 2002 and 2004
       
        /*
        modularity Q = summation_{over i}( eii - ai^2 ) = Tr e - ||e^2||
       where ||e^2|| is the sum of the square of elements of matrix e.
       where matrix e is a k x k matrix with each 
           element eij being the fraction of network edges connecting community i to j.
           matrix e is upper right triangle or eij is split in half between 
               eij and eji for symmetry (either way ensures i:j gets counted only once).
       
       where Tr e = summation_{over i}( eii )
           Tr e gives fraction of network edges in same community connecting.
           this should be large for good community divisions, but is not
           a good indicator of the quality of the division.
       
       where ai = summation_{over j}( eij ) for columns
           and represents the fraction of edges that connect to community i.

       where bj = summation_{over i}( eij ) for rows
           and represents the fraction of edges that connect to community i.

   in a network in which edges fall between vertices without regard for the 
   communities they belong to, we would have
       eij = ai*aj.
   
   If the number of within-community edges is no better than random, we will get
       Q = 0.
     value approaching the maximum,
       Q = 1, indicate strong community structure [50].
     values tend to be 0.3 to 0.7 and higher values are rare.
      */
        throw new UnsupportedOperationException("not yet implemented");
    }
}
