package thirdparty.fr.inria.optimization.cmaes.fitness;

/** Minimalistic interface of a single-objective function 
 * (fitness function) to be minimized. 
 * 
 * see CMAEvolutionStrategy.java for copyright.
*/
public interface IObjectiveFunction {
    /*@param x  a point (candidate solution) in the pre-image of the objective function 
        @return  objective function value of the input search point  
     */

    /**
     *
     @param x
     @return
     */

    double valueOf(double x[]);

    /**
     *
     @param x
     @return
     */
    boolean isFeasible(double x[]);
}

