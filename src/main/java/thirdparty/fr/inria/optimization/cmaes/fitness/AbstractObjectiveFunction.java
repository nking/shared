package thirdparty.fr.inria.optimization.cmaes.fitness;

/** extending from this abstract class implements a generic isFeasible method and adds the  
 * IObjectiveFunctionParallel interface to a class that implements 
 * the interface IObjectiveFunction 
 * see CMAEvolutionStrategy.java for copyright.
 */
public abstract class AbstractObjectiveFunction implements 
IObjectiveFunction,
IObjectiveFunctionParallel  { 

    /**
     *
     @param x
     @return
     */
    abstract public double valueOf(double[] x);
    public double [] valuesOf(double[][] pop) {
        double [] res = new double[pop.length];
        for (int i = 0; i < pop.length; ++i)
            res[i] = valueOf(pop[i]);
        return res;
    }

    /**
     *
     @param x
     @return
     */
    public boolean isFeasible(double[] x) {
    	return true;
    }
}
