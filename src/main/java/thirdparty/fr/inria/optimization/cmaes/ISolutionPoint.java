package thirdparty.fr.inria.optimization.cmaes;

    /** solution point in search space, single-objective case
     *  
     * see copyright in CMAEvolutionStrategy
     * 
     * */
public interface ISolutionPoint {
    /** objective function value (fitness) of the search point x
     @return  */ 
    public double getFitness();  
    /** count at what evaluation number the search point x was evaluated
     @return  */
    public long getEvaluationNumber();
    /** value of the point in search space, that is in the 
     * preimage of the objective function to be optimized 
     @return  */ 
    public double[] getX();
    
    /** objective function value (fitness) of the search point x
     @param fitness */ 
    public void setFitness( double fitness); // TODO better FunctionValue than Fitness ? 
    /** count at what evaluation number the search point x was evaluated
     @param evaluation */
    public void setEvaluationNumber( long evaluation);
    /** value of the solution point in search space, the 
     * preimage of the objective function to be optimized 
     @param x */ 
    public void setX(double[] x);
}
