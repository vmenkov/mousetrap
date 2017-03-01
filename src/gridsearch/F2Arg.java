package gridsearch;

import java.util.*;
import util.*;

/** Class used for grid-based minimax optimization of an arbitrary
    function of two arguments. The arguments are of the ParVec type.
*/
abstract class F2Arg {
    abstract double f(ParVec alpha, ParVec beta);

    double f(ParVec[] ab) { return f(ab[0], ab[1]); }

    protected static class Res {
	ParVec[] ab;
	double val;
	Res(ParVec[] _ab, double _val) { ab = _ab; val = _val; }
	public String toString() {
	    return "{alpha=" + ab[0] + "; beta=" + ab[1] + "; f=" + val+"}";
	}
	
    }

    enum LookFor {
	MIN, MAX;
	boolean min() { return this==MIN; }
	LookFor other()  { return min()? MAX: MIN; }
    };

    /** Parameters controlling the grid search computational optimization. The default values can be overridden by initParams()
     */
    static class Parameters {
	int mfactor = 3; //10;
	int maxlevel = 4;
	public String toString() {
	    return "Gridsearch parameters: range is divided into mfactor=" + mfactor + " sections at each level; total of " + maxlevel + " levels";
	}
    }

    static Parameters params = new Parameters();

    /** Sets parameter values from command line options.
	@param ht Represents command-line options (-Dname=value)
    */
    static void initParams(ParseConfig ht) {
	params.mfactor = ht.getOption("grid.mfactor", params.mfactor);
	params.maxlevel = ht.getOption("grid.maxlevel", params.maxlevel);	
    }

    static boolean debug = false;//true;
    /** When a finer grid is created, how many cells of the coarser grid, in
	each direction, does it create? */
    static final int L=1;

    /** Optimization (min or max) over one variable 
	@param fixedPar This variable stays constant
	@param minOver This is the variable we optimize over
	@param lookFor Look for min or max?
	@return The values of arguments at which the min or max is reached, and the function value at them.
     */
    Res optimizeOverOneVar(ParVec fixedPar, Constraint cons, int dim, LookFor lookFor, int minOver) {
	Grid g = Grid.cubeGrid(dim, params.mfactor, cons);
	return optimizeOverOneVarLoop(fixedPar, g, lookFor, minOver, params.mfactor, params.maxlevel);	
    }

    /** @param minOver which param one varies? 0 alpha, 1 beta. The other is fixed.
	@param lookFor Does "optimize" mean "minimize" or "maximize"? 
     */
    private Res optimizeOverOneVarLoop(ParVec fixedPar, Grid g, LookFor lookFor, int minOver, int mfactor, int maxlevel) {

	for(int level = 0; ; level++) {

	    Res best = null;
	    for(Iterator<ParVec> it = g.getParVecIterator(); it.hasNext(); ){

		ParVec[] args = new ParVec[2];
		args[ 1-minOver ] = fixedPar;
		args[ minOver ] = it.next();

		if (debug && !g.constraint.holds(args[ minOver ])) {
		    System.out.println("Error context:\ngrid=\n" +g +
				       "\nit.ci=" + ((Grid.ParVecIterator)it).ci + 
				       "\npoint=" + args[ minOver ]);
		    throw new IllegalArgumentException("Iterator over the grid produced a point value that does not satisfy the constraint");
		}


		double val = f(args);
		//if (debug) System.out.println("f("+args[minOver]+")=" + val);
		if (best == null ||
		    (lookFor.min()? val<best.val : val>best.val)) best=new Res(args,val);
	    }
	    	    if (debug) System.out.println("At level=" + level + ", " +
	    				      lookFor + " at " + best);
				      	
	    if (level == maxlevel) return best;
	    g = g.vicinityGrid(best.ab[minOver], mfactor, L);
	}
    }

    /** Finds the min max or max min of a function of 2 arguments.
	@param dim (for the two arguments of f[])
	@param cons Geometrical constraints for the two domains
	@param outerLookFor: the outer optimization is min or max
	@param outerMinOver: the outer optimization is for variable 0 or 1 (i.e. min_{x[0]} max_{x[1]} or max_{x[1]} min_{x[0]} 
     */
    Res findSaddlePoint(int[] dim, Constraint cons[], 
			LookFor outerLookFor, int outerMinOver) {       
	if (cons==null) cons = new Constraint[2];

	Grid og = Grid.cubeGrid(dim[outerMinOver], params.mfactor, cons[outerMinOver]);
	final int inner  = 1 - outerMinOver;

	for(int level = 0; ; level++) {

	    Res best = null;
	    for(Iterator<ParVec> it = og.getParVecIterator(); it.hasNext(); ){
		ParVec fixedPar = it.next();
	
		Res r = optimizeOverOneVar(fixedPar, cons[inner], dim[inner], outerLookFor.other(), inner);
		
		//if (debug) System.out.println(r);
		if (best == null ||
		    (outerLookFor.min()? r.val<best.val: r.val>best.val)) best=r;
	    }
	    if (debug) System.out.println("Outer level=" + level + ", " +
				      outerLookFor + " at " + best);
				      	
	    if (level == params.maxlevel) return best;
	    og = og.vicinityGrid(best.ab[outerMinOver], params.mfactor, L);
	}
    }

    

}
 
