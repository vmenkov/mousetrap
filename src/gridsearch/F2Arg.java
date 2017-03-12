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
	@param hints If not null, this parameter tells us that sometimes there is no need to optimize beyond a certain point in the inner loop of max_a min_b. The inner loop can return at once if it knows that its "winner" won't be as good at the winner-so-far of the outer loop.
	@return The values of arguments at which the min or max is reached, and the function value at them.
     */
    Res optimizeOverOneVar(ParVec fixedPar, Constraint cons, int dim, LookFor lookFor, int minOver, Hints hints) {
	Grid g = Grid.cubeGrid(dim, params.mfactor, cons);
	return optimizeOverOneVarLoop(fixedPar, g, lookFor, minOver, params.mfactor, params.maxlevel, hints);	
    }

    /** @param minOver which param one varies? 0 alpha, 1 beta. The other is fixed.
	@param lookFor Does "optimize" mean "minimize" or "maximize"? 
     */
    private Res optimizeOverOneVarLoop(ParVec fixedPar, Grid g, LookFor lookFor, int minOver, int mfactor, int maxlevel, Hints hints) {

	Res best = null;

	if (hints!=null && hints.bestSaddleSoFar!=null) {
	    ParVec[] args = new ParVec[2];
	    args[ 1-minOver ] = fixedPar;
	    args[ minOver ] = hints.bestSaddleSoFar.ab[minOver];
	    double val = f(args);
	    best=new Res(args,val);
	    if (hints.willNotWin(val)) return best;
	}


	for(int level = 0; ; level++) {

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
		    (lookFor.min()? val<best.val : val>best.val)) {
		    best=new Res(args,val);
		    if (hints!=null && hints.willNotWin(val)) return best;
		}
	       

	    }
	    if (debug) System.out.println("At level=" + level + ", " +
					  lookFor + " at " + best);
				      	
	    if (level == maxlevel) return best;
	    g = g.vicinityGrid(best.ab[minOver], mfactor, L);
	}
    }

    /** This is used to short-circuit computations in the inner loop
	(e.g.  "Min_b" in "Max_a Min_b", if we know that the result of
	the inner loop ("Min_b") won't be as good as the best outer loop
	result (the "Max_a" so far).
     */
    static class Hints {
	Res bestSaddleSoFar;
	LookFor innerLookFor;
	Hints(Res r, LookFor _innerLookFor) {
	    bestSaddleSoFar = r;
	    innerLookFor = _innerLookFor;
	}
	/** This inner loop (say, for min_b) won't beat the current max_a min_b,
	    because we know that the current min_b will be smaller than the
	    currently found max_a
	 */
	boolean willNotWin(double f) {
	    if (bestSaddleSoFar==null) return false;
	    return (innerLookFor==LookFor.MIN) ? 
		f <= bestSaddleSoFar.val : f >= bestSaddleSoFar.val;
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
	Res best = null;
	final LookFor innerLookFor =  outerLookFor.other();
		
	for(int level = 0; ; level++) {

	    for(Iterator<ParVec> it = og.getParVecIterator(); it.hasNext(); ){
		ParVec fixedPar = it.next();
	
		Hints hints = new Hints(best, innerLookFor);
		Res r = optimizeOverOneVar(fixedPar, cons[inner], dim[inner], innerLookFor, inner, hints);
		
		//if (debug) System.out.println(r);
		if (best == null ||
		    (outerLookFor.min()?r.val<best.val: r.val>best.val)) best=r;
	    }
	    if (debug) System.out.println("Outer level=" + level + ", " +
				      outerLookFor + " at " + best);
				      	
	    if (level == params.maxlevel) return best;
	    og = og.vicinityGrid(best.ab[outerMinOver], params.mfactor, L);
	}
    }

    

}
 
