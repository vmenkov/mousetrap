package gridsearch;

import java.util.*;

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


    static boolean debug = true;
    /** When a finer grid is created, how many cells of the coarser grid, in
	each direction, does it create? */
    static final int L=1;

    /** Optimization (min or max) over one variable 
	@param fixedPar This variable stays constant
	@param minOver This is the variable we optimize over
	@param lookFor Look for min or max?
	@return The values of arguments at which the min or max is reached, and the function value at them.
     */
    Res optimizeOverOneVar(boolean simplex, ParVec fixedPar, LookFor lookFor, int minOver) {
	final int mfactor = 3; //10;
	final int maxlevel = 4;
	Grid g = simplex? Grid.simplexGrid(fixedPar.dim(), mfactor) :
	    Grid.cubeGrid(fixedPar.dim(), mfactor);
	return optimizeOverOneVarLoop(fixedPar, g, lookFor, minOver, mfactor, maxlevel);	
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

    /** @param outerLookFor:  the outer optimization is min or max
	@param outerMinOver: the outer optimization is for variable 0 or 1
     */
    Res findSaddlePoint( boolean simplex, int dim, LookFor outerLookFor, int outerMinOver) {

	final int mfactor = 3; //10;
	final int maxlevel = 4;

	Grid og = simplex?
	    Grid.simplexGrid(dim, mfactor):
	    Grid.cubeGrid(dim, mfactor);

	for(int level = 0; ; level++) {

	    Res best = null;
	    for(Iterator<ParVec> it = og.getParVecIterator(); it.hasNext(); ){
		ParVec fixedPar = it.next();
		Res r = optimizeOverOneVar(simplex, fixedPar, outerLookFor.other(), 1-outerMinOver);
		
		//if (debug) System.out.println(r);
		if (best == null ||
		    (outerLookFor.min()? r.val<best.val: r.val>best.val)) best=r;
	    }
	    if (debug) System.out.println("Outer level=" + level + ", " +
				      outerLookFor + " at " + best);
				      	
	    if (level == maxlevel) return best;
	    og = og.vicinityGrid(best.ab[outerMinOver], mfactor, L);
	}
    }

    

}
 
