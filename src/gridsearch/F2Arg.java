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
    Res optimizeOverOneVar(ParVec fixedPar, Constraint cons, int dim, LookFor lookFor, int minOver) {
	final int mfactor = 3; //10;
	final int maxlevel = 4;
	//	Constraint cons=simplex? SingleConstraint.simplex(fixedPar.dim()) :null;
	Grid g = Grid.cubeGrid(dim, mfactor, cons);
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

		if (!g.constraint.holds(args[ minOver ])) {
		    System.out.println("Error context:\ngrid=\n" +g +"\npoint=" + args[ minOver ]);
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

    /** @param dimnension (for the two arguments of f[])
	@param outerLookFor:  the outer optimization is min or max
	@param outerMinOver: the outer optimization is for variable 0 or 1
     */
    Res findSaddlePoint( //boolean simplex, 
			int[] dim, Constraint cons[], 
			LookFor outerLookFor, int outerMinOver) {       
	if (cons==null) cons = new Constraint[2];

	final int mfactor = 3; //10;
	final int maxlevel = 4;

	//Constraint cons=simplex? SingleConstraint.simplex(dim) :null;

	Grid og = Grid.cubeGrid(dim[outerMinOver], mfactor, cons[outerMinOver]);
	final int inner  = 1 - outerMinOver;

	for(int level = 0; ; level++) {

	    Res best = null;
	    for(Iterator<ParVec> it = og.getParVecIterator(); it.hasNext(); ){
		ParVec fixedPar = it.next();

		//		Constraint cons2=simplex? SingleConstraint.simplex(fixedPar.dim()) :null;

		Res r = optimizeOverOneVar(fixedPar, cons[inner], dim[inner], outerLookFor.other(), inner);
		
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
 
