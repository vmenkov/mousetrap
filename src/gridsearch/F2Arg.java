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
    };

    Res optimizeOverOneVar(ParVec fixedPar, LookFor lookFor, int minOver) {

	final int mfactor = 3; //10;
	Grid g = Grid.cubeGrid(fixedPar.dim(), mfactor);
	final int maxlevel = 4;
	
	return optimizeOverOneVarLoop(fixedPar, g, lookFor, minOver, mfactor, maxlevel);
	
    }

    static boolean debug = true;
    /** When a finer grid is created, how many cells of the coarser grid, in
	each direction, does it create? */
    final int L=1;
    
    /** @param minOver which param one varies? 0 alpha, 1 beta. The other is fixed.
	@param lookFor Does "optimize" mean "minimize" or "maximize"? 
     */
    Res optimizeOverOneVarLoop(ParVec fixedPar, Grid g, LookFor lookFor, int minOver, int mfactor, int maxlevel) {

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

}
 
