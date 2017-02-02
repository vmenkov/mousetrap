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

	final int mfactor = 10;
	Grid g = Grid.cubeGrid(fixedPar.dim(), mfactor);
	final int maxlevel = 4;
	
	return optimizeOverOneVarRecursive(fixedPar, g, lookFor, minOver, mfactor, 0, maxlevel);
	
    }

    static boolean debug = true;
    
    /** @param minOver which param one varies? 0 alpha, 1 beta. The other is fixed.
	@param lookFor Does "optimize" mean "minimize" or "maximize"? 
     */
    Res optimizeOverOneVarRecursive(ParVec fixedPar, Grid g, LookFor lookFor, int minOver, int mfactor, int level, int maxlevel) {

	Res best = null;
	for(Iterator<ParVec> it = g.getParVecIterator(); it.hasNext(); ){
	    ParVec q = it.next();
	    ParVec[] args = new ParVec[2];
	    args[ minOver ] = q;
	    args[ 1-minOver ] = fixedPar;
	    double val = f(args);
	    if (best == null ||
		(lookFor.min()? val<best.val : val>best.val)) best=new Res(args,val);
	}
	if (debug) System.out.println("At level=" + level + ", " +
				      lookFor + " at " + best);
				      
	
	if (level == maxlevel) return best;
	final int L=1;
	Grid g1 = g.vicinityGrid(best.ab[minOver], mfactor, L);

	return  optimizeOverOneVarRecursive( fixedPar, g1, lookFor,  minOver, mfactor, level+1,  maxlevel);

    }

}
 
