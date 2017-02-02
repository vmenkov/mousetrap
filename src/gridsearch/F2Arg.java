package gridsearch;

import java.util.*;

/** Class used for grid-based minimax optimization of an arbitrary
    function of two arguments. The arguments are of the ParVec type.
*/
abstract class F2Arg {
    abstract double f(ParVec alpha, ParVec beta);

    double f(ParVec[] ab) { return f(ab[0], ab[1]); }

    private static class Res {
	ParVec[] ab;
	double val;
	Res(ParVec[] _ab, double _val) { ab = _ab; val = _val; }
    }


    Res minOverOneVar(ParVec fixedPar, boolean lookForMin, int minOver) {

	final int mfactor = 10;
	Grid g = Grid.cubeGrid(fixedPar.dim(), mfactor);
	final int maxlevel = 4;
	
	return optimizeOverOneVarRecursive(fixedPar, g, lookForMin, minOver, mfactor, 0, maxlevel);
	
    }
    
    /** @param minOver which param one varies? 0 alpha, 1 beta. The other is fixed.
	@param lookForMin Does "optimize" mean "minimize" or "maximize"? If this is true, we are minimizing f; otherwise, we're maximizing f
     */
    Res optimizeOverOneVarRecursive(ParVec fixedPar, Grid g, boolean lookForMin, int minOver, int mfactor, int level, int maxlevel) {

	Res best = null;
	for(Iterator<ParVec> it = g.getParVecIterator(); it.hasNext(); ){
	    ParVec q = it.next();
	    ParVec[] args = new ParVec[2];
	    args[ minOver ] = q;
	    args[ 1-minOver ] = fixedPar;
	    double val = f(args);
	    if (best == null ||
		lookForMin? val<best.val : val>best.val) best=new Res(args,val);
	}
	if (level == maxlevel) return best;
	final int L=1;
	Grid g1 = g.vicinityGrid(best.ab[minOver], mfactor, L);

	return  optimizeOverOneVarRecursive( fixedPar, g1, lookForMin,  minOver, mfactor, level+1,  maxlevel);

    }
}
 
