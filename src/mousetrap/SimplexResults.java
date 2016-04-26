package mousetrap;

import java.util.*;
//import java.text.*;
//import java.io.*;

/** Our interface to the Simplex algorithm code. */
class SimplexResults {
    /** Will contain the result of the optimization, vector p.
     */
    double [] p;
    /** Will contain f(p) for the p found by optimization.
    */
    double maxval=0;

    /** For an L-by-H matrix A, finds the L-dimensional vector p such
	that f(p)=max<sub>p' in P</sub>f(p),
 where P is the set of
	vectors whose components are non-negative and sum to 1
	(i.e. e<sup>T</sup>p=1).

	<center>
	f(p) = min<sub>q</sub> p<sup>T</sup> A q.
	</center>

	<p>This is done by means of describing the problem in the form
	<center>
	f(p) = min<sub>j</sub> &Sigma;<sub>j</sub> a<sub>ij</sub> p<sub>i</sub>	
	</center>
	For each <em>j</em>, we consider the set P<sub>j</sub>, defined as the
	subset of P on which the minimum in the above formula is reached at that
	j. 

	@param payoffMatrix Matrix A.
     */
    SimplexResults(double[][] payoffMatrix) {
	final int L=payoffMatrix.length;
	final int H=payoffMatrix[0].length;
	final int L1 = L-1;
	int ibest = 0;
	Rational[] xbest = null;
	// maximization within each of H subsets
	for(int i=0; i<H; i++) { 
	    Rational[][] A = new Rational[H][];
	    Rational[] b = new Rational[H], c= new Rational[L1], x = new Rational[L1];
	    
	    for(int j=0; j<H; j++) {
		A[j] = new Rational[L1];
		if (j==i) {
		    // the last inequality
		    for(int k=0; k<L1; k++) {
			A[j][k] = Rational.ONE;
		    }
		    b[j] = Rational.ONE;
		} else {
		    double s = payoffMatrix[L1][i] - payoffMatrix[L1][j];
		    for(int k=0; k<L1; k++) {
			A[j][k] = Mousetrap.approxRational( payoffMatrix[k][i] - payoffMatrix[k][j] - s);
		    }
		    b[j] =Mousetrap.approxRational( -s );
		}
	    }
	    for(int k=0; k<L1; k++) {
		c[k] = Mousetrap.approxRational( payoffMatrix[k][i] - payoffMatrix[L1][i]);
	    }
	    // omega = max_{x : such that A*x <= b, x >= 0 }( c*x)
	    Rational omega = Simplex.simplex(A, b, c, x);
	    if (omega==null) {
		// The feasible space is empty
		continue;		
		//throw new IllegalArgumentException("Simplex.simplex returned null!");
	    }
	    double val = omega.doubleValue() + payoffMatrix[L1][i];
	    if (xbest==null || val>maxval) {
		maxval = val;
		ibest = i;
		xbest = x;
	    }
	}
	if (xbest==null) throw new IllegalArgumentException("All simplexes were empty!?");
	p = new double[L];
	double s = 0;
	for(int k=0; k<L1; k++) {
	    s += p[k] = xbest[k].doubleValue();
	}
	p[L1] = 1 - s;
   }

}


