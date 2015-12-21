package mousetrap;

import java.util.*;
//import java.text.*;
//import java.io.*;

class SimplexResults {
    double [] p;
    double maxval=0;

    SimplexResults(double[][] payoffMatrix) {
	final int L=payoffMatrix.length;
	final int H=payoffMatrix[0].length;
	final int L1 = L-1;
	int ibest = 0;
	Rational[] xbest = null;
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


