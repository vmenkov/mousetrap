package gridsearch;

import java.util.*;
import mousetrap.*;

/** Tools to convert a parameter vector (from the space over which we optimize)
    to an actual matrix */
public class ParametrizedMatrix  {

    /** Describes some form of symmetry of the graph (a
	homomorphism). When optimizing for the players'
	strategies, we only look for strategies that are
	similarly symmetric (i.e. invariant with respect to
	the graph's mapping onto itself by the homomorphism).
     */
    static class Symmetry {
	static final int NONE = -1, REST= -2;
	int mapsto[];
	/** Which other value in aptr[][]   corresponds to aptr[k][j]? */
	int lookup(int w[][], int aptr [][], int k, int j) {
	    // edge (k,m) is mapped to (k1,m1) by the symmetry
	    int m = w[k][j];
	    int k1 = mapsto[k];
	    int m1 = mapsto[m];

	    if (k1==NONE || m1==NONE) return NONE;

	    for(int i=0; i<w[k1].length; i++) {
		if (w[k1][i] == m1) {
		    if (aptr[k1]==null) return NONE;
		    int a1 = aptr[k1][i];
		    return  (a1 == NONE) ? NONE : a1;
		}
	    }
	    return NONE;
	}

	static Symmetry none(int n) {
	    Symmetry s = new Symmetry();
	    s.mapsto = new int[n];
	    for(int i=0; i<n; i++) s.mapsto[i] = NONE;
	    return s;
	}
      	
	static Symmetry mirror(int n) {
	    int z[] = new int[n];
	    Symmetry s = new Symmetry();
	    for(int i=0; i<n; i++) s.mapsto[i] = n-1-i;
	    return s;	    
	}

	static int indexOf(int[] a, int z) {
	    for(int i=0; i<a.length; i++) {
		if (a[i]==z) return i;
	    }
	    return -1;
	}

	void verify(int[][] w) {
	    final String msg = "The Symmetry map is not consistent with the graph structure";
	    for(int k=0; k<w.length; k++) {
		int k1 = mapsto[k];
		if (k1==NONE || k1==k) continue;
		if (w[k].length != w[k1].length) throw new IllegalArgumentException(msg + ": w["+k+"].length="+w[k].length+", w["+k1+"].length="+w[k1].length);
		for(int i=0; i<w[k].length; i++) {
		    int r0 = w[k][i];
		    int r = mapsto[r0];
		    if (indexOf(w[k1],r) < 0) throw new IllegalArgumentException(msg + ": w["+k+"] has "+r0+", w["+k1+"] has no "+r);
		}
	    }
	}

    }

    private static class AsgMap {
	int [][] aptr;
	int apos;

	/** Fills aptr[][] */
	AsgMap(int w[][], Symmetry sym, int apos0) {
	    apos = apos0;
	    aptr =   arrayStructureCopy(w,  Symmetry.NONE);
	    for(int k=0; k<w.length; k++) {
		for(int i=0; i<w[k].length; i++) {
		    int a1 = sym.lookup( w,  aptr, k, i);
		    if (w[k][i] == k) {
			// the diagonal element = 1 - sum(others)
			aptr[k][i] = Symmetry.REST;
		    } else if (a1 != Symmetry.NONE)  {		    // use symmetry...
			aptr[k][i] = a1;
		    } else {
			aptr[k][i] = apos ++;
		    }
		}
	    }
	}
    }

    final int [][] w;
    int nvar;
    int [][] aposUnseen, aposSeen;
    

    ParametrizedMatrix(int _w[][], Symmetry sym) {
	w = _w;
	AsgMap map1 = new AsgMap(w, sym, 0);
	AsgMap map2 = new AsgMap(w, sym, map1.apos);
	aposUnseen = map1.aptr;
	aposSeen = map2.aptr;
	nvar = map2.apos;
    }

    /** An object of this class describes an actual transition matrix for 
	one player */
    static class MatrixData {
	final int[][] w;
	double [][]  aUnseen, aSeen;

	private static double[][] fillData(final int apos[][], double [] q) {
	    double [][] a = new double[apos.length][];
	    for(int k=0; k<apos.length; k++) {
		double s = 0;
		int diagPos = Symmetry.NONE;
		a[k] = new double[apos[k].length];
		for(int i=0; i<apos[k].length; i++) {
		    int p = apos[k][i];
		    if (p == Symmetry.REST) {
			if (diagPos!=Symmetry.NONE) throw new  IllegalArgumentException("Two diagonal values?!");
			diagPos=i;
		    } else {
			a[k][i] = q[p];
			s += 	a[k][i];
		    }
		}
		if (diagPos==Symmetry.NONE) throw new  IllegalArgumentException("No diagonal values found for k="+k);
		a[k][diagPos] = 1.0 - s;
	    }
	    return a;
	}

	MatrixData(ParametrizedMatrix mi, double [] q) {
	    w = mi.w;
	    if (q.length != mi.nvar) throw new IllegalArgumentException("var cnt mismatch");
	    aUnseen =  fillData(mi.aposUnseen, q);
	    aSeen   =  fillData(mi.aposSeen, q);
	}
    }

    /** A JointProbVector describes the probability distribution over
	all possible states of the two-player system.  xUnseen[i][j]
	is the probability of Player A being at i and Player B at j,
	and them not seeing each other; xSeen[i] is the probability of
	Players A and B both being at i, and seeing each other. All
	probabilities should sum to 1.0.
     */
    static class JointProbVector {
	double [][] xUnseen;
	double [] xSeen;
	int n() { return xSeen.length; }
	JointProbVector(int n) {
	    xUnseen = zeroMat(n,n);
	    xSeen = new double[n];
	}

	/** Checks if the values sum to 1.0 (within a computational error) 
	 */
	void validate() {
	    double ss=sumSeen(), su=0;
	    for(int i=0; i< xUnseen.length; i++) {
		for(int j=0; j< xUnseen[i].length; j++) {
		    su += xUnseen[i][j];
		}
	    }
	    double s=ss+su;
	    if (Math.abs(s - 1.0) > 1e-6) throw new IllegalArgumentException("Probabilities don't sum to 1.0. ss=" + ss+", su=" + su+", s=" + s);
	}


	double sumSeen() {
	    double ss=0;
	    for(int i=0; i< xSeen.length; i++) {
		ss += xSeen[i];
	    }
	    return ss;
	}


	static double[][] zeroMat(int nrow, int ncol) {
	    double [][] x = new double[nrow][];
	    for(int i=0; i<x.length; i++)  x[i] = new double[ncol];
	    return x;
	}

	/**
	   xi'_i = phi * ( \sum_k aSeen_{ik} bSeen_{ik} * xi_k  +
	       \sum_{kl} aUnseen_{ik} bUnseen_{ik} * x_{kl} ).

	   Remember that sparse matrices A and B are stored by column.
	 */
	JointProbVector apply( MatrixData a, MatrixData b, double phi) {
	    JointProbVector res = new JointProbVector(n()); 

	    for(int k=0; k<a.w.length; k++) {
		for(int pi=0; pi<a.w[k].length; pi++) {
		    int i = a.w[k][pi];

		    for(int pj=0; pj<b.w[k].length; pj++) {
			int j = b.w[k][pj];
			
			double r = a.aSeen[k][i] * b.aSeen[k][j] * xSeen[k];
			res.xUnseen[i][j] += r;
			if (i==j) res.xSeen[i] += r;
		    }
		}
	    }

	    double [][] u = zeroMat(b.w.length, b.w.length);

	    for(int k=0; k<a.w.length; k++) {
		for(int l=0; l<b.w.length; l++) {
		    for(int pj=0; pj<b.w[l].length; pj++) {
			int j = b.w[l][pj];
			u[j][k] += b.aUnseen[l][j] * xUnseen[k][l];
		    }
		}
	    }

	    for(int j=0; j<b.w.length; j++) {
		for(int k=0; k<a.w.length; k++) {
		    for(int pi=0; pi<a.w[k].length; pi++) {
			int i = a.w[k][pi];
			double r = a.aSeen[k][pi] * u[j][k];
			res.xUnseen[i][j] += r;
			if (i==j) res.xSeen[i] += r;
		    }
		}
	    }

	    for(int i=0; i<a.w.length; i++) {
		res.xSeen[i] *= phi;
		res.xUnseen[i][i] -= 	res.xSeen[i];
	    }
	
	    return res;
	}


    }
    

    /** Creates a new 2-dimensional array with the same structure as
	w[][], filled with a specified value (such as NONE) */
    static int[][] arrayStructureCopy(int w[][], int val) {
	int [][] aptr = new int[w.length][];
	for(int k=0; k<w.length; k++) {
	    aptr[k] = new int[w[k].length];
	    for(int i=0; i<w[k].length; i++) {
		aptr[k][i] = val;
	    }
	}
 	return aptr;
    }

    static String a2str(int a) {
	return a==Symmetry.NONE ? "?" :
	    a==Symmetry.REST ? "R" :  "" + a;
    }

    static String report(int w[][], int aptr[][], Symmetry sym) {
	StringBuffer b = new StringBuffer();
	for(int k=0; k<w.length; k++) {
	    b.append("[");
	    for(int i=0; i<w[k].length; i++) {
		b.append(" " + w[k][i] + ":" + a2str(aptr[k][i]));
	    }
	    b.append("]\n");
	}
	return b.toString();
    }

    public class F2ArgImmediatePayoff extends F2Arg {
  	ParametrizedMatrix aScheme, bScheme;
	JointProbVector jpv0;
	final double phi;
	F2ArgImmediatePayoff(Mousetrap2 mo, Symmetry sym, 
			     JointProbVector _jpv0) {
	    phi = mo.phi;
	    jpv0 = _jpv0;
	    aScheme = new ParametrizedMatrix(mo.w, sym);
	    bScheme = new ParametrizedMatrix(mo.w2, sym);
	}

	double f(ParVec alpha, ParVec beta) {
	    MatrixData amat = new MatrixData(aScheme, alpha.getX());
	    MatrixData bmat = new MatrixData(bScheme, beta.getX());
	    JointProbVector jpv = jpv0.apply(amat, bmat, phi);
	    return jpv.sumSeen();
	}
    }

    static public void main() {
	Mousetrap2 mo = Mousetrap2.mo1();
	Symmetry sym = Symmetry.mirror(mo.h);
	JointProbVector jpv = new JointProbVector(mo.h);
	int startPoint = mo.h/2;
	jpv.xUnseen[startPoint][startPoint] = 1.0;  // both players start at the center point



    }


}
