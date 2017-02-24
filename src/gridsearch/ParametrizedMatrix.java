package gridsearch;

import java.io.*;
import java.util.*;
import mousetrap.*;

/** Tools to convert a parameter vector (from the space over which we optimize)
    to an actual matrix */
public class ParametrizedMatrix  {

    private static class AsgMap {
	int [][] aptr;
	int apos;

	/** Fills aptr[][]; adds constraints to mc */
	AsgMap(int w[][], Symmetry sym, int apos0, MultiConstraint mc) {
	    apos = apos0;
	    aptr =   arrayStructureCopy(w,  Symmetry.NONE);
	    for(int k=0; k<w.length; k++) {
		HashSet<Integer> h = new HashSet<Integer>();
		for(int i=0; i<w[k].length; i++) {
		    int a1 = sym.lookup( w,  aptr, k, i);
		    if (w[k][i] == k) {
			// the diagonal element = 1 - sum(others)
			aptr[k][i] = Symmetry.REST;
		    } else if (a1 != Symmetry.NONE)  {	    // use symmetry...
			aptr[k][i] = a1;
			h.add(aptr[k][i] );
		    } else {
			aptr[k][i] = apos ++;
			h.add(aptr[k][i] );
		    }
		}		
		mc.addSimplexConstraintIfUnique(h);

	    }
	}
    }

    /** The matrix structure (same as in Mousetrap class) */
    final int [][] w;
    /** The number of parameters */
    int nvar;
    /** How matrix elements are based on parameters. The structure of
	these arrays is the same as w[][]. */
    int [][] aposUnseen, aposSeen;
    /** Constraints used to restrict the space of legal parameter combinations */
    Constraint constraint;

    ParametrizedMatrix(int _w[][], Symmetry sym) {
	w = _w;
	MultiConstraint mc = new MultiConstraint();

	System.out.println("Building map for unseen, start apos=0");
	AsgMap map1 = new AsgMap(w, sym, 0, mc);
	System.out.println("Building map for unseen, start apos=" + map1.apos);
	AsgMap map2 = new AsgMap(w, sym, map1.apos, mc);
	aposUnseen = map1.aptr;
	aposSeen   = map2.aptr;
	nvar = map2.apos;
	constraint = mc;
    }

    /** An object of this class describes an actual transition matrix for 
	one player */
    static class MatrixData {
	final int[][] w;
	/** Sparse matrices whose structure (by column) is described by w[][] */
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
		if (diagPos==Symmetry.NONE) throw new  IllegalArgumentException("No diagonal value found for k="+k);
		a[k][diagPos] = 1.0 - s;
		if (a[k][diagPos]<0) {
		    throw new  IllegalArgumentException("Negative diagonal value ("+a[k][diagPos]+") computed for k="+k);
		}
	    }
	    return a;
	}


	MatrixData(ParametrizedMatrix mi, ParVec qv) {
	    this(mi, qv.getX());
	}

	MatrixData(ParametrizedMatrix mi, double [] q) {
	    w = mi.w;
	    if (q.length != mi.nvar) throw new IllegalArgumentException("var cnt mismatch");
	    aUnseen =  fillData(mi.aposUnseen, q);
	    aSeen   =  fillData(mi.aposSeen, q);
	}

	/** @param s aUnseen[][] or aSeen[][] */
	double [][] toDenseMatrix(double [][] s) {
	    double[][] a = JointProbVector.zeroMat(w.length, w.length);
	    for(int k=0; k<w.length; k++) {
		for(int i=0; i<w[k].length; i++) {
		    a[k][w[k][i]] = s[k][i];
		}
	    }
	    return a;
	}

	double [][][] toDenseMatrices() {
	    return new double [][][] { toDenseMatrix(aUnseen),
				       toDenseMatrix(aSeen) };
	}

	public String toString() {
	    StringBuffer b = new StringBuffer();
	    double [][][] z = toDenseMatrices();
	    for(int k=0; k<2;k++) {
		b.append(k==0? "Unseen\n": "Seen\n");
		for(double[] q: z[k]) {
		    b.append("[" + Arrays.toString(q) + "]\n");
		}
	    }
	    return b.toString();
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
	    a==Symmetry.REST ? "X" :  "" + a;
    }

    /** Displays the matrix structure and the way matrix elements are
	based on the parameters */
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

    public static class F2ArgPayoff extends F2Arg {
  	ParametrizedMatrix aScheme, bScheme;
	JointProbVector jpv0;
	final double phi;

  	ParametrizedMatrix[] schemes() {
	    return new  ParametrizedMatrix[] { aScheme, bScheme};
	}

	final boolean ev;

	F2ArgPayoff(Mousetrap2 mo, Symmetry sym, 
		    JointProbVector _jpv0, boolean _ev) {
	    ev = _ev;
	    phi = mo.phi;
	    jpv0 = _jpv0;
	    jpv0.validate();
	    aScheme = new ParametrizedMatrix(mo.w, sym);
	    bScheme = new ParametrizedMatrix(mo.w2, sym);

	    System.out.println("Player A parametrization map:");
	    System.out.println( report(mo.w, aScheme.aposUnseen, sym));

	    System.out.println("Player B parametrization map:");
	    System.out.println( report(mo.w2, bScheme.aposUnseen, sym));

	}

	double f(ParVec alpha, ParVec beta) {
	    return ev? f_longTerm(alpha,beta) :  
		f_ImmediatePayoff(alpha,beta);
	}

	double f_ImmediatePayoff(ParVec alpha, ParVec beta) {
	    MatrixData amat = new MatrixData(aScheme, alpha);
	    MatrixData bmat = new MatrixData(bScheme, beta);
	    JointProbVector jpv = jpv0.apply(amat, bmat, phi);
	    return jpv.sumSeen();
	}

	int callsTotal=0, callsNoConv=0;
	long tTotal=0;

	String statsReport() {
	    return ev?
		"EV: " + callsTotal + " calls so far; failed to achieve convergence in " +  callsNoConv + " calls; <T> = "+((double)tTotal)/callsTotal:
		"Imediate payoff: " + callsTotal + " calls so far";
	}

	double f_longTerm(ParVec alpha, ParVec beta) {
	    if (callsTotal % 10000 ==0) {
		System.out.println(statsReport());
	    }


	    callsTotal++;
	    MatrixData amat = new MatrixData(aScheme, alpha);
	    MatrixData bmat = new MatrixData(bScheme, beta);

	    final int maxT = 100;
	    int t=0;
	    double r, r0=0;
	    JointProbVector jpv = jpv0;
	    while(true) {
		t++;
		jpv = jpv.apply(amat, bmat, phi);		
		r = jpv.sumSeen() / (jpv.sumSeen() + jpv.sumUnseen());
		if (Double.isInfinite(r)) {
		    System.out.println("Infinity encountered: t=" + t+", jpv.sumSeen() = " + jpv.sumSeen() + ", jpv.sumUnseen()=" + jpv.sumUnseen());
		    System.out.println("jpv=" + jpv);
		    System.out.println("amat=" + amat);
		    System.out.println("bmat=" + bmat);
		    
		    throw new IllegalArgumentException();
		}
		if (t>1) {
		    if (Math.abs(r-r0) < 1e-4) break;
		    r0 = r;
		}
		if (t>maxT) {
		    callsNoConv ++;
		    break;
		}
	    }
	    tTotal += t;
	    return r;

	}


    }

    /** The function being optimized is the payoff for the defender (i.e. the
	second player). Thus we go for min_A max_B  and max_B min_A
    */
    static public void main(String [] argv) {
	Mousetrap2 mo = Mousetrap2.mo1();
	//Mousetrap2 mo = Mousetrap2.moChain(4, false, 1, 1);


	Symmetry sym = Symmetry.mirror(mo.h);
	JointProbVector jpv = new JointProbVector(mo.h);
	int startPoint = mo.h/2;
	jpv.xUnseen[startPoint][startPoint] = 1.0;  // both players start at the center point

	PrintStream out = System.out;
	out.println("MODEL " + mo.modelName);
	out.println("Defender's efficiency phi="+ mo.phi);
	out.println("Attacker's allowed movement map:");
	out.println(mo.wMatrixToString(mo.w));
	out.println("Defender's allowed movement map:");
	out.println(mo.wMatrixToString(mo.w2));

	F2ArgPayoff test  = new F2ArgPayoff(mo, sym, jpv, true);

	int dim[] = { test.aScheme.nvar, test.bScheme.nvar};
	Constraint cons[] = {test.aScheme.constraint, test.bScheme.constraint};

	F2Arg.Res res = test.findSaddlePoint(dim, cons, F2Arg.LookFor.MIN, 0);
	System.out.println(test.statsReport());

	out.println("A=mouse, B=cat");
	out.println("min_a max_b at: " + res);

	String[] labels = {"A", "B"};
	for(int l=0; l<2; l++) {
	    MatrixData amat = new MatrixData(test.schemes()[l], res.ab[l]);
	    double [][][] z = amat.toDenseMatrices();
	    out.println(labels[l] + ".unseen=\n" + mo.matrixToString(z[0]));
	    out.println(labels[l] + ".seen  =\n" + mo.matrixToString(z[1]));
	}

    }

}
