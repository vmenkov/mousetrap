package gridsearch;

import java.io.*;
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
	    s.mapsto = new int[n];
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

	AsgMap map1 = new AsgMap(w, sym, 0, mc);
	AsgMap map2 = new AsgMap(w, sym, map1.apos, mc);
	aposUnseen = map1.aptr;
	aposSeen = map2.aptr;
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
		if (diagPos==Symmetry.NONE) throw new  IllegalArgumentException("No diagonal values found for k="+k);
		a[k][diagPos] = 1.0 - s;
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
	    return s;
	}

	double [][][] toDenseMatrices() {
	    return new double [][][] { toDenseMatrix(aUnseen),
				       toDenseMatrix(aSeen) };
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

    public static class F2ArgImmediatePayoff extends F2Arg {
  	ParametrizedMatrix aScheme, bScheme;
	JointProbVector jpv0;
	final double phi;

  	ParametrizedMatrix[] schemes() {
	    return new  ParametrizedMatrix[] { aScheme, bScheme};
	}
	F2ArgImmediatePayoff(Mousetrap2 mo, Symmetry sym, 
			     JointProbVector _jpv0) {
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
	    MatrixData amat = new MatrixData(aScheme, alpha);
	    MatrixData bmat = new MatrixData(bScheme, beta);
	    JointProbVector jpv = jpv0.apply(amat, bmat, phi);
	    return jpv.sumSeen();
	}
    }

    /** The function being optimized is the payoff for the defender (i.e. the
	second player). Thus we go for min_A max_B  and max_B min_A
    */
    static public void main(String [] argv) {
	Mousetrap2 mo = Mousetrap2.mo1();
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

	F2ArgImmediatePayoff test  = new F2ArgImmediatePayoff(mo, sym, jpv);



	int n = test.aScheme.nvar;
	F2Arg.Res res = test.findSaddlePoint(true, n, F2Arg.LookFor.MIN, 0);
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
