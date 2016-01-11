package mousetrap;

import java.util.*;
import java.text.*;
import java.io.*;

/** Cat and Mouse simulation for Paul Kantor. An instance of this 
    class describes a particular geometry. The two players are 
    described as "constrained" and "mobile".
*/
public class Mousetrap {
    
    /** The human-readable model name */
    final String modelName;
    
    /** The number of holes */
    final int h;
    /** The names of holes (by default, simply "0", "1", etc */
    final String names[];
    /** The model's geometry: for each hole X, w[X] contains the list
	of holes that can be played by the "constrained" player at the
	next step after X */
    final int[][] w;

    /** Cat's (defender's) efficiency */
    final double phi = 1.0;
    /** Discount factor for adding expected future benefits to the current
	round's immediate payoff */
    final double r = 1.0;

    Mousetrap(int [][]_w) {
	this("unnamed", null, _w);
    }

    /**
       @param _names List of hole names. If null is given, use default names.
       @param _w the model's geometry
     */
    Mousetrap(String _modelName, String[] _names, int [][]_w) {
	modelName = _modelName;
	w = _w;
	if (w==null || w.length ==0) throw new IllegalArgumentException("Invalid w");
	h = w.length;

	if (_names==null)  {
	    names= new String[h];
	    for(int i=0; i<h; i++) names[i] = "" + i;
	} else if (_names.length != h) {
	    throw new IllegalArgumentException("names[], w[] length mismatch");
	} else {
	    names = _names;
	}

	int i=0;
	for(int[] z: w) {
	    if (z==null) throw new IllegalArgumentException("null in row " + i);
	    for(int k=0; k<z.length; k++) {
		if (z[k] < 0 || z[k]>=h || k>0 && z[k]<=z[k-1]) throw new IllegalArgumentException("Invalid data in w[" + i + "]: values are not in order");
	    }
	    i++;
	}
    }

    
   /** Descending order with respect to the values in the array. */
    static class DescendingComparator implements Comparator<Integer> {
        private final double[] f;
        DescendingComparator(double[] _f) { f = _f; }
        public int compare(Integer o1, Integer o2) {
            double x = f[o2.intValue()] - f[o1.intValue()];
            return (x>0) ? 1 : (x<0) ? -1 : 0;
        }

	private static Integer[] list(int n)  {
	    Integer[] v = new Integer[n];
	    for(int i=0; i<n; i++) v[i] = new Integer(i);
	    return v;
	}
	/** Returns a list of integers which indicate the positions
	    of largest, 2nd largest etc elements of f[w[*]], in descending
	    order */
 	int[] sortIndexes(int w[])  {
	    Integer[] v = new Integer[w.length];
	    for(int i=0; i<w.length; i++) v[i] = new Integer(w[i]);
	    Arrays.sort(v, this);
	    int q[] = new int[v.length];
	    for(int i=0; i<v.length; i++) q[i] = v[i].intValue();
	    return q;
	}
   }


    /** Represents a pair of mixed strategies (the Mouse's p and the Cat's 
	q), as obtained by an optimization method.
     */
    static class OptResults {
	double[] p, q;
	OptResults(	double[] _p,	double[] _q) {
	    p = _p;
	    q = _q;
	}
	//	boolean qPat[];
	static double[][] assembleP(OptResults po[]) {
	    double p[][] = new double[po.length][];
	    for(int i=0; i<po.length; i++) p[i] = po[i].p;
	    return p;
	}
  	static double[][] assembleQ(OptResults po[]) {
	    double q[][] = new double[po.length][];
	    for(int i=0; i<po.length; i++) q[i] = po[i].q;
	    return q;
	}
   }

    /** Find Nash equilibrium mixed strategies for the situation
	when the mouse can play any of the holes listed in w[],
	and the payoff consists of the immediate payoff
	and the expected future benefits (in f[]).
	This method uses the peculiar structure of the
	payoff matrix ((I - phi * E) + e^T * w) to simplify
	computatoins.
	@param f future benefits 	
     */
    OptResults pOptimize(double [] f, int w[]) {
	int[] fi = (new DescendingComparator(f)).sortIndexes(w);
	double sumF = 0;
	final int h1 = fi.length;
	boolean[] setHi = new boolean[h];
	int l=0, setCnt=0;
	while( l<h1 && sumF < phi + setCnt * f[ fi[l] ]) {
	    sumF += f[ fi[l] ];
	    setHi[fi[l]] = true;
	    setCnt ++;
	    l++;
	}

	double[] p = new double[h];
	final double hi = 1.0/setCnt;
	for(int i=0;i<h; i++) p[i] = setHi[i] ? hi : 0;

	// apportion q
	double fStar = (sumF - phi) / setCnt;
	double[] q = new double[h];
	for(int i=0;i<h; i++) q[i] = setHi[i] ? (f[i]-fStar)/phi : 0;
		

	OptResults res = new OptResults(p,q);
	return res;
    }

    /** Computes the Mouse's aggregate expected payoff of an n-round game. 
	@param p The Mouse's play at the 1st round
	@param q The Cat's play at the 1st round
	@param oldF the Mouse's aggregate expected payoff of an (n-1)-round game (the last n-1 rounds of the n-round game), already multiplied by r
     */
    double[] newF(double[][] p, double[][] q, double [] oldF) {
	double f[] = new double[h];
	for(int i = 0; i<h; i++) {
	    double s = 1.0;
	    for(int j=0; j<h; j++) s += p[i][j] * (- phi*q[i][j] + oldF[j]);
	    f[i] = s;
	}
	return f;
    }

    static DecimalFormat fmt = new DecimalFormat("0.0000");
    private String format(double x) {
	return (x==0.0) ? "0" : fmt.format(x);
    }

    String wMatrixToString(int [][]w) {
	StringBuffer s = new	StringBuffer();
	for(int i = 0; i<h; i++) {
	    s.append(names[i] + " :");
	    char[] z = new char[h];

	    for(int j=0; j<h; j++) z[j] = '-';
	    for(int j=0; j<w[i].length; j++) z[w[i][j]] = '*';
	    s.append(new String(z) + "\n");
	}
	return s.toString();
    }

    String matrixToString(double [][]p) {
	StringBuffer s = new	StringBuffer();
	for(int i = 0; i<h; i++) {
	    s.append(names[i] + " :");
	    for(int j=0; j<h; j++) s.append("\t" + format(p[i][j]));
	    s.append("\n");
	}
	return s.toString();
    }

    String matrixToString(Rational [][]p) {
	StringBuffer s = new	StringBuffer();
	for(int i = 0; i<h; i++) {
	    s.append(names[i] + " :");
	    for(int j=0; j<h; j++) s.append("\t" + p[i][j]);
	    s.append("\n");
	}
	return s.toString();
    }


    static double infNormDiff(double a[], double b[]) {
	double d = 0;
	for(int i=0; i<a.length; i++) {
	    double x = Math.abs(a[i]-b[i]);
	    if (x>d) d=x;
	}
	return d;
    }

    static double infNormDiff(double a[][], double b[][]) {
	double d = 0;
	for(int i=0; i<a.length; i++) {
	    double x = infNormDiff(a[i],b[i]);
	    if (x>d) d=x;
	}
	return d;
    }


    /** Find Nash equilibrium mixed strategies for the situation
	when the mouse can play any of the holes listed in w[],
	and the payoff consists of the immediate payoff
	and the expected future benefits (in f[]).
	This method works creating a generic payoff matrix,
	and optimizing using the Simplex algorithm.
	@param f future benefits 	
     */
    OptResults pOptimize2(double [] f, int w[]) {    
	final  int L = w.length;
	double[][] payoffMatrix = new double[L][];
	for(int i = 0; i<L; i++) {
	    payoffMatrix[i] =  new double[L];
	    for(int j = 0; j<L; j++) {
		payoffMatrix[i][j] = 1 - (i==j? phi : 0.0) + f[w[i]];
	    }
	}
	double[][] trans = transpose( payoffMatrix);
	mult(trans, -1);
	SimplexResults mouseRes = new SimplexResults(payoffMatrix);
	SimplexResults catRes = new SimplexResults(trans);
	//	System.out.println("M/C: " + mouseRes.maxval + " : " + (-catRes.maxval));
	double [] p = new double[h], q=new double[h];
	for(int i = 0; i<L; i++) {
	    p[ w[i]] = mouseRes.p[i];
	    q[ w[i]] = catRes.p[i];
	}
	OptResults res = new OptResults(p,q);
	return res;
 	
    }

    static double[][] transpose(double[][] a) {
	final int L=a.length; 
	final int H=a[0].length;
	double b[][] = new double[H][];
	for(int i=0; i<H; i++) {
	    b[i] = new double[L];
	    for(int k=0; k<L; k++) {
		b[i][k] = a[k][i];
	    }	   
	}	
	return b;
    }

    static void mult(double[][] a, double c) {
	for(int i=0; i<a.length; i++) {
	    for(int k=0; k<a[i].length; k++) {
		a[i][k] *= c;
	    }
	}
    }

    static final boolean useSimplex = true;

    void optimize() {
	OptResults[] po = new OptResults[h];
	int n=0;
	double f[] = new double[h];
	double[][] p0=null, q0=null;
	double[] avgF = new double[h];

	System.out.println("MODEL " + modelName);
	System.out.println("Defender's efficiency phi="+ phi+", discount rate=" +r);
	System.out.println("Constrained player's allowed movement map:");
	System.out.println(wMatrixToString(w));


	boolean conv = false;

	for(; n<100 && !conv; n++) {
	    System.out.println("---- " + (n+1) + "-round game: ------------------");
	    for(int i = 0; i<h; i++) {
		if (useSimplex) {
		    po[i] = pOptimize2(f, w[i]);
		} else {
		    po[i] = pOptimize(f, w[i]);
		}
	    }
	    double p[][] = OptResults.assembleP(po);
	    double q[][] = OptResults.assembleQ(po);
	    System.out.print("Mouse plays P=\n" + matrixToString(p));
	    System.out.print("Cat   plays Q=\n" + matrixToString(q));	    
	    double[] f1 = newF( p,  q, f);
	    for(int i = 0; i<h; i++) {
		avgF[i] = f1[i] / (n+1);
		f[i] = f1[i]*r;
	    }
	    System.out.println("Mouse's avg payoff per round=");
	    for(int i = 0; i<h; i++) {
		System.out.print("\t" +format(avgF[i]));

	    }
	    System.out.println();

	    final double eps=1e-5;
	    if (p0!=null && q0!=null && infNormDiff(p0,p)<eps  && infNormDiff(q0,q)<eps) {
		System.out.println("Convergence on P and Q achieved within eps=" + eps);
		conv = true;
	    }
	    p0=p;
	    q0=q;
	}
	System.out.println("===== Approximating with rational numbers: ======");
	System.out.print("Approx P=\n" + matrixToString(approxRational(p0)));
	System.out.print("Approx Q=\n" + matrixToString(approxRational(q0)));
	Rational[] ravgF = approxRational(avgF);
	System.out.println("Mouse's approx avg payoff per round=");
	for(int i = 0; i<h; i++) {
	    System.out.print("\t" +ravgF[i].toString());
	    
	}
	System.out.println();

    }

    /** Paul's 3-wall model */
    static Mousetrap mo1() {
	Mousetrap mo = new
	    Mousetrap("Three holes",
		      new String[] {"L", "C", "R"},
		      new int [][] {new int[] {0, 1},
				    new int[] {0, 1, 2}, 
				    new int[] {1, 2}});
	return mo;
    }

    /** Long chain of holes */
    static Mousetrap mo2() {
	int h = 10;
	int w[][] = new int [h][];
	for(int i=0; i<h; i++) {
	    w[i] = (i==0)? new int[] {i, i+1} :
	    (i==h-1)?  new int[] {i-1, i} :
	    new int[] {i-1, i, i+1};
	}
	return  new	    Mousetrap("Chain of " + h + " holes", null, w);
    }


    /** One dead-end hole next to 9 all-connected holes */
    static Mousetrap mo3() {
	int h = 10;
	int w[][] = new int [h][];
	for(int i=0; i<h; i++) {
	    if (i==0) {
		w[i] = new int[] {i, i+1};
	    } else if (i==1) {
		w[i] = new int[h];
		for(int k=0; k<h; k++) w[i][k] = k;
	    } else {
		w[i] = new int[h-1];
		for(int k=0; k<h-1; k++) w[i][k] = k+1;		
	    }
	}
	return  new  Mousetrap("One dead-end node + a clique", null, w);
    }

    /* Four-pointed star */
    static Mousetrap mo4() {
	Mousetrap mo = new
	    Mousetrap("4-pointed star",
		      new String[] {"C", "N", "E", "S", "W"},
		      new int [][] { {0, 1, 2, 3, 4},
				     {0, 1}, 
				     {0, 2},
				     {0, 3},
				     {0, 4}});
	return mo;
    }

    /* Four-pointed star, longer arms */
    static Mousetrap mo5() {
	Mousetrap mo = new
	    Mousetrap("5-pointed star, arm length=2",
		      new String[] {"C", "N1", "N2", "E1", "E2", "S1","S2", "W1", "W2"},
		      new int [][] { {0, 1, 3, 5, 7},
				     {0, 1, 2}, 
				     {1, 2},
				     {0, 3, 4},
				     {3, 4},
				     {0, 5, 6},
				     {5, 6},
				     {0, 7, 8},
				     {7, 8}});
	return mo;
    }

    /** Star with a "clique" heart. Each of n rays consists of r holes. */
    static Mousetrap mo6() {
       final int n = 5, r=3;
       String modelName = "" + n + "-pointed star, arm length=" + r + " with a clique heart";

       int h=n*r;
       String[] names = new String [h];
       int w[][] = new int [h][];
       int j=0;
       for(int k=0; k<r; k++) { // radial var
	   for(int i=0; i<n; i++) {   // angular var
	       char x = (char)('A' + i);
	       names[j] = ("" + x) + k;
	       if (k==0) {
		   w[j] = new int[n+1];
		   int l=0;
		   for(; l<n; l++) w[j][l] = l;
		   w[j][l] = j + n;
	       } else if (k<r-1) {
		   w[j] = new int[]{ j-n, j, j+n};
	       } else {
		   w[j] = new int[]{ j-n, j};
	       }
	       j++;
	   }
       }
       return  new   Mousetrap(modelName, names, w);
     }

    /** Approximately convert a floating-point value to a rational, for use
	in the simplex code. A kludge, of course.
     */       
    static Rational approxRational(double x) {
	if (x==0) return Rational.ZERO;
	final int n0 = 360 * 3;
	int bestn = 0;
	double minerr = 0;
	for(int n=1;n<= n0; n++) {
	    double q=n*x;
	    int iq = (int)Math.round(q);
	    double err = Math.abs( q - iq)/Math.abs(q);
	    if (err < 1e-6) return new Rational( iq, n);
	    if (bestn == 0 || err < minerr) {
		bestn = n;
		minerr = err;
	    }
	}
	int iq = (int)Math.round(bestn*x);	
	return new Rational( iq, bestn);
    }

    static Rational[] approxRational(double[] x) {
	Rational[] r = new Rational[x.length];
	for(int i=0; i<x.length; i++) r[i] = approxRational(x[i]);
	return r;
    }
    static Rational[][] approxRational(double[][] x) {
	Rational[][] r = new Rational[x.length][];
	for(int i=0; i<x.length; i++) r[i] = approxRational(x[i]);
	return r;
    }


    static public void main(String argv[]) {

	Mousetrap mos[] = {mo1(),
		mo2(),
		mo3(),
		mo4(),
		mo5(),
		mo6()};

	//Mousetrap mo = mo2();
	//Mousetrap mo = mo3();
	//Mousetrap mo = mo4();
	//Mousetrap mo = mo5();
	//	Mousetrap mo = mo6();

	
	for(int i=0; i<mos.length; i++) {
	    System.out.println("============= System no. " + (i+1) + "=======================================");
	    Mousetrap mo = mos[i];
	    mo.optimize();
	}

    }
   
}
