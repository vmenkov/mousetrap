package mousetrap;

import java.util.*;
import java.text.*;
import java.io.*;

/** Cat and Mouse simulation for Paul Kantor. An instance of this 
    class describes a particular geometry. */
public class Mousetrap {
    
    /** The number of holes */
    final int h;
    /** The names of holes (by default, simply "0", "1", etc */
    final String names[];
    /** The model's geometry: for each hole X, w[X] contains the list
	of holes that can be played at the next step after X */
    final int[][] w;

    /** Cat's efficiency */
    final double phi = 1.0;
    /** Discount factor for adding expected future benefits to the current
	round's immediate payoff */
    final double r = 1.0;

    Mousetrap(int [][]_w) {
	this(null, _w);
    }

    /**
       @param _names List of hole names. If null is given, use default names.
       @param _w the model's geometry
     */
    Mousetrap(String[] _names, int [][]_w) {
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
	q)
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

    /** Find Nash equilibrium mixed strategies for the situatio
	when the mouse can play any of the holes listed in w[],
	and the payoff consists of the immediate payoff
	and the expected future benefits (in f[])
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

    String matrixToString(double [][]p) {
	StringBuffer s = new	StringBuffer();
	for(int i = 0; i<h; i++) {
	    s.append(names[i] + " :");
	    for(int j=0; j<h; j++) s.append("\t" + format(p[i][j]));
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


    void optimize() {
	OptResults[] po = new OptResults[h];
	int n=0;
	double f[] = new double[h];
	double[][] p0=null, q0=null;

	System.out.println("Cat's efficiency phi="+ phi+", discount rate=" +r);

	for(; n<100; n++) {
	    System.out.println("---- " + (n+1) + "-round game: ------------------");
	    for(int i = 0; i<h; i++) {
		po[i] = pOptimize(f, w[i]);
	    }
	    double p[][] = OptResults.assembleP(po);
	    double q[][] = OptResults.assembleQ(po);
	    System.out.print("Mouse plays P=\n" + matrixToString(p));
	    System.out.print("Cat   plays Q=\n" + matrixToString(q));	    
	    double[] f1 = newF( p,  q, f);
	    double[] avgF = new double[h];
	    for(int i = 0; i<h; i++) {
		avgF[i] = f1[i] / (n+1);
		f[i] = f1[i]*r;
	    }
	    System.out.println("Mouse's avg payoff per round=");
	    for(int i = 0; i<h; i++) {
		System.out.print("\t" +format(avgF[i]));

	    }
	    System.out.println();

	    final double eps=0.0001;
	    if (p0!=null && q0!=null && infNormDiff(p0,p)<eps  && infNormDiff(q0,q)<eps) {
		System.out.println("Convergence on P and Q achieved within eps=" + eps);
		break;
	    }
	    p0=p;
	    q0=q;
	}
    }

    /** Paul's 3-wall model */
    static Mousetrap mo1() {
	Mousetrap mo = new
	    Mousetrap(new String[] {"L", "C", "R"},
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
	return  new	    Mousetrap(w);
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
	return  new	    Mousetrap(w);
    }

    /* Four-pointed star */
    static Mousetrap mo4() {
	Mousetrap mo = new
	    Mousetrap(new String[] {"C", "N", "E", "S", "W"},
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
	    Mousetrap(new String[] {"C", "N1", "N2", "E1", "E2", "S1","S2", "W1", "W2"},
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

    /** Star with a "clique" heart */
    static Mousetrap mo6() {
       final int n = 5, r=3;
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
       return  new   Mousetrap(names, w);
     }


    static public void main(String argv[]) {

	//Mousetrap mo = mo1();
	//Mousetrap mo = mo2();
	//Mousetrap mo = mo3();
	//Mousetrap mo = mo4();
	//Mousetrap mo = mo5();
	Mousetrap mo = mo6();

	mo.optimize();

    }
   
}