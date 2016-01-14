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

	constrainedPlayerIsTheAttacker = true;
    
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


    /** Represents a pair of mixed strategies (the constrained
	player's p and the mobile platyer's q), as obtained by an
	optimization method.
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

    boolean constrainedPlayerIsTheAttacker = true;

    String nameOfConstrainedPlayer() {
	return  constrainedPlayerIsTheAttacker ? "attacker" : "defender";
    }

    String nameOfMobilePlayer() {
	return  constrainedPlayerIsTheAttacker ?  "defender" : "attacker";
    }


    /** Find Nash equilibrium mixed strategies for the situation
	when the mouse (invader) can play any of the holes listed in w[],
	and the payoff consists of the immediate payoff
	and the expected future benefits (in f[]).
	This method uses the peculiar structure of the
	payoff matrix ((I - phi * E) + e^T * w) to simplify
	computatoins.
	@param f future benefits 	
     */
    OptResults pOptimize(double [] f, int w[]) {
	if (!constrainedPlayerIsTheAttacker) throw new IllegalArgumentException("Special optimization is only supported when the constrained player is the attacker");
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

    /** Computes the constrained player's aggregate expected payoff of
	an n-round game, based on the two player's given mixed
	strategies (p and q) in the 1st round and the expected payoff
	for the (n-1)-round game (i.e., during the rest of the game.
	@param p The constrained player's play at the 1st round
	@param q The mobile player's play at the 1st round
	@param oldF the constrained player's aggregate expected payoff of an (n-1)-round game (the last n-1 rounds of the n-round game), already multiplied by r
     */
    double[] newF(double[][] p, double[][] q, double [] oldF) {
	double f[] = new double[h];
	for(int i = 0; i<h; i++) {


	    //		double caught = (i==j? phi : 0.0);
	    //	payoffMatrix[i][j] = 
	    //	    (constrainedPlayerIsTheAttacker? 1-caught : caught)+f[w[i]];

	    double s=0;
	    if (constrainedPlayerIsTheAttacker) {
		s = 1.0;
		for(int j=0; j<h; j++) s += p[i][j] * (- phi*q[i][j] + oldF[j]);
	    } else {
		for(int j=0; j<h; j++) s += p[i][j] * (phi*q[i][j] + oldF[j]);
	    }
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

    /** Computes the infinity-norm difference of two vectors */
    static double infNormDiff(double a[][], double b[][]) {
	double d = 0;
	for(int i=0; i<a.length; i++) {
	    double x = infNormDiff(a[i],b[i]);
	    if (x>d) d=x;
	}
	return d;
    }

    /** Converts sparse representation of a vector to dense.
	@param x[] an array of L elements
	@param h size of the dense-representation array to produce.
	@param w[] an array of L indexes, indicating where the values from x[] should go inthe "dense" array of h elements
	@return  An array that will contain h elements (among which L non-zeros)

     */
    private static double[] spreadArray(double[] x, int h, int[] w) {
	double[] y = new double[h];
	if (x.length != w.length) throw new IllegalArgumentException();
	for(int i=0; i<w.length; i++) {
	    y[ w[i]] = x[i];
	}
	return y;
    }

    /** Find Nash equilibrium mixed strategies for the situation
	when the constrained player can play any of the holes listed in w[],
	and the payoff consists of the immediate payoff
	and the expected future benefits (in f[]).
	This method works creating a generic payoff matrix,
	and optimizing using the Simplex algorithm.

	<p>
	There is a difference how this method operates when the constrained 
	player is the attacker vs. when it is the defender. In both cases
	the constrained player can only play some holes (the geometrically allowed
	once); but for the other (mobile) player the choice differs: in the former
	case the other player is the defender and it makes sense for it to
	only play the holes available to the attacker; in the latter case,
	the mobile player (the defender) will want to also play the holes *not* available
	to the constrained one (the attacker). (In fact, it only will play them,
	if such holes exist!).


	@param f future benefits 	
     */
    OptResults pOptimize2(double [] f, int w[]) {    
	final  int L = w.length;
	double[][] payoffMatrix = new double[L][];
	for(int i = 0; i<L; i++) {

	    if (constrainedPlayerIsTheAttacker) {
		payoffMatrix[i] =  new double[L];
		Arrays.fill( payoffMatrix[i], f[w[i]]);
		for(int j = 0; j<L; j++) {
		    double caught = (i==j? phi : 0.0);
		    payoffMatrix[i][j] += 1-caught;
		}
	    } else {
		payoffMatrix[i] =  new double[h];
		Arrays.fill( payoffMatrix[i], f[w[i]]);
		for(int j = 0; j<L; j++) {
		    double caught = (i==j? phi : 0.0);
		    payoffMatrix[i][w[j]] += caught;
		}
	    }
	}
	double[][] trans = transpose( payoffMatrix);
	mult(trans, -1);
	SimplexResults mouseRes = new SimplexResults(payoffMatrix);
	SimplexResults catRes = new SimplexResults(trans);
	//	System.out.println("M/C: " + mouseRes.maxval + " : " + (-catRes.maxval));
	double [] p = spreadArray(mouseRes.p, h, w);
	double [] q = constrainedPlayerIsTheAttacker ? spreadArray(catRes.p, h, w) :
	    catRes.p;
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

    void optimize(PrintStream out) {
	OptResults[] po = new OptResults[h];
	int n=0;
	double f[] = new double[h];
	double[][] p0=null, q0=null;
	double[] avgF = new double[h];

	out.println("MODEL " + modelName);
	out.println("Constrained player is the "+ nameOfConstrainedPlayer() +", mobile player is the " +
		    nameOfMobilePlayer() );

	out.println("Defender's efficiency phi="+ phi+", discount rate=" +r);
	out.println("Constrained player's allowed movement map:");
	out.println(wMatrixToString(w));

	String lab1 = "Constrained player (" + nameOfConstrainedPlayer() +")";
	String lab2 = "Mobile player ("+    nameOfMobilePlayer()+")";

	boolean conv = false;
	for(; n<100 && !conv; n++) {
	    out.println("---- " + (n+1) + "-round game: ------------------");
	    for(int i = 0; i<h; i++) {
		if (useSimplex) {
		    po[i] = pOptimize2(f, w[i]);
		} else {
		    po[i] = pOptimize(f, w[i]);
		}
	    }
	    double p[][] = OptResults.assembleP(po);
	    double q[][] = OptResults.assembleQ(po);
	    out.print(lab1 + " plays P=\n" + matrixToString(p));
	    out.print(lab2 + " plays Q=\n" + matrixToString(q));	    
	    double[] f1 = newF( p,  q, f);
	    for(int i = 0; i<h; i++) {
		avgF[i] = f1[i] / (n+1);
		f[i] = f1[i]*r;
	    }
	    out.println(lab1 + "'s avg payoff per round=");
	    for(int i = 0; i<h; i++) {
		out.print("\t" +format(avgF[i]));

	    }
	    out.println();

	    final double eps=1e-5;
	    if (p0!=null && q0!=null && infNormDiff(p0,p)<eps  && infNormDiff(q0,q)<eps) {
		out.println("Convergence on P and Q achieved within eps=" + eps);
		conv = true;
	    }
	    p0=p;
	    q0=q;
	}
	out.println("===== Approximating with rational numbers: ======");
	out.print("Approx P=\n" + matrixToString(approxRational(p0)));
	out.print("Approx Q=\n" + matrixToString(approxRational(q0)));
	Rational[] ravgF = approxRational(avgF);
	out.println(lab1 + "'s approx avg payoff per round=");
	for(int i = 0; i<h; i++) {
	    out.print("\t" +ravgF[i].toString());
	    
	}
	out.println();

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


    static public void main(String argv[]) throws FileNotFoundException {


	ParseConfig ht = new ParseConfig();
	boolean mobileCat = ht.getOption("mobileCat", true);	
	String fname = ht.getOption("out", "mousetrap.out");

	// let's always produce UNIX-style output, even when running under
	// MS Windows
	String sep = System.setProperty("line.separator", "\n");
	//System.out.println("Separator was ["+sep+"]");
	//sep = System.getProperty("line.separator");
	//System.out.println("Separator is ["+sep+"]");

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

	

	File f=new File(fname);
	System.out.println("Output will go to file " + f);
	PrintStream out = new PrintStream(new FileOutputStream(f));

	for(int i=0; i<mos.length; i++) {	    
	    out.println("============= System no. " + (i+1) + "=======================================");
	    Mousetrap mo = mos[i];
	    mo.constrainedPlayerIsTheAttacker = mobileCat;
	    mo.optimize(out);
	}
	out.close();

    }
   
}
