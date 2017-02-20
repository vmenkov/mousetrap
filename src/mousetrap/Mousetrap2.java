package mousetrap;

import java.util.*;
import java.text.*;
import java.io.*;

/** Cat and Mouse simulation for Paul Kantor. An instance of this 
    class describes a particular geometry. The two players are 
    described as "constrained" and "mobile".

    In this class both players are viewed as constrained. The first player
    is the attacker, the second, the defender.
*/
public class Mousetrap2 extends Mousetrap {
    
    /** The second part of the model's geometry: for each hole X, w2[X] contains the list
	of holes that can be played by the 2nd player at the
	next step after X */
    public final int[][] w2;

    /** Cat's (defender's) efficiency */
    //   final double phi = 1.0;
    /** Discount factor for adding expected future benefits to the current
	round's immediate payoff */
    //    final double r = 1.0;

    Mousetrap2(int [][]_w, int[][] _w2) {
	this("unnamed", null, _w, _w2);
    }

    /**
       @param _names List of hole names. If null is given, use default names.
       @param _w the model's geometry
     */
    Mousetrap2(String _modelName, String[] _names, int [][]_w, int[][] _w2) {
	super(_modelName, _names, _w);
	w2 = _w2;

	if (w2==null) throw new IllegalArgumentException("Invalid w2=null");
	if (w2.length !=h) throw new IllegalArgumentException("Invalid w2 (len=" + w2.length+")");

	int i=0;
	for(int[] z: w2) {
	    if (z==null) throw new IllegalArgumentException("w2: null in row " + i);
	    for(int k=0; k<z.length; k++) {
		if (z[k] < 0 || z[k]>=h || k>0 && z[k]<=z[k-1]) throw new IllegalArgumentException("Invalid data in w2[" + i + "]: values are not in order");
	    }
	    i++;
	}
    }

    /** Computes the first player's (attacker's) aggregate expected payoff of
	an n-round game, based on the two player's given mixed
	strategies (p and q) in the 1st round and the expected payoff
	for the (n-1)-round game (i.e., during the rest of the game.

	<center>
	f_{ij}^{(n)} = \Sum_{kl} {  p_{ij,k} q_{ij,l} (1 - phi*delta_{kl} + f_{kl}^{(n-1)} =
         1 - phi * \Sum_k  {  p_{ij,k} q_{ij,k} } + 
	 \Sum_{kl} {  p_{ij,k} q_{ij,l}  f_{kl}^{(n-1)}                  
	</center>

	@param p The constrained player's play at the 1st round
	@param q The mobile player's play at the 1st round
	@param oldF the constrained player's aggregate expected payoff of an (n-1)-round game (the last n-1 rounds of the n-round game), already multiplied by r
     */
    double[][] newF(double[][][] p, double[][][] q, double [][] oldF) {
	double f[][] = alloc2(h,h);
	for(int i = 0; i<h; i++) {
	    for(int j = 0; j<h; j++) {
		double caught = 0;
		for(int k=0; k<h; k++) caught +=  p[i][j][k] * q[i][j][k];
		double s = 1.0 - phi * caught;

		for(int k=0; k<h; k++) {
		    for(int l=0; l<h; l++) {
			s += p[i][j][k] * q[i][j][l] * oldF[k][l];
		    }
		}
		f[i][j] = s;
	    }
	}
	return f;
    }

    static class OptResults2 {
	//	boolean qPat[];
	static double[][][] assembleP(OptResults po[][]) {
	    double p[][][] = new double[po.length][][];
	    for(int i=0; i<po.length; i++) {
		p[i] = new double[po[i].length][];
		for(int j=0; j<po[i].length; j++) {	
		    p[i][j] = po[i][j].p;
		}
	    }
	    return p;
	}
  	static double[][][] assembleQ(OptResults po[][]) {
	    double q[][][] = new double[po.length][][];
	    for(int i=0; i<po.length; i++) {
		q[i] = new double[po[i].length][];
		for(int j=0; j<po[i].length; j++) {	
		    q[i][j] = po[i][j].q;
		}
	    }
	    return q;
	}
    }



    /** Find Nash equilibrium mixed strategies for two constrained
	players.  The 1st player (the attacker) can play any of the
	holes listed in w[], the 2nd player (the defender) can play
	any of the holes listed in w2[].  The payoff consists of the
	immediate payoff and the expected future benefits (in f[]).
	This method works creating a generic payoff matrix, and
	optimizing using the Simplex algorithm.

	@param f f[j][k] contains future benefits for the first player (attacker) if the two players play (j,k)	
     */
    OptResults pOptimize2(double [][] f, int w[], int w2[]) {    
	final  int L1 = w.length, L2 = w2.length;
	double[][] payoffMatrix = new double[L1][];
	for(int i = 0; i<L1; i++) {
	    payoffMatrix[i] = arrayExtract( f[w[i]], w2);
	    for(int j = 0; j<L2; j++) {
		double caught = (w[i]==w2[j] ? phi : 0.0);
		payoffMatrix[i][j] += 1-caught;
	    }
	}
	double[][] trans = transpose( payoffMatrix);
	mult(trans, -1);
	SimplexResults mouseRes = new SimplexResults(payoffMatrix);
	SimplexResults catRes = new SimplexResults(trans);
	//	System.out.println("M/C: " + mouseRes.maxval + " : " + (-catRes.maxval));
	double [] p = spreadArray(mouseRes.p, h, w);
	double [] q = spreadArray(catRes.p, h, w2);
	OptResults res = new OptResults(p,q);
	return res;
 	
    }



    String matrixToString2(double [][][]p) {
	StringBuffer s = new	StringBuffer();
	for(int i = 0; i<h; i++) {
	    for(int j = 0; j<h; j++) {
		s.append(names[i] + "\t" + names[j] + ":");
		for(int k=0; k<h; k++) s.append("\t" + format(p[i][j][k]));
		s.append("\n");
	    }
	}
	return s.toString();
    }

    String matrixToString2(Rational [][][]p) {
	StringBuffer s = new	StringBuffer();
	for(int i = 0; i<h; i++) {
	    for(int j = 0; j<h; j++) {
		s.append(names[i] + "\t" + names[j] + ":");
		for(int k=0; k<h; k++) s.append("\t" + p[i][j][k]);
		s.append("\n");
	    }
	}
	return s.toString();
    }

    double[][] alloc2(int n1, int n2) {
	double[][] f = new double[n1][];
	for(int i = 0; i<n1; i++) f[i] = new double[n2];
	return f;
    }

    void optimize(PrintStream out, double eps) {
	OptResults[][] po = new OptResults[h][];
	int n=0;
	double[][] f = alloc2(h,h);
	double[][][] p0=null, q0=null;
	double[][] avgF = alloc2(h,h);

	for(int i = 0; i<h; i++) {
	    po[i] = new OptResults[h];
	}

	out.println("MODEL " + modelName);
	out.println("Defender's efficiency phi="+ phi+", discount rate=" +r);
	out.println("Attacker's allowed movement map:");
	out.println(wMatrixToString(w));
	out.println("Defender's allowed movement map:");
	out.println(wMatrixToString(w2));

	String lab1 = "Attacker";
	String lab2 = "Defender";

	boolean conv = false;
	double avgAvgF=0;
	final int T=100;
	for(; n<T && !conv; n++) {
	    out.println("---- " + (n+1) + "-round game: ------------------");
	    for(int i = 0; i<h; i++) {
		for(int j = 0; j<h; j++) {
		    po[i][j] = pOptimize2(f, w[i], w2[j]);
		}
	    }
	    double p[][][] = OptResults2.assembleP(po);
	    double q[][][] = OptResults2.assembleQ(po);

	    double[][] f1 = newF( p, q, f);
	    avgAvgF=0;
	    for(int i = 0; i<h; i++) {
		for(int j = 0; j<h; j++) {
		    avgF[i][j] = f1[i][j] / (n+1);
		    f[i][j] = f1[i][j]*r;
		    avgAvgF += avgF[i][j];
		}
	    }
	    avgAvgF /= (h*h);


	    if (p0!=null && q0!=null && infNormDiff(p0,p)<eps  && infNormDiff(q0,q)<eps) {
		out.println("Convergence on P and Q achieved within eps=" + eps);
		conv = true;
	    }


	    if (n<5 || n<100 && (n+1) % 10 == 0 || (n+1)% 100 == 0 || conv) {
		out.print(lab1 + " plays P=\n" + matrixToString2(p));
		out.print(lab2 + " plays Q=\n" + matrixToString2(q));	    
	    }

	    out.println(lab1+"'s avg payoff per round (hole avg="+avgAvgF+")=");
	    for(int i = 0; i<h; i++) {
		out.print(names[i]);
		for(int j = 0; j<h; j++) {
		    out.print("\t" +format(avgF[i][j]));
		}
		out.println();
	    }
	    
	    p0=p;
	    q0=q;
	}
	if (!conv) out.println("NO CONVERGENCE ACHIEVED in " + T + " rounds");
	out.println("===== Approximating with rational numbers: ======");
	out.print("Approx P=\n" + matrixToString2(approxRational(p0)));
	out.print("Approx Q=\n" + matrixToString2(approxRational(q0)));
	Rational[][] ravgF = approxRational(avgF);
	out.println(lab1 + "'s approx avg payoff per round=");
	for(int i = 0; i<h; i++) {
	    out.print(names[i]);
	    for(int j = 0; j<h; j++) {
		out.print("\t" +ravgF[i][j].toString());
	    }
	    out.println();	    
	}
	out.println("Avg payoff per round, averaged for all starting holes = " + format(avgAvgF) + " ~= " + approxRational(avgAvgF));
    }

    /** Paul's 3-wall model */
    public static Mousetrap2 mo1() {
	int [][] w=
	    new int [][] {new int[] {0, 1},
			  new int[] {0, 1, 2}, 
			  new int[] {1, 2}};
	Mousetrap2 mo = new
	    Mousetrap2("Three holes",
		      new String[] {"L", "C", "R"},
		      w, w);
	return mo;
    }

 
    /** Creates an open or close chain of holes. Both players are constrained,
	but may move with different speeds.
     */
    static Mousetrap2 moChain(int h, boolean cyclic, int speed1, int speed2) {
	//int [] speeds = { speed1, speed2};
	//int ws[][][] = { makeChainW(h, cyclic, speed1), makeChainW(h, cyclic, speed2)};
	return  new  Mousetrap2((cyclic? "Circular chain" : "Open chain") +
				" of " + h + " holes", null, 
				makeChainW(h, cyclic, speed1), 
				makeChainW(h, cyclic, speed2));
    }


    /** Long chain of holes */
    //    static Mousetrap2 mo2(int h) {
    //	return  moChain( h, false, 1, 1);
    //    }
    

    /** One dead-end hole next to 9 all-connected holes */
    static Mousetrap2 mo3() {
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
	return  new  Mousetrap2("One dead-end node + a clique", null, w, w);
    }

    /* Four-pointed star */
    static Mousetrap2 mo4() {
	int [][] w=
		      new int [][] { {0, 1, 2, 3, 4},
				     {0, 1}, 
				     {0, 2},
				     {0, 3},
				     {0, 4}};
	Mousetrap2 mo = new
	    Mousetrap2("4-pointed star",
		      new String[] {"C", "N", "E", "S", "W"},
		      w,w);
	return mo;
    }

    /* Four-pointed star, longer arms */
    static Mousetrap mo5() {
	int w[][]=
	     { {0, 1, 3, 5, 7},
				     {0, 1, 2}, 
				     {1, 2},
				     {0, 3, 4},
				     {3, 4},
				     {0, 5, 6},
				     {5, 6},
				     {0, 7, 8},
				     {7, 8}};
	Mousetrap2 mo = new
	    Mousetrap2("5-pointed star, arm length=2",
		      new String[] {"C", "N1", "N2", "E1", "E2", "S1","S2", "W1", "W2"}, w,w);
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


    /**Usage:
       java -classpath lib/mousetrap.jar mousetrap.Mousetrap2
     */
    static public void main(String argv[]) throws FileNotFoundException {


	ParseConfig ht = new ParseConfig();
	boolean mobileCat = ht.getOption("mobileCat", true);	
	String fname = ht.getOption("out", "mousetrap.out");

	final double eps=ht.getOptionDouble("eps", 1e-5);

	// let's always produce UNIX-style output, even when running under
	// MS Windows
	String sep = System.setProperty("line.separator", "\n");
	//System.out.println("Separator was ["+sep+"]");
	//sep = System.getProperty("line.separator");
	//System.out.println("Separator is ["+sep+"]");

	Mousetrap2 mos[] = {mo1(),
			    moChain( 10, false, 1, 1),
			    moChain( 10, false, 1, 2),
			    moChain( 10, false, 2, 1),
			    moChain( 10, false, 2, 2),

			    moChain( 10, true, 1, 1),
			    moChain( 10, true, 1, 2),
			    moChain( 10, true, 2, 1),
			    moChain( 10, true, 2, 2),
			    mo3(),
			    mo4(),
			    //mo5(),
			    //mo6()};
	
	};
	

	File f=new File(fname);
	System.out.println("Output will go to file " + f);
	PrintStream out = new PrintStream(new FileOutputStream(f));

	for(int i=0; i<mos.length; i++) {	    
	    out.println("============= System no. " + (i+1) + "=======================================");
	    Mousetrap2 mo = mos[i];
	    mo.optimize(out, eps);
	}
	out.close();

    }
   
}
