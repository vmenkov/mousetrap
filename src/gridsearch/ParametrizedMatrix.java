package gridsearch;

import java.io.*;
import java.util.*;

import mousetrap.*;
import util.ParseConfig;


/** Tools to convert a parameter vector (from the space over which we optimize)
    to an actual matrix.

    A ParametrizedMatrix object contains information about the
    transition matrix structure and the way the matrix elements can be
    reconstructed from a specified set of parameters. The matrix
    structure is represented by array w[][], with exactly the same
    semantics as in the Mousetrap class: w[j][] lists graph nodes that
    a player can reach in one step from node j. The parametrization
    mapping is in aposUnseen[][] and aposSeen[][]. Elements in
    aposUnseen[j][] contain the numeric IDs of parameters that control
    the matrix lements (transition probabilities) for transitions
    corresponding to w[j][], when the players have not seen each other
    at the previous step; elements in aposSeen[j][] are the parameter
    indexes for the matrix elements describing transition
    probabilities when the players have seen each other.
   
 */
public class ParametrizedMatrix  {

    static String sepline = "------------------------------------------------------------------------";

    /** An auxiliary class, used to package the resuts of one of the steps
	of creating the parametrization map. 
    */
    private static class AsgMap {
	int [][] aptr;
	int apos;

	/** Figures out how many new parameters are needed to describe
	a transition matrix with the sparsity pattern described by
	w[][]. Fills aptr[][], and adds constraints to mc. 
	@param sym The symmetry rules for the matrix. May be null (for no rules).
	@param apos0 How many parameters have already been used. The numbering
	of new parameters will start with this number.
	@param mc Will add new constraints (describing rules for new parameters)
	to this constraint set.
	*/
	AsgMap(int w[][], Symmetry sym, int apos0, MultiConstraint mc) {
	    apos = apos0;
	    aptr =   arrayStructureCopy(w,  Symmetry.NONE);
	    for(int k=0; k<w.length; k++) {
		int con[] = new int[w[k].length-1], conp=0;
		for(int i=0; i<w[k].length; i++) {
		    int a1 = (sym==null)? Symmetry.NONE: 
			sym.lookup( w,  aptr, k, i);
		    if (w[k][i] == k) {
			// the diagonal element = 1 - sum(others)
			aptr[k][i] = Symmetry.REST;
		    } else if (a1 != Symmetry.NONE)  {	    // use symmetry...
			aptr[k][i] = a1;
			con[conp++] = aptr[k][i];
		    } else {
			aptr[k][i] = apos ++;
			con[conp++] = aptr[k][i];
		    }
		}		
		mc.addSimplexConstraintIfUnique(con);

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

	System.out.println(sepline);
	System.out.println("Building map for unseen, start apos=0");
	AsgMap map1 = new AsgMap(w, sym, 0, mc);
	System.out.println("Building map for seen, start apos=" + map1.apos);
	AsgMap map2 = new AsgMap(w, sym, map1.apos, mc);
	aposUnseen = map1.aptr;
	aposSeen   = map2.aptr;
	nvar = map2.apos;
	constraint = mc;
    }

    /** An object of this class describes an actual pair of transition
	matrix for one player */
    static class MatrixData {
	final int[][] w;
	/** Sparse matrices whose structure (by column) is described by w[][] */
	double [][]  aUnseen, aSeen;

	/** Fills one transition matrix ("unseen" or "seen").
	    @param q parameter values
	    @param apos the map that explains how matrix elements are
	    computed from the parameters
	 */
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
		    throw new IllegalArgumentException("Negative diagonal value ("+a[k][diagPos]+") computed for k="+k );
		}
	    }
	    return a;
	}


	MatrixData(ParametrizedMatrix mi, ParVec qv) {
	    this(mi, qv.getX());
	}

	/** Fills the matrix based on the assignment map in mi, and
	    parameter values in q[] */
	MatrixData(ParametrizedMatrix mi, double [] q) {
	    w = mi.w;
	    if (q.length != mi.nvar) throw new IllegalArgumentException("var cnt mismatch");
	    //	    try {
		aUnseen =  fillData(mi.aposUnseen, q);
		aSeen   =  fillData(mi.aposSeen, q);
		/*} catch (IllegalArgumentException ex) {
		  System.out.println(sepline+"\n"+
				   "Error context: mi=\n" + mi +
				   "\n,q=" + Arrays.toString(q));
		throw ex;
		}*/
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
    static String report(int w[][], int aptr[][]) {
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

   
    public String toString() {
	return 
	    "(" + nvar + " parameters)\n" +
	    "Unseen:\n" + report(w, aposUnseen) + 
	    "Seen:\n" + report(w, aposSeen) +
	    "Constraints:\n" + constraint	    ;
    }

    public static class F2ArgPayoff extends F2Arg {
  	ParametrizedMatrix aScheme, bScheme;
	JointProbVector jpv0;
	final double phi;

  	ParametrizedMatrix[] schemes() {
	    return new  ParametrizedMatrix[] { aScheme, bScheme};
	}

	final boolean ev;
	final int maxT;

	/**
	   @param ev If true, measure long-term payoff, rather than 1-step immediate payoff
	 */
	F2ArgPayoff(Mousetrap2 mo, Symmetry sym, 
		    JointProbVector _jpv0, boolean _ev, int _maxT) {
	    ev = _ev;
	    maxT = _maxT;
	    phi = mo.phi;
	    jpv0 = _jpv0;
	    jpv0.validate();
	    aScheme = new ParametrizedMatrix(mo.w, sym);
	    bScheme = new ParametrizedMatrix(mo.w2, sym);
	    

	    System.out.println(sepline);
	    System.out.println("Player A parametrization map:");
	    System.out.println(aScheme);

	    System.out.println("Player B parametrization map:");
	    System.out.println(bScheme);
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

	long callsTotal=0, callsNoConv=0;
	long tTotal=0;

	String statsReport() {
	    return ev?
		"EV: " + callsTotal + " calls so far; failed to achieve convergence in " +  callsNoConv + " calls; <T> = "+((double)tTotal)/callsTotal:
		"Imediate payoff: " + callsTotal + " calls so far";
	}

	double f_longTerm(ParVec alpha, ParVec beta) {
	    //if (callsTotal % 1000000 ==0) System.out.println(statsReport());

	    callsTotal++;
	    MatrixData amat = new MatrixData(aScheme, alpha);
	    MatrixData bmat = new MatrixData(bScheme, beta);

	    final int avgT = 10; // averaging interval
	    int t=0;
	    double r=0, r0=0;
	    JointProbVector jpv = jpv0;
	    while(true) {
		t++;
		jpv = jpv.apply(amat, bmat, phi);		
		double sumSeen = jpv.sumSeen(), sumUnseen = jpv.sumUnseen();
		double ratio = sumSeen / (sumSeen + sumUnseen);
		if (Double.isInfinite(ratio)) {
		    System.out.println("Infinity encountered: t=" + t+", jpv.sumSeen() = " + jpv.sumSeen() + ", jpv.sumUnseen()=" + jpv.sumUnseen());
		    System.out.println("jpv=" + jpv);
		    System.out.println("amat=" + amat);
		    System.out.println("bmat=" + bmat);
		    
		    throw new IllegalArgumentException();
		}

		final double eps=1e-4;

		if (avgT==1) {
		    r = ratio;
		    if (t>1) {		    		    
			if (Math.abs(r-r0) < eps) break;
		    }
		    r0 = r;
		} else {
		    r += ratio;
		    if (t % avgT == 0) {
			if (t / avgT > 1) {
			    if (Math.abs(r-r0) < eps*avgT) break;
			}
			r0 = r;
			r = 0;
		    }
		}
		
		if (t/avgT>1 && t>maxT) {
		    callsNoConv ++;
		    break;
		}
	    }
	    tTotal += t;
	    return r / avgT;

	}
    }

    private static void reportResults(PrintStream out, Mousetrap2 mo, F2ArgPayoff test, F2Arg.Res res) {
	String[] labels = {"A", "B"};
	for(int l=0; l<2; l++) {
	    MatrixData amat = new MatrixData(test.schemes()[l], res.ab[l]);
	    double [][][] z = amat.toDenseMatrices();
	    out.println(labels[l] + ".unseen=\n" + mo.matrixToString(z[0]));
	    out.println(labels[l] + ".seen  =\n" + mo.matrixToString(z[1]));
	}
    }


    /** The function being optimized is the payoff for the defender (i.e. the
	second player). Thus we go for min_A max_B  and max_B min_A
    */
    static public void main(String [] argv) {
	PrintStream out = System.out;
	ParseConfig ht = new ParseConfig();
	F2Arg.initParams(ht);
	out.println(sepline);
	out.println(F2Arg.params);

	out.println(sepline);
	int maxT = ht.getOption("f.T", 2000);
	boolean ev = ht.getOption("f.ev", true);


	int h =ht.getOption("graph.h", 3);
	boolean cyclic =ht.getOption("graph.cyclic", false);
	if (h<3) throw new IllegalArgumentException("h=" + h + "; need h>=3");

	out.println(sepline);
	out.println("Creating " + (cyclic? "cyclic" : "linear") + " graph with " + h + " nodes");
	Mousetrap2 mo = (h==3 && !cyclic) ? Mousetrap2.mo1() :
	    Mousetrap2.moChain(4, cyclic, 1, 1);

	boolean dosym = ht.getOption("graph.sym", true);

	Symmetry sym = dosym? Symmetry.mirror(mo.h) : null;

	out.println(sepline);
	out.println("Optimizing for " +
		    (ev?"long-term payoff with T=" + maxT: "immediate payoff"));

	
	JointProbVector jpv = new JointProbVector(mo.h);
	final boolean uniform =  ht.getOption("f.uniform", true);;
	if (uniform) {
	    // both players start at the same point... equal prob for any such point
	    jpv.setUniformDiagUnseen();
	    out.println("The two players start at the same random position");
	} else {
	    int startPoint = mo.h/2;
	    jpv.xUnseen[startPoint][startPoint] = 1.0;  // both players start at the center point
	    out.println("The two players start at position " + startPoint);
	}

	out.println(sepline);
	out.println("MODEL " + mo.modelName);
	out.println("Defender's efficiency phi="+ mo.phi);
	out.println("Attacker's allowed movement map:");
	out.println(mo.wMatrixToString(mo.w));
	out.println("Defender's allowed movement map:");
	out.println(mo.wMatrixToString(mo.w2));

	out.println(sepline);
	out.println("Solutions' required symmetry: " + sym);


	F2ArgPayoff test  = new F2ArgPayoff(mo, sym, jpv, ev, maxT);

	int dim[] = { test.aScheme.nvar, test.bScheme.nvar};
	Constraint cons[] = {test.aScheme.constraint, test.bScheme.constraint};

	//------------------------- 
	F2Arg.Res res1 = test.findSaddlePoint(dim, cons, F2Arg.LookFor.MIN, 0);
	System.out.println(test.statsReport());

	out.println(sepline);
	out.println("A=mouse, B=cat");
	out.println("min_a max_b f(a,b) = f(a1, b1), at\n" + res1);
	reportResults(out,  mo, test,  res1);

/*
	F2Arg.Res resRev = 
	    test.optimizeOverOneVar(res1.ab[1],cons[0],dim[0],F2Arg.LookFor.MIN, 0, null);
	out.println("Reverse: if we fix the above B strategy (b1), and look for min_a f(a, b1), the results are:\n" + resRev);
*/
	
	//------------------------- 
	F2Arg.Res res2 = test.findSaddlePoint(dim, cons, F2Arg.LookFor.MAX, 1);
	System.out.println(test.statsReport());

	out.println(sepline);
	out.println("A=mouse, B=cat");
	out.println("max_b min_a  = f(a2, b2), at:\n" + res2);
	reportResults(out,  mo, test,  res2);

	/*
	resRev = 
	    test.optimizeOverOneVar(res2.ab[0],cons[1],dim[1],F2Arg.LookFor.MAX, 1, null);
	out.println("Reverse: if we fix the above A strategy (a2), and look for max_b f(a2, b), the results are:\n" + resRev);
	*/

	//------------------------- a1,b2
	
	double fCombo = test.f(res1.ab[0], res2.ab[1]);
	out.println("Combination result: f(a1,b2) = " + fCombo);

    }

}
