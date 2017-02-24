package gridsearch;

//import java.io.*;
import java.util.*;
//import mousetrap.*;

    /** A JointProbVector describes the probability distribution over
	all possible states of the two-player system.  xUnseen[i][j]
	is the probability of Player A being at i and Player B at j,
	and them not seeing each other; xSeen[i] is the probability of
	Players A and B both being at i, and seeing each other. All
	probabilities should sum to 1.0.
     */
public class JointProbVector {
    double [][] xUnseen;
    double [] xSeen;
    int n() { return xSeen.length; }
    
    /** Creates a vector of all zeros */
    JointProbVector(int n) {
	xUnseen = zeroMat(n,n);
	xSeen = new double[n];
    }

    /** Checks if the values sum to 1.0 (within a computational error) 
     */
    void validate() {
	double ss=sumSeen(), su=sumUnseen();
	double s=ss+su;	
	if (Math.abs(s - 1.0) > 1e-6) throw new IllegalArgumentException("Probabilities don't sum to 1.0. ss=" + ss+", su=" + su+", s=" + s);
    }

    double sumUnseen() {
	double su=0;
	for(int i=0; i< xUnseen.length; i++) {
	    for(int j=0; j< xUnseen[i].length; j++) {
		su += xUnseen[i][j];
	    }
	}
	return su;
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


	   x'_{ij} = (\sum_{kl} aUnseen_{ik} bUnseen_{jl} x_{kl} +
                  \sum_k aSeen_{ik} bSeen_{jk} xi_k ) - \delta_{ij} xi'_i.

	  Remember that sparse matrices A and B are stored by column.
    */
    JointProbVector apply( ParametrizedMatrix.MatrixData a, ParametrizedMatrix.MatrixData b, double phi) {
	JointProbVector res = new JointProbVector(n()); 

	validate();

	for(int k=0; k<a.w.length; k++) {
	    for(int pi=0; pi<a.w[k].length; pi++) {
		int i = a.w[k][pi];
		
		for(int pj=0; pj<b.w[k].length; pj++) {
		    int j = b.w[k][pj];
		    
		    double r = a.aSeen[k][pi] * b.aSeen[k][pj] * xSeen[k];
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
		    u[j][k] += b.aUnseen[l][pj] * xUnseen[k][l];
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
	    res.xUnseen[i][i] -= res.xSeen[i];
	    if (res.xUnseen[i][i]<0) throw new IllegalArgumentException("seen["+i+"]=" + res.xSeen[i] + ", unseen=" + res.xUnseen[i][i]);
	}
	
	res.validate();
	
	return res;
    }

	public String toString() {
	    StringBuffer b = new StringBuffer();
	    b.append("Unseen\n");
	    for(double[] q: xUnseen) {
		b.append("[" + Arrays.toString(q) + "]\n");
	    }
	    b.append("Seen\n");
	    b.append("[" + Arrays.toString(xSeen) + "]\n");
	    return b.toString();
	}


    }
    
