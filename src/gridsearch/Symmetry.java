package gridsearch;

//import java.io.*;
//import java.util.*;

/** Describes some form of symmetry of the graph (a homomorphism).     
    When optimizing for the players' strategies, we only look for
    strategies that are similarly symmetric (i.e. invariant with
    respect to the graph's mapping onto itself by the homomorphism).
*/
 class Symmetry {
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

