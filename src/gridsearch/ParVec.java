package gridsearch;

import java.util.Arrays;

public class ParVec {
    private double [] x;
    public int dim() { return x.length; }
    public double x(int i) { return x[i]; }
    public double[] getX() { return x;}
    public ParVec(int n) { x = new double[n]; }
    public ParVec(double c[]) {
	x = Arrays.copyOf(c, c.length);
    }
    public static ParVec zero(int n) { return new ParVec(n); }
    public static ParVec ones(int n) {
	ParVec a = new ParVec(n);
	for(int i=0; i<n; i++) a.x[i] = 1.0;
	return a;
    }

    public String toString() {
	StringBuffer b = new StringBuffer("[");
	for(double q: x) {
	    b.append(" " + q);
	}
	b.append("]");
	return b.toString();
    }

}
