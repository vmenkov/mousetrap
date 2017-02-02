package gridsearch;

public class F2ArgTest1 extends F2Arg {

    /** minimum at x=(0.5, 0.25, 0.125, ...) */
    double f(ParVec alpha, ParVec beta) {
	double w=0;
	double r=1;
	for(int i=0; i<beta.dim(); i++) {
	    r *= 0.5;
	    double d = beta.x(i) - r;
	    w += d*d;
	}
	return w;
    }

    
    static public void main(String[] argv) {
	F2Arg f = new F2ArgTest1();
	int n=5;
	Res res = f.optimizeOverOneVar(ParVec.zero(n), LookFor.MIN, 1);
	System.out.println("Minimized at: " + res);
	res = f.optimizeOverOneVar(ParVec.zero(n), LookFor.MAX, 1);
	System.out.println("Maximized at: " + res);
    }
}
	    

