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

	r = 1;
	for(int i=0; i<alpha.dim(); i++) {
	    r /= 3.0;
	    double d = alpha.x(i) - r;
	    w -= d*d;
	}
	
	return w;
    }

    
    static public void test1(String[] argv) {
	F2Arg f = new F2ArgTest1();
	int n=5;
	//	boolean simplex = false;
	Res res;
	boolean ss[] = {false, true};
	for(boolean simplex: ss) {
	    String lab = simplex? "simplex" : "cube";
	    res = f.optimizeOverOneVar(simplex, ParVec.zero(n), LookFor.MIN, 1);
	    System.out.println("Minimized on "+lab+" at: " + res);
	    res = f.optimizeOverOneVar(simplex, ParVec.zero(n), LookFor.MAX, 1);
	    System.out.println("Maximized on "+lab+" at: " + res);
	}


    }

    static public void test2(String[] argv) {
	F2Arg f = new F2ArgTest1();
	int n=5;
	
	boolean ss[] = {false, true};
	for(boolean simplex: ss) {
	    Res res = f.findSaddlePoint(simplex, n, LookFor.MAX, 0);
	    System.out.println("max_a min_b at: " + res);
	    res = f.findSaddlePoint(simplex, n, LookFor.MIN, 1);
	    System.out.println("min_b max_a at: " + res);
	}
    }

    static public void main(String[] argv) {
	test1(argv);
    }


    
}
	    

