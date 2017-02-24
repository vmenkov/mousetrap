package gridsearch;

import java.util.*;

/** Represents a constraint of the form sum_j (a[j]*x[j]) &le; b,
    with non-negative a[]
*/
class SingleConstraint extends Constraint {
    double a[];
    double b=1.0;
    /** Checks that the constraint is all-non-negative */
    public void validate() {
	for(int i=0; i<a.length; i++) {
	    if (a[i]<0) throw new IllegalArgumentException("Constraints with negative coefficients not supported");
	}
    }
    static SingleConstraint simplex(int n) {
	SingleConstraint c=new SingleConstraint();
	c.a = new double[n];
	for(int i=0; i<n; i++) c.a[i] = 1.0;
	c.validate();
	return c;
    }
    public boolean holds(double x[]) {
	double sum = 0;
	for(int i=0; i<a.length; i++) {
	    sum += a[i] * x[i];
	}
	return sum <= b;
    }

    public boolean holds(int x[]) {
	double sum = 0;
	for(int i=0; i<a.length; i++) {
	    sum += a[i] * x[i];
	}
	return sum <= b;
    }

    public Constraint constraintInt(Grid g) {
	if (a.length != g.dim()) throw new IllegalArgumentException();
	SingleConstraint q = new 	SingleConstraint();
	q.a = new double[a.length];
	q.b = 1.0;
	for(int i=0; i<a.length; i++) {
	    q.b -= g.corners[0].x(i);
	    q.a[i] = a[i] * g.cellWidth(i);
	}
	q.validate();
	return q;
    }

    public boolean equals(Object o) {
	if (!(o instanceof SingleConstraint)) return false;
	SingleConstraint c = (SingleConstraint)o;
	return Arrays.equals(c.a, a) && c.b==b;
    }

}

