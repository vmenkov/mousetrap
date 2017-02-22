package gridsearch;

/** Represents a constraint of the form sum_j (a[j]*x[j]) &le; b,
    with non-negative a[]
*/
public interface Constraint {
    /** Checks that the constraint is all-non-negative */
    public void validate();
    public boolean holds(double x[]);
    public boolean holds(int x[]);
    /** This constraint is applied to the integer coordinates of grid points (on
	the [0..m] scale, rather than [0..1])  */
    public Constraint constraintInt(Grid g);
    public boolean equals(Object x);
}

