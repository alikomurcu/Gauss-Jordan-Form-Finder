double* Gauss(int n)
{
    double** a = new double* [n];  // The size of the matrix n*n
    for(int i=0; i<n; i++) a[i] = new double[n];
    double b[n];  // The result matrix in the form X*a=b
    for(int k=0; k<n; k++)
    {
        // Pivot is the a[k][k] we need to make it 1.
        for(int j=0; j<n; j++)
        {
            a[k][j] = a[k][j]/a[k][k];  // Divide each row by pivot, so that the diagonal entries will be all 0s.
        }
        b[k] = b[k] / a[k][k];  // Divide b matrix also by pivot.

        // Make all entries 0 in corresponding column except the pivot.
        /* Making 0 algorithm is as follows:
            1) Iterate each row.
            2) Select the entry in corresponding column and row.
            3) Multiply that entry by the row of pivot (which is made 1 in above) and negate it. 
            (e.g. Say k=2, and we know a[2][2] = 1 by above. The pivot row is R(2).
            Multiply the row R(1) by -R(2) and add it to R(1)
            Skip multiplying R(2) by -R(2) .................  (since we don't want to make pivot 0)
            Multiply the row R(3) by -R(2) and add it to R(3)
            Multiply the row R(4) by -R(2) and add it to R(4)
            ...
            ...
            ...
            Multiply the row R(n-1) by -R(2) and add it to R(n-1))
        */
        for(int i=0; i<n; i++)
        {
            if(k==i) continue;

            for(int j=0; j<n; j++)
            {
                a[i][j] -= a[i][k]*a[k][j];  // This line is for -R(i)*R(k) -> R(i) and make the upper and lower entries 0 in pivot column.
            }

            b[i] -= a[i][k]*b[k];  // Do the same for the result matrix.
        }
        // In the end b matrix contains the solution matrix.
    }
    return b;
}