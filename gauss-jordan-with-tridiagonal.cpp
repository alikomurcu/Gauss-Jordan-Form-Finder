double* Gauss(double alpha, int n)
{
	double up[n-1], lo[n-1], di[n-1], *b, lo1, di1, up1, b1;  // the size of the matrix is (n-1)*(n-1)
    b = new double[n-1];
	for(int i=0; i<n-1; i++)  // Initializing the tri diagonal matrix by 3 seperate array.
	{
		up[i] = -1+alpha;  // Upper diagonal (entries a(1,2) a(2,3) a(3,4) ... a(n-1, n))
		di[i] = 1;  // Diagonal (entries a(1,1) a(2,2) a(3,3) ... a(n,n))
		lo[i] = -alpha;  // Lower diagonal (entries a(2,1) a(3,2) ... a(n, n-1))
		b[i] = 0;  // The result matrix which is in form X*a = b
	}

	b[0] = alpha;
	lo[0] = 0;
	lo[n-1] = 0;
	up[n-2] = 0;
	for(int i=0; i<n-2; i++) // n-2 because we will do last step as another part
	{
    // Divide the first row by the entry in corresponding diagonal.
		di[i] = 1;
		up[i] = up[i]/di[i];
		lo[i] = lo[i]/di[i];
		b[i] = b[i]/di[i];

    // Multiply a row by the negative value of the entry in lower diagonal and add it to corresponding entry in diagonal. In other words, make lower diagonal 0.
		lo1 = -di[i]*lo[i+1] + lo[i+1];
		if(i != n-3) di1 = -up[i]*lo[i+1] + di[i+1];
		b[i+1] = -lo[i+1]*b[i] + b[i+1];
		
		lo[i+1] = lo1;
		if(i != n-3) di[i+1] = di1;
	}
    // Multiply a row by the negative value of the entry in upper diagonal and add it to corresponding entry in diagonal. In other words, make upper diagoanl 0.
	for (int i = n-3; i >= 0; i--)
	{
		up1 = -di[i+1]*up[i] + up[i];
		b1 = -up[i]*b[i+1] + b[i];
		
		up[i] = up1;
		b[i] = b1;
	}
	return b;
}