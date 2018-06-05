void subMatrizes(int n, double matA[n][n], double matB[n][n], double matC[n][n])
{
	int i, j;

	for ( i = 0; i < n; i++ )
	{  
		for (j = 0; j < n; j++)
		{
			matC[i][j] = matA[i][j] - matB[i][j];
		}
	}
}

void multMatrizes(int n1, int n2, int n3,double matA[n1][n2], double matB[n2][n3], double matC[n1][n3])
{
	int i, j, k;

	for ( i = 0; i < n1; i++ )
	{  
		double matAik = matA[i][0];
		for ( j = 0; j < n3; j++ )
			matC[i][j] = matAik * matB[0][j];
		for ( k = 1; k < n2; k++ )
		{  
			matAik = matA[i][k];
			for ( j = 0; j < n3; j++ )
				matC[i][j] += matAik * matB[k][j];
		}
	}
}

void multMatrizVetor(int n1, int n2,double matA[n1][n2], double vetB[n2], double matC[n2])
{
	int i, j;

	for ( i = 0; i < n1; i++ )
	{  
		matC[i] = 0.0;
		for (j = 0; j < n2; j++)
		{
			matC[i] += matA[i][j] * vetB[j];
		}
	}
}



void inversaMatriz(int ordem, double a[ordem][ordem], double inv [ordem][ordem]){

	int i,j,k;

	int auxI,auxI1,auxI2;

	
	double id [ordem][ordem];

	double pivo, auxD;

	for (i = 0; i < ordem; i++)
	{
		for (j = 0; j < ordem; j++)
		{
			if (i == j)
			{
				id[i][j] = 1.0;
			}
			else
			{
				id[i][j] = 0.0;
			}
		}
	}

	for (i = 0; i < ordem; i++)
	{
		
		for(j = i+1; j < ordem; j++){

			pivo = a[i][i];
			auxD = a[j][i]/pivo;
			if (fabs(a[j][i]) > pow(10,-6))
			{
				for(k = i; k < ordem; k++)
				{
					a[j][k] = a[j][k] - a[i][k]*auxD;
					if (fabs(a[j][k]) < pow(10,-6))
					{
						a[j][k] = 0.0;
					}
				}

				a[j][i] = auxD;

			} else
			{
				a[j][i] = 0.0;
			}

		}
	}


	
	for(i = 0; i < ordem; i++)
	{
		for (j = 0; j < ordem; j++)
		{

			inv[i][j] = 0.0;
		}

	}

	for(i = 0; i < ordem - 1; i++)
	{
		for (j = i+1; j < ordem; j++)
		{
			for(k = 0; k < ordem - 1; k++)
			{
				id[j][k] = id[j][k] - a[j][i]*id[i][k];
				if (fabs(id[j][k]) < pow(10,-6))
				{
					id[j][k] = 0.0;
				}
			}
		}
	}

	for(i = ordem - 1; i >= 0; i--)
	{
		

		for (k = 0; k < ordem; k++)
		{
			auxD = 0.0;
			for(j = i; j < ordem; j++)
			{
				auxD = auxD + a[i][j]*inv[j][k];
			}

			inv[i][k] = (id[i][k]-auxD)/a[i][i];
		}

	}


}