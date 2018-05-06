
void gauss_parcial(int ordem, double a[ordem][ordem], double b [ordem],double x [ordem]){
	
	int ordemEqs [ordem];

	int i,j,k;

	int auxI,auxI1,auxI2;

	double temp [ordem];

	double pivo, auxD;


	for (i = 0; i < ordem; i++)
	{
		ordemEqs[i] = i;
	}

	for (i = 0; i < ordem; i++)
	{
		auxI1 = -1;
		auxD = fabs(a[i][i]);
		for (k = i+1; k < ordem; k++)
		{
			if (auxD < fabs(a[k][i]))
			{
				auxI1 = k;
				auxD = fabs(a[k][i]);
			}
		}


		if (auxI1 != -1)
		{
			auxI2 = ordemEqs[i];
			ordemEqs[i] = ordemEqs[auxI1];
			ordemEqs[auxI1] = auxI2;

			for (k = 0; k < ordem; k++)
			{
				auxD = a[i][k];
				a[i][k] = a[auxI1][k];
				a[auxI1][k] = auxD;
				
			}
		}

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
		temp[i] = b[i];

		x[i] = 0.0;

	}

	for(i = 0; i < ordem; i++)
	{
		auxI = ordemEqs[i];
		b[i] = temp[auxI];

	}

	for(i = 0; i < ordem - 1; i++)
	{
		for(j = i+1; j < ordem; j++)
		{
			b[j] = b[j] - a[j][i]*b[i];
			if (fabs(b[j]) < pow(10,-6))
			{
				b[j] = 0.0;
			}
		}
	}

	for(i = ordem - 1; i >= 0; i--)
	{
		auxD = 0.0;

		for(j = i; j < ordem; j++)
		{
			auxD = auxD + a[i][j]*x[j];
		}

		x[i] = (b[i]-auxD)/a[i][i];

	}


}