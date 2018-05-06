
void gauss_parcial(int nlinhas, int ncolunas, double sistema[nlinhas][ncolunas], double x [nlinhas]){
	int i, j, k, p, auxId;
	double aux, pivo;

	int indices [nlinhas];

	for (i = 0; i < nlinhas; i++)
	{
		indices[i] = i;
	}

	for (k = 0; k < nlinhas - 1; k++){
		p = k;
		aux = fabs(sistema[k][k]);
		for (i = k+1; i < nlinhas; i++)
		{
			if (fabs(sistema[i][k]) > aux)
			{
				p = i;
				aux = fabs(sistema[i][k]);
			}
		}

		if (p != k)
		{
			auxId = indices[k];
			indices[k] = indices[p];
			indices[p] = auxId;
			for (j = k; j <= ncolunas; j++){
				aux = sistema[p][i];
				sistema[p][j] = sistema[k][j];
				sistema[k][j] = aux;
			}
		}

		for (i = k+1; i < nlinhas; i++){
			aux = sistema[i][k]/sistema[k][k];
			if (fabs(sistema[i][k]) > pow(10,-5))
			{
				for (j = k; j <= nlinhas ; j++){
					sistema[i][j] = sistema[i][j] - aux*sistema[k][j];
					if (fabs(sistema[i][j]) < pow(10,-5))
					{
						sistema[i][j] = 0.0;
					}
				}
			} else {
				sistema[i][k] = 0.0;
			}
			
			
		}
	}
	

	x[nlinhas-1]=sistema[nlinhas-1][ncolunas-1]/sistema[nlinhas-1][nlinhas-1];

	for(k = nlinhas-2;k >= 0;k--)
	{
		aux = 0.0;
		for(j = k+1;j < nlinhas; j++)
		{
			aux=aux + sistema[k][j]*x[j];
		}
		x[k]=(sistema[k][ncolunas-1]-aux)/sistema[k][k];
	}

}