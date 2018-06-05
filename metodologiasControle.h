#include "opMatrizes.h"

void desvioTensao(int nB, barra barras [nB], int * maiorDesvio, double dVc [nB])
{
	
	int k;
	double erro = 0;
	double eK;
	for(k = 0; k < nB; k++)
	{
		if (barras[k].tipo != 3 && barras[k].tipo != 2)
		{
			
			if (barras[k].v < barras[k].vMin)
			{
				eK = barras[k].vMin - barras[k].v;
			}
			else if (barras[k].v > barras[k].vMax)
			{
				eK = barras[k].v - barras[k].vMax;
			}
			else
				eK = 0.0;

			dVc[k] = eK;

			if (eK > erro)
			{
				erro = eK;
				*maiorDesvio = k;
			}
		}
		else
			dVc[k] = 0.0;
		
	}

}

void calcJ1(int nB,barra barras [nB], ligacao ligacoes [nB], double j1 [nB-1][nB-1], lista listaPQPV)
{

	lista * li, * lj;

	int k, m;

	int i, j;

	i = 0;
	for (li = listaPQPV.prox; li != NULL; li = li->prox)
	{
		k = li->m;
		j = 0;
		for (lj = listaPQPV.prox; lj != NULL; lj = lj->prox)
		{
			m = lj->m;
			j1[i][j] = 0.0;

			if (k == m){

				j1[i][j] = calc_Hkk(k,barras,ligacoes);

			} else {

				j1[i][j] = calc_Hkm(k,m,barras,ligacoes);
			}

			j++;
		}

		i++;
	}

}

void calcJ2(int nB,barra barras [nB], ligacao ligacoes [nB], double j2 [nB-1][nB], lista listaPQPV)
{

	lista * li;

	int k, m;

	int i, j;

	i = 0;
	for (li = listaPQPV.prox; li != NULL; li = li->prox)
	{
		k = li->m;
		j=0;
		for (m = 0; m < nB; m++)
		{
			j2[i][j] = 0.0;

			if (k == m){

				j2[i][j] = (2*(barras[k].gsh*barras[k].v) + calc_Nkk(k,barras,ligacoes));

			} else {

				j2[i][j] = calc_Nkm(k,m,barras,ligacoes);
			}
		
			j++;	
		}

		i++;
	}

}

void calcJ3(int nB,barra barras [nB], ligacao ligacoes [nB], double j3 [nB][nB-1], lista listaPQPV)
{

	lista  * lj;

	int k, m;

	int i, j;

	i = 0;

	for (k = 0; k < nB; k++)
	{
		j = 0;
		for (lj = listaPQPV.prox; lj != NULL; lj = lj->prox)
		{
			m = lj->m;

			j3[i][j] = 0.0;

			if (k == m){

				j3[i][j] = calc_Mkk(k,barras,ligacoes);

			} else {

				j3[i][j] = calc_Mkm(k,m,barras,ligacoes);
			}


			j++;
		}

		i++;
	}

}

void calcJ4(int nB,barra barras [nB], ligacao ligacoes [nB], double j4 [nB][nB])
{

	int k, m;

	int i, j;

	i = 0;
	for (k = 0; k < nB; k++)
	{
		j = 0;
		for (m = 0; m < nB; m++)
		{
			j4[i][j] = 0.0;

			if (k == m){

				j4[i][j] = (2*(barras[k].bsh*barras[k].v) + calc_Lkk(k,barras,ligacoes));

			} else {

				j4[i][j] = calc_Lkm(k,m,barras,ligacoes);
			}
			
			j++;
		}

		i++;
	}

}

void calcJqv(int nB,barra barras [nB],ligacao ligacoes [nB],lista listaPQPV,double jqv [nB][nB])
{
	double j1 [nB-1][nB-1];
	double j2 [nB-1][nB];
	double j3 [nB][nB-1];
	double j4 [nB][nB];
	double invj1 [nB-1][nB-1];
	double j31 [nB][nB-1];
	double j312 [nB][nB];
	
	calcJ1(nB,barras,ligacoes,j1,listaPQPV);
	calcJ2(nB,barras,ligacoes,j2,listaPQPV);
	calcJ3(nB,barras,ligacoes,j3,listaPQPV);
	calcJ4(nB,barras,ligacoes,j4);


	//
	// CALCULANDO J1⁻¹
	//
	inversaMatriz(nB-1,j1,invj1);

	//
	// CALCULANDO J3*J1⁻¹
	//
	multMatrizes(nB, nB-1, nB-1,j3, invj1, j31);

	//
	// CALCULANDO (J3*J1⁻¹)*J2
	//
	multMatrizes(nB, nB-1, nB,j31, j2, j312);

	//
	// CALCULANDO J4 - (J3*J1⁻¹*J2)
	//
	subMatrizes(nB,j4,j312,jqv);

}

void calcDeltaQg(int nB, barra barras [nB],double jqv [nB][nB], double dVc [nB], double dQ [nB])
{
	multMatrizVetor(nB,nB, jqv, dVc, dQ);

	int i;
	for(i = 0; i < nB; i++)
	{
		if (barras[i].tipo != 3 && barras[i].tipo != 2)
		{
			dQ[i] = 0.0;
		}
	}
}

int maxDQ(int nB, barra barras [nB], double dQ [nB])
{
	int mInd = -1;

	int i;
	for (i = 0; i < nB; i++)
	{
		if (i == 1)
		{
			printf("DQ[1] %lf\n", dQ[i]);
		}
		if (mInd == -1 && (barras[i].tipo == 3 || barras[i].tipo == 2))
		{
			if (barras[i].v > barras[i].vMin && (barras[i].v < barras[i].vMax))
			{
				mInd = i;
			}
		}
		else if (dQ[i] > dQ[mInd] && (barras[i].tipo == 3 || barras[i].tipo == 2))
		{
			if (barras[i].v > barras[i].vMin && (barras[i].v < barras[i].vMax))
			{
				mInd = i;
			}
		}
		
	}

	printf("DQMAX[%d] %lf\n", mInd,dQ[mInd]);

	return mInd;

}

int maxQ(int nB, barra barras [nB], double dVc [nB], double dQ [nB], double jqv [nB][nB], double dVg [nB])
{
	int j = maxDQ(nB, barras, dQ);
	if (j == -1)
	{
		return -1;
	}

	gauss_parcial(nB,jqv,dQ,dVg);

	return j;
}

int maxS(int nB, int coluna, barra barras [nB], double s[nB][nB])
{
	int i ,j = -1;
	for (i = 0; i < nB; i ++)
	{
		
		if (j == -1 && (barras[i].tipo == 3 || barras[i].tipo == 2))
		{
			if (barras[i].v > barras[i].vMin && (barras[i].v < barras[i].vMax))
			{
				j = i;
			}
		}
		else if (s[i][coluna] > s[j][coluna])
		{
			if (barras[i].v > barras[i].vMin && (barras[i].v < barras[i].vMax))
			{
				j = i;
			}
		}
	}

	return j;
}

int vsf(int nB, barra barras [nB], int maxDVc, double dQ [nB], double jqv [nB][nB], double dVg [nB])
{
	double s [nB][nB];

	inversaMatriz(nB,jqv,s);

	int j = maxS(nB, maxDVc, barras, s);

	multMatrizVetor(nB, nB,s, dQ, dVg);

	return j;
}

void atualizarV(barra barras [], int j, double dVg [])
{
	double nV = barras[j].v + dVg[j];
	if (nV > barras[j].vMax)
	{
		barras[j].v = barras[j].vMax;
	}
	else if (nV < barras[j].vMin)
	{
		barras[j].v = barras[j].vMin;
	}
	else
	{
		barras[j].v = nV;
	}

}