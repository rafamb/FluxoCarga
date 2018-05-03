#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "lerCDF.h"
#include "gauss.h"


//
// CONSP: FUNCAO QUE RETORNA O VALOR DE Pk
//
double consP(int k,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	int m;
	ligacao * lig = &ligacoes[k];
	while(lig != NULL){
		m = lig->j ;
		tKM = barras[k].theta - barras[m].theta;

		expr += barras[k].v * (barras[m].v*(lig->g*cos(tKM) + lig->b*sin(tKM)));

		lig = lig->prox;
	}


	return expr;
}


//
// CONSQ: FUNCAO QUE RETORNA O VALOR DE Qk
//
double consQ(int k,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	int m;
	ligacao * lig = &ligacoes[k];
	while(lig != NULL){
		m = lig->j ;
		tKM = barras[k].theta - barras[m].theta;

		expr += barras[k].v * (barras[m].v*(lig->g*sin(tKM) - lig->b*cos(tKM)));

		lig = lig->prox;
	}


	return expr;
}


//
// FINDLIGACAO: RETORNA A LIGACAO DA BARRA K PARA M
//
ligacao * findLigacao(int k,int m, ligacao ligacoes[]){
	ligacao * lig = &ligacoes[k];
	while(lig != NULL){
		if (lig->j == m)
		{
			return lig;
		}
		lig = lig->prox;

	}
	return NULL;

}


//
// MATRIZH: RETORNA O VALOR CALCULADO PARA Hkm
//
double matrizH(int k,int m,barra barras [],ligacao ligacoes [],ligacao * lig){
	if(m != k){
		return  (barras[k].v*barras[m].v*(lig->g*sin(barras[k].theta-barras[m].theta) - lig->b*cos(barras[k].theta-barras[m].theta)));
	}
	else{
		return -consQ(k,barras,ligacoes) - lig->b*pow(barras[k].v,2);
	}
}


//
// MATRIZN: RETORNA O VALOR CALCULADO PARA Nkm
//
double matrizN(int k,int m,barra barras [],ligacao ligacoes [],ligacao * lig){
	if(m != k){
		return  (barras[k].v*(lig->g*cos(barras[k].theta-barras[m].theta) + lig->b*sin(barras[k].theta-barras[m].theta)));
	}
	else{
		return pow(barras[k].v,-1)*(consP(k,barras,ligacoes) + pow(barras[k].v,2)*lig->g) ;
	}
}


//
// MATRIZM: RETORNA O VALOR CALCULADO PARA Mkm
//
double matrizM(int k,int m,barra barras [],ligacao ligacoes [],ligacao * lig){
	if(m != k){
		return  (-barras[k].v*barras[m].v*(lig->g*cos(barras[k].theta-barras[m].theta) + lig->b*sin(barras[k].theta-barras[m].theta)));
	}
	else{
		return consP(k,barras,ligacoes) - lig->g*pow(barras[k].v,2);
	}
}


//
// MATRIZL: RETORNA O VALOR CALCULADO PARA Lkm
//
double matrizL(int k,int m,barra barras [],ligacao ligacoes [],ligacao * lig){
	if(m != k){
		return  (barras[k].v*(lig->g*sin(barras[k].theta-barras[m].theta) - lig->b*cos(barras[k].theta-barras[m].theta)));
	}
	else{
		return pow(barras[k].v,-1)*(consQ(k,barras,ligacoes) - pow(barras[k].v,2)*lig->b) ;
	}
}
	
//
// CONSSISTEMA: CONSTROI O SISTEMA LINEAR A SER RESOLVIDO PELO METODO DE GAUSS
//
void consSistema(int nlinhas,int ncolunas, double jac [nlinhas][ncolunas], barra barras [], ligacao ligacoes [], double deltaP_Q [], lista listaPQPV, lista listaPQ){
	int i,j;
	lista * li, * lj;
	ligacao * lig;
	

	i = 0;


	// PARA CADA BARRA (li->m => k) PQ OU PV
	li = listaPQPV.prox;
	while(li != NULL){
		j = 0;

		// PARA CADA BARRA (lj->m => m) PQ OU PV
		lj = listaPQPV.prox;
		while(lj != NULL){

			//ENCONTRA LIGACAO ENTRE k e m
			lig = findLigacao(li->m,lj->m,ligacoes);

			//SE HA LIGACAO CALCULA Hkm
			if (lig != NULL)
				jac[i][j] = -matrizH(li->m,lj->m,barras,ligacoes,lig);
			else
				jac[i][j] = 0.0;

			lj = lj->prox;
			j++;
		}

		// PARA CADA BARRA (lj->m => m) PQ
		lj = listaPQ.prox;
		while(lj != NULL){

			//ENCONTRA LIGACAO ENTRE k e m
			lig = findLigacao(li->m,lj->m,ligacoes);

			//SE HA LIGACAO CALCULA Nkm
			if (lig != NULL)
				jac[i][j] = -matrizN(li->m,lj->m,barras,ligacoes,lig);
			else
				jac[i][j] = 0.0;

			lj = lj->prox;
			j++;
		}

		// SETA NA ULTIMA COLUNA DO SISTEMA O ELEMENTO deltaP_Q
		jac [i][j] = deltaP_Q[i];
		li = li->prox;
		i++;
	}


	// PARA CADA BARRA (li->m => k) PQ
	li = listaPQ.prox;
	while(li != NULL){
		j = 0;

		// PARA CADA BARRA (lj->m => m) PQ OU PV
		lj = listaPQPV.prox;
		while(lj != NULL){

			//ENCONTRA LIGACAO ENTRE k e m
			lig = findLigacao(li->m,lj->m,ligacoes);

			//SE HA LIGACAO CALCULA Mkm
			if (lig != NULL)
				jac[i][j] = -matrizM(li->m,lj->m,barras,ligacoes,lig);
			else
				jac[i][j] = 0.0;

			lj = lj->prox;
			j++;
		}

		// PARA CADA BARRA (lj->m => m) PQ
		lj = listaPQ.prox;
		while(lj != NULL){

			//ENCONTRA LIGACAO ENTRE k e m
			lig = findLigacao(li->m,lj->m,ligacoes);

			//SE HA LIGACAO CALCULA Lkm
			if (lig != NULL)
				jac[i][j] = -matrizL(li->m,lj->m,barras,ligacoes,lig);
			else
				jac[i][j] = 0.0;

			lj = lj->prox;
			j++;
		}

		// SETA NA ULTIMA COLUNA DO SISTEMA O ELEMENTO deltaP_Q
		jac [i][j] = deltaP_Q[i];
		li = li->prox;
		i++;
	}

}


int main(void)
{
	char arquivo[]="ieee14cdf.txt";
	FILE *arq;
	double baseMVA;
	int nB;


	//
	// CARREGANDO PONTEIRO PARA O ARQUIVO, LE A BASEMVA E DEPOIS CALCULA O NUMERO DE BARRAS (nB)
	//
	arq = carregarArquivo(arquivo);
	baseMVA = carregarBaseMVA(arq);
	nB = carregarnB(arq);


	//
	// CRIADO UM VETOR DE NB BARRAS, QUE VAI GUARDAR AS INFORMACOES LIDAS DAS BARRAS,
	// ALEM DE UM VETOR DE NB LIGACOES, QUE VAI ARMAZENAR AS INFORMACOES DAS LIGACOES DE CADA BARRA
	//
	barra barras[nB] ;
	ligacao ligacoes[nB] ;

	//
	// DECLARADOS OS VALORES NPQ, NPV E REF, QUE SAO O NUMERO DE BARRAS PQ, NUMERO DE BARRAS PV E 
	// O NUMERO DA BARRA DE REFERENCIA, RESPECTIVAMENTE
	//
	int nPQ=0;
	int nPV=0;
	int ref;

	//
	// LISTAPQPV: LISTA DAS BARRAS PQ E PV, O CAMPO LISTAPQPV.PROX APONTA PARA O PRIMEIRO ELEMENTO DA LISTA
	// E O CAMPO LISTAPQPV.ANT APONTA PARA O ULTIMO ELEMENTO DA LISTA
	//
	lista listaPQPV;

	listaPQPV.prox = NULL;
	listaPQPV.ant = NULL;

	//
	// LISTAPQ: LISTA DAS BARRAS PQ, O CAMPO LISTAPQ.PROX APONTA PARA O PRIMEIRO ELEMENTO DA LISTA
	// E O CAMPO LISTAPQ.ANT APONTA PARA O ULTIMO ELEMENTO DA LISTA
	//
	lista listaPQ;

	listaPQ.prox = NULL;
	listaPQ.ant = NULL;
	
	//
	// LE AS INFORMACOES DAS BARRAS DO ARQUIVO
	//
	carregarBarras(arq,barras,baseMVA,&nPQ,&nPV,&ref,&listaPQPV,&listaPQ);


	//
	// APLICA A SOLUCAO INICIAL (V0 = 1.0 E THETA0 = THETA_REF) PARA AS BARRAS PQ (V E THETA) E PV (THETA)
	//
	solucaoInicial(barras,nB,ref);


	//
	// INICIALIZACAO DO VETOR DE LIGACOES DAS BARRAS, ONDE O PRIMEIRO ELEMENTO DIZ RESPEITO A UMA LIGACAO FICTICIA
	// PARA A PROPRIA BARRA, APENAS PARA MAIOR FACILIDADE DE OBTENCAO DOS VALORES G E B QUANDO NECESSARIO CALCULAR
	// G11, B22 POR EXEMPLO
	//
	inicializarLigacoes(barras,ligacoes,nB);

	//
	// LE AS INFORMACOES DAS LIGACOES DO ARQUIVO
	//
	carregarLigacoes(arq,barras,ligacoes);


	fclose(arq);


	//
	// ERRO DE 10^-4
	//
	double erro = pow(10,-4);

	int nIteracoes = 0;



	while(1){
		
		//printf("%d\n", nIteracoes);

		if (nIteracoes >= 10)
		{
			printf("NUMERO DE ITERACOES MAXIMO ATINGIDO, PROVAVELMENTE O SISTEMA NAO CONVERGE!\n");
			break;
		}


		int continua = 0;


		double deltaP_Q [2*nPQ+nPV];

		int i = 0;
		lista * l;
		ligacao * lig;


		l = listaPQPV.prox;
		while(l != NULL){

			int k = l->m;

			deltaP_Q[i] = consP(k,barras,ligacoes) - (barras[k].pg - barras[k].pc);
			

			lig = ligacoes[l->m].prox;

			while(lig != NULL){
				int m = lig->j;

				double tKM = barras[k].theta - barras[m].theta;
				if (lig->info->tipo != 0)
				{
					deltaP_Q[i] -= (lig->info->g * pow(lig->info->tap,2) * pow(barras[k].v,2) - lig->info->tap * barras[k].v * barras[m].v * (lig->info->g * cos(tKM) + lig->info->b * sin(tKM)));
				} else {
					deltaP_Q[i] -= (lig->info->g * pow(barras[k].v,2) - lig->info->tap * barras[k].v * barras[m].v * (lig->info->g * cos(tKM) + lig->info->b * sin(tKM)));
				}


				lig = lig->prox;
			}

			if (continua == 0 && (fabs(deltaP_Q[i]) > erro ))
			{
				continua = 1;
			}

			l = l->prox;
			i++;
		}

		l = listaPQ.prox;
		while(l != NULL){
			int k = l->m;

			deltaP_Q[i] = consQ(l->m,barras,ligacoes) - (barras[l->m].qg - barras[l->m].qc);

			lig = ligacoes[l->m].prox;

			while(lig != NULL){
				int m = lig->j;

				double tKM = barras[k].theta - barras[m].theta;
				if (lig->info->tipo != 0)
				{
					deltaP_Q[i] -= ( -(lig->info->b * pow(lig->info->tap,2) + lig->info->bsh) * pow(barras[k].v,2) + barras[k].v * barras[m].v * (lig->info->b * cos(tKM) - lig->info->g * sin(tKM)));
				} else {
					deltaP_Q[i] -= ( -(lig->info->b) * pow(barras[k].v,2) + barras[k].v * barras[m].v * (lig->info->b * cos(tKM) - lig->info->g * sin(tKM)));
				}

				lig = lig->prox;
			}

			printf("%lf\n", deltaP_Q[i]);

			if (continua == 0 && (fabs(deltaP_Q[i]) > erro ))
			{
				continua = 1;
			}

			l = l->prox;
			i++;
		}



		if (continua == 0)
		{
			break;
		}

		double sistema [2*nPQ+nPV][2*nPQ+nPV+1];
		double x [2*nPQ+nPV];
		consSistema(2*nPQ+nPV,2*nPQ+nPV+1,sistema,barras,ligacoes,deltaP_Q, listaPQPV,listaPQ);


		gauss_parcial(2*nPQ+nPV,2*nPQ+nPV+1,sistema,x);
		
		printf("\n");

		lista * li;
		li = listaPQPV.prox;
		i = 0;

		printf("VETOR X\n");
		while(li != NULL){
			printf("%lf\n", x[i]);
			barras[li->m].theta = barras[li->m].theta + x[i];
			i++;
			li = li->prox;
		}
		li = listaPQ.prox;
		while(li != NULL){
			printf("%lf\n", x[i]);
			barras[li->m].v = barras[li->m].v + x[i];
			i++;
			li = li->prox;
		}

		break;
	
		

		nIteracoes++;

	}

	barras[ref].p = consP(ref,barras,ligacoes);

	

	int i;
	for(i=0; i<nB; i++){
		if (barras[i].tipo == 3)
		{
			barras[i].q = consQ(i,barras,ligacoes);
		} else if (barras[i].tipo == 2 )
		{
			barras[i].q = consQ(i,barras,ligacoes);
			barras[i].p = barras[i].pg - barras[i].pc;
		} 
		else 
		{
			barras[i].p = barras[i].pg - barras[i].pc;
			barras[i].q = barras[i].qg - barras[i].qc;
		}
	}

	for (i = 0; i < nB; i++)
	{
		//printf("%d: V: %lf Theta: %lf P: %lf Q: %lf\n", i+1,barras[i].v,barras[i].theta,barras[i].p,barras[i].q);
	}




	liberarMemoriaLigacoes(ligacoes,nB);
	liberarLista(&listaPQPV);
	liberarLista(&listaPQ);

	//printLigacoes(ligacoes,nB);	
	
	return 0;
}