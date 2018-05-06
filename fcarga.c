#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "lerCDF.h"
#include "gaussprof.h"


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


//FUNCAO AUXILIAR PARA CALCULO DO ELEMENTO Lkk
double auxL(int k,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	int m;
	ligacao * lig = ligacoes[k].prox;
	while(lig != NULL){
		m = lig->j ;

		tKM = barras[k].theta - barras[m].theta;

		expr +=2*(lig->info->bsh*barras[k].v) - (barras[m].v*(lig->g*sin(tKM) - lig->b*cos(tKM)));

		lig = lig->prox;
	}


	return expr;
}

	
//
// CONSSISTEMA: CONSTROI O SISTEMA LINEAR A SER RESOLVIDO PELO METODO DE GAUSS
//
void consSistema(int nlinhas,double jac [nlinhas][nlinhas], barra barras [], ligacao ligacoes [], lista listaPQPV, lista listaPQ){
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
			int k = li->m;
			int m = lj->m;

			jac[i][j] = 0.0;

			
			//AQUI BUSCA TODAS AS LIGACOES ENTRE K E M 
			//APOS ISSO CALCULA-SE Hkm
			if (k == m){

				jac[i][j] = -matrizH(k,k,barras,ligacoes,&ligacoes[k]);

			} else {

				lig = &ligacoes[k];
				while(lig != NULL){
					if (lig->j == m)
					{
						jac[i][j] += -matrizH(li->m,lj->m,barras,ligacoes,lig);
					}
					lig = lig->prox;

				}
			}

			lj = lj->prox;
			j++;
		}

		// PARA CADA BARRA (lj->m => m) PQ
		lj = listaPQ.prox;
		while(lj != NULL){

			int k = li->m;
			int m = lj->m;

			jac[i][j] = 0.0;

			
			//AQUI BUSCA TODAS AS LIGACOES ENTRE K E M 
			//APOS ISSO CALCULA-SE Nkm
			if (k == m){

				jac[i][j] = -matrizN(k,k,barras,ligacoes,&ligacoes[k]);

			} else {

				lig = &ligacoes[k];
				while(lig != NULL){
					if (lig->j == m)
					{
						jac[i][j] += -matrizN(li->m,lj->m,barras,ligacoes,lig);
					}
					lig = lig->prox;

				}
			}

			lj = lj->prox;
			j++;
		}

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

			int k = li->m;
			int m = lj->m;

			jac[i][j] = 0.0;

			
			//AQUI BUSCA TODAS AS LIGACOES ENTRE K E M 
			//APOS ISSO CALCULA-SE Mkm
			if (k == m){

				jac[i][j] = -matrizM(k,k,barras,ligacoes,&ligacoes[k]);

			} else {

				lig = &ligacoes[k];
				while(lig != NULL){
					if (lig->j == m)
					{
						jac[i][j] += -(-barras[k].v*barras[m].v*(lig->g*cos(barras[k].theta-barras[m].theta) + lig->b*sin(barras[k].theta-barras[m].theta)));
					}
					lig = lig->prox;

				}
			}

			lj = lj->prox;
			j++;
		}

		// PARA CADA BARRA (lj->m => m) PQ
		lj = listaPQ.prox;
		while(lj != NULL){

			int k = li->m;
			int m = lj->m;

			jac[i][j] = 0.0;

			
			//AQUI BUSCA TODAS AS LIGACOES ENTRE K E M 
			//APOS ISSO CALCULA-SE Lkm
			if (k == m){

				jac[i][j] = (2*(ligacoes[k].b*barras[k].v) + auxL(k,barras,ligacoes));

			} else {

				lig = &ligacoes[k];
				while(lig != NULL){
					if (lig->j == m)
					{
						jac[i][j] += -(barras[k].v*(lig->g*sin(barras[k].theta-barras[m].theta) - lig->b*cos(barras[k].theta-barras[m].theta)));
					}
					lig = lig->prox;

				}
			}

			lj = lj->prox;
			j++;
		}

		li = li->prox;
		i++;
	}

}

/*//
// CONSSISTEMA: CONSTROI O SISTEMA LINEAR A SER RESOLVIDO PELO METODO DE GAUSS
//
void consSistema(int nlinhas, int ncolunas,double jac [nlinhas][ncolunas], barra barras [], ligacao ligacoes [],double deltaP_Q [], lista listaPQPV, lista listaPQ){
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
			int k = li->m;
			int m = lj->m;

			jac[i][j] = 0.0;

			
			//AQUI BUSCA TODAS AS LIGACOES ENTRE K E M 
			//APOS ISSO CALCULA-SE Hkm
			if (k == m){

				jac[i][j] = -matrizH(k,k,barras,ligacoes,&ligacoes[k]);

			} else {

				lig = &ligacoes[k];
				while(lig != NULL){
					if (lig->j == m)
					{
						jac[i][j] += -matrizH(li->m,lj->m,barras,ligacoes,lig);
					}
					lig = lig->prox;

				}
			}

			lj = lj->prox;
			j++;
		}

		// PARA CADA BARRA (lj->m => m) PQ
		lj = listaPQ.prox;
		while(lj != NULL){

			int k = li->m;
			int m = lj->m;

			jac[i][j] = 0.0;

			
			//AQUI BUSCA TODAS AS LIGACOES ENTRE K E M 
			//APOS ISSO CALCULA-SE Nkm
			if (k == m){

				jac[i][j] = -matrizN(k,k,barras,ligacoes,&ligacoes[k]);

			} else {

				lig = &ligacoes[k];
				while(lig != NULL){
					if (lig->j == m)
					{
						jac[i][j] += -matrizN(li->m,lj->m,barras,ligacoes,lig);
					}
					lig = lig->prox;

				}
			}

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

			int k = li->m;
			int m = lj->m;

			jac[i][j] = 0.0;

			
			//AQUI BUSCA TODAS AS LIGACOES ENTRE K E M 
			//APOS ISSO CALCULA-SE Mkm
			if (k == m){

				jac[i][j] = -matrizM(k,k,barras,ligacoes,&ligacoes[k]);

			} else {

				lig = &ligacoes[k];
				while(lig != NULL){
					if (lig->j == m)
					{
						jac[i][j] += -(-barras[k].v*barras[m].v*(lig->g*cos(barras[k].theta-barras[m].theta) + lig->b*sin(barras[k].theta-barras[m].theta)));
					}
					lig = lig->prox;

				}
			}

			lj = lj->prox;
			j++;
		}

		// PARA CADA BARRA (lj->m => m) PQ
		lj = listaPQ.prox;
		while(lj != NULL){

			int k = li->m;
			int m = lj->m;

			jac[i][j] = 0.0;

			
			//AQUI BUSCA TODAS AS LIGACOES ENTRE K E M 
			//APOS ISSO CALCULA-SE Lkm
			if (k == m){

				jac[i][j] = (2*(ligacoes[k].b*barras[k].v) + auxL(k,barras,ligacoes));

			} else {

				lig = &ligacoes[k];
				while(lig != NULL){
					if (lig->j == m)
					{
						jac[i][j] += -(barras[k].v*(lig->g*sin(barras[k].theta-barras[m].theta) - lig->b*cos(barras[k].theta-barras[m].theta)));
					}
					lig = lig->prox;

				}
			}

			lj = lj->prox;
			j++;
		}

		// SETA NA ULTIMA COLUNA DO SISTEMA O ELEMENTO deltaP_Q
		jac [i][j] = deltaP_Q[i];

		li = li->prox;
		i++;
	}

}*/


int main(void)
{
	char arquivo[]="ieee57cdf.txt";
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
	int convergiu = 1;



	while(1){

		if (nIteracoes >= 10)
		{
			printf("NUMERO DE ITERACOES MAXIMO ATINGIDO, PROVAVELMENTE O SISTEMA NAO CONVERGE!\n");
			convergiu = 0;
			break;
		}


		int continua = 0;


		double deltaP_Q [2*nPQ+nPV];

		int i = 0;
		lista * l;

		printf("P\n");

		l = listaPQPV.prox;
		while(l != NULL){

			int k = l->m;

			deltaP_Q[i] = -(barras[k].pg - barras[k].pc - consP(k,barras,ligacoes));

			printf("%.5lf\n", -deltaP_Q[i]);


			if (continua == 0 && (fabs(deltaP_Q[i]) > erro ))
			{
				continua = 1;
			}

			l = l->prox;
			i++;
		}

		printf("Q\n");

		l = listaPQ.prox;
		while(l != NULL){
			int k = l->m;

			deltaP_Q[i] = -((barras[k].qg - barras[k].qc) - consQ(k,barras,ligacoes));

			printf("%.5lf\n", -deltaP_Q[i]);
			

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

		double jac [2*nPQ+nPV][2*nPQ+nPV];
		//double sistema [2*nPQ+nPV][2*nPQ+nPV+1];
		double x [2*nPQ+nPV];
		//consSistema(2*nPQ+nPV,2*nPQ+nPV+1,sistema,barras,ligacoes,deltaP_Q, listaPQPV,listaPQ);
		consSistema(2*nPQ+nPV,jac,barras,ligacoes,listaPQPV,listaPQ);

		/*int j;
		for (i = 0; i < 2*nPQ+nPV; i++)
		{
			for (j = 0; j <= 2*nPQ+nPV; j++ ){
				printf("%.5lf;", sistema[i][j]);
			}
			printf("\n");
		}*/


		gauss_parcial(2*nPQ+nPV,jac,deltaP_Q,x);
		//gauss_parcial(2*nPQ+nPV,2*nPQ+nPV+1,sistema,x);
		

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

		/*for (i = 0; i < nB; i++)
		{
			printf("%d: V: %lf Theta: %lf\n", i+1,barras[i].v,barras[i].theta);
		}*/

		//break;
	
		

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
		printf("%d: V: %lf Theta: %lf P: %lf Q: %lf\n", i+1,barras[i].v,barras[i].theta,barras[i].p,barras[i].q);
	}




	liberarMemoriaLigacoes(ligacoes,nB);
	liberarLista(&listaPQPV);
	liberarLista(&listaPQ);
	
	return 0;
}