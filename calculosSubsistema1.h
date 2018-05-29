#include "consPQ.h"
#include "jacobiana.h"
#include "gaussprof.h"
#include "qlim.h"

int calculosSubsistema1(int nB, int * nPQ, int * nPV, double erro, int nMAX, barra barras [], ligacao ligacoes [], lista * listaPQ, lista * listaPQPV)
{

	int nIteracoes = 0;


	while(1){


		if (nIteracoes >= nMAX)
		{
			printf("NUMERO DE ITERACOES MAXIMO ATINGIDO, PROVAVELMENTE O SISTEMA NAO CONVERGE!\n");
			return 0;
		}


		int continua = 0;


		//
		// VERIFICA QUAIS BARRAS QUE VIERAM A SER TORNAR PQ PODEM VOLTAR A SER PV
		//
		qLimInicio(nB, nPQ, barras,ligacoes,listaPQ);


		//
		// CRIACAO DO VETOR DE MISMATCHES
		//
		double deltaP_Q [2*(*nPQ)+(*nPV)];

		int i = 0;
		lista * l;


		//
		// CALCULO DO VETOR DE MISMATCHES E ANALISE SE ESTA ABAIXO DO ERRO
		//
		l = listaPQPV->prox;
		while(l != NULL){

			int k = l->m;

			deltaP_Q[i] = -consP(k,barras,ligacoes);


			if (continua == 0 && (fabs(deltaP_Q[i]) > erro ))
			{
				continua = 1;
			}

			l = l->prox;
			i++;
		}

		l = listaPQ->prox;
		while(l != NULL){
			int k = l->m;

			deltaP_Q[i] = -consQ(k,barras,ligacoes);
			

			if (continua == 0 && (fabs(deltaP_Q[i]) > erro ))
			{
				continua = 1;
			}

			l = l->prox;
			i++;
		}


		//
		// SE TODOS OS VALORES ESTAO MENORES QUE O ERRO O METODO DE NEWTON E ENCERRADO
		//
		if (continua == 0)
		{
			break;
		}


		//
		// CRIACAO DA MATRIZ JACOBIANA E DO VETOR X
		//
		double jac [2*(*nPQ)+(*nPV)][2*(*nPQ)+(*nPV)];		
		double x [2*(*nPQ)+(*nPV)];
		
		//
		// CALCULO DOS ELEMENTOS DA MATRIZ JACOBIANA
		//
		consSistema(2*(*nPQ)+(*nPV),jac,barras,ligacoes,*listaPQPV,*listaPQ);

		//
		// RESOLUCAO DO SISTEMA LINEAR PELO METODO DE GAUSS, COM O VETOR RESPOSTA SENDO X
		//
		gauss_parcial(2*(*nPQ)+(*nPV),jac,deltaP_Q,x);

		lista * li;
		li = listaPQPV->prox;
		i = 0;

		//
		// ATUALIZACAO DOS VALORES DE Vk E THETAk
		//
		while(li != NULL){
			barras[li->m].theta = barras[li->m].theta + x[i];
			i++;
			li = li->prox;
		}
		li = listaPQ->prox;
		while(li != NULL){
			barras[li->m].v = barras[li->m].v + x[i];
			i++;
			li = li->prox;
		}


		//
		// VERIFICA QUAIS BARRAS DE GERACAO VIOLARAM OS LIMITES DE GERACAO DE POTENCIA REATIVA,
		// TRANSFORMANDO-AS ENTAO EM BARRAS PQ
		//
		qLimFinal(nB, nPQ, barras,ligacoes,listaPQ);
		

		nIteracoes++;

	}

	return 1;

}