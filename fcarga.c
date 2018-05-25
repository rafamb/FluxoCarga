#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "lerCDF.h"
#include "consPQ.h"
#include "jacobiana.h"
#include "gaussprof.h"
#include "qlim.h"


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
	// INICIALIZACAO DO VETOR DE LIGACOES DAS BARRAS
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


		if (nIteracoes >= 10)
		{
			printf("NUMERO DE ITERACOES MAXIMO ATINGIDO, PROVAVELMENTE O SISTEMA NAO CONVERGE!\n");
			goto fimProgama;
		}


		int continua = 0;


		//
		// VERIFICA QUAIS BARRAS QUE VIERAM A SER TORNAR PQ PODEM VOLTAR A SER PV
		//
		qLimInicio(nB, &nPQ, barras,ligacoes,&listaPQ);


		//
		// CRIACAO DO VETOR DE MISMATCHES
		//
		double deltaP_Q [2*nPQ+nPV];

		int i = 0;
		lista * l;


		//
		// CALCULO DO VETOR DE MISMATCHES E ANALISE SE ESTA ABAIXO DO ERRO
		//
		l = listaPQPV.prox;
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

		l = listaPQ.prox;
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
		double jac [2*nPQ+nPV][2*nPQ+nPV];		
		double x [2*nPQ+nPV];
		
		//
		// CALCULO DOS ELEMENTOS DA MATRIZ JACOBIANA
		//
		consSistema(2*nPQ+nPV,jac,barras,ligacoes,listaPQPV,listaPQ);

		//
		// RESOLUCAO DO SISTEMA LINEAR PELO METODO DE GAUSS, COM O VETOR RESPOSTA SENDO X
		//
		gauss_parcial(2*nPQ+nPV,jac,deltaP_Q,x);

		lista * li;
		li = listaPQPV.prox;
		i = 0;

		//
		// ATUALIZACAO DOS VALORES DE Vk E THETAk
		//
		while(li != NULL){
			barras[li->m].theta = barras[li->m].theta + x[i];
			i++;
			li = li->prox;
		}
		li = listaPQ.prox;
		while(li != NULL){
			barras[li->m].v = barras[li->m].v + x[i];
			i++;
			li = li->prox;
		}


		//
		// VERIFICA QUAIS BARRAS DE GERACAO VIOLARAM OS LIMITES DE GERACAO DE POTENCIA REATIVA,
		// TRANSFORMANDO-AS ENTAO EM BARRAS PQ
		//
		qLimFinal(nB, &nPQ, barras,ligacoes,&listaPQ);
		

		nIteracoes++;

	}


	fimProgama:

	liberarMemoriaLigacoes(ligacoes,nB);
	liberarLista(&listaPQPV);
	liberarLista(&listaPQ);
	
	return 0;
}