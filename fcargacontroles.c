#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "estruturas.h"
#include "lerCDF.h"
#include "fc.h"
#include "metodologiasControle.h"


#define MAXQ 1
#define VSF 2



int main(int argc, char *argv[])
{
	char arquivo[50];

	int print = 0;

	int metodologia = MAXQ;

	switch(argc){
		case 1:
			strcpy(arquivo,"ieee14cdf.txt");
			break;
		default:
			strcpy(arquivo,argv[1]);
			break;
	}
	if (argc == 3)
	{
		if (strcmp(argv[2],"vsf") == 0)
		{
			printf("vsf\n");
			metodologia = VSF;
		}
	}

	
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
	int nPQ = 0;
	int nPV = 0;
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

	//
	// NUMERO MAXIMO DE ITARACOES
	//
	int nMax = 10;

	int i = 0;

	while(1)
	{
		printf("%d\n", i);

		//
		// RESOLUCAO DO FLUXO DE CARGA
		//
		fc(nB, ref, &nPQ, &nPV, erro, nMax, barras, ligacoes, &listaPQ, &listaPQPV);

		int barraViolada = -1;

		double dVc [nB];

		//
		// CALCULO DE dVc
		//
		desvioTensao(nB, barras, &barraViolada, dVc);

		//
		// SE NAO HOUVER BARRAS DE CARGA COM VIOLACOES DE TENSAO O ALGORITMO SE ENCERRA
		//
		if (barraViolada == -1)
		{
			printf("SEM MAIS BARRAS VIOLADAS\n");
			break;
		}

		printf("Barra %d\n", barraViolada+1);

		double jqv [nB][nB];
		
		//
		// CALCULO DA MATRIZ JQV
		//
		calcJqv(nB, barras, ligacoes,listaPQPV, jqv);

		double dQ [nB];

		//
		// CALCULO DE dQ
		//
		calcDeltaQg(nB, barras,jqv, dVc, dQ);

		double dVg [nB];

		int j;

		if (metodologia == MAXQ)
		{
			j = maxQ( nB, barras, dVc , dQ , jqv , dVg );
		}
		else if (metodologia == VSF)
		{

			j = vsf( nB, barras, barraViolada , dQ , jqv , dVg );
		}

		atualizarV(barras, j, dVg);

		printf("NOVO VALOR V%d = %lf\n", j+1,barras[j].v);

		i++;
		if (i == 4)
		{
			break;
		}

	}

	liberarMemoriaLigacoes(ligacoes,nB);
	liberarLista(&listaPQPV);
	liberarLista(&listaPQ);
	
	return 0;
}