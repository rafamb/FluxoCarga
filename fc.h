#include "calculosSubsistema1.h"
#include "calculosSubsistema2.h"
#include "perdas.h"

double fc(int nB, int ref, int * nPQ, int * nPV, double erro, int nMax, barra barras [], ligacao ligacoes [], lista * listaPQ, lista * listaPQPV){
	//
	// RESOLUCAO DO SUBSISTEMA1 PELO METODO DE NEWTON-RAPHSON
	//
	int convergiu = calculosSubsistema1(nB, nPQ, nPV, erro, nMax, barras,ligacoes, listaPQ, listaPQPV);

	if (convergiu == 0)
	{
		return -999;
	}

	//
	// RESOLUCAO DO SUBSISTEMA2
	//
	calculosSubsistema2(nB,ref,barras,ligacoes);

	printSolucao(nB, barras);

	double perdas = calcPerdas(nB,barras,ligacoes);

	return perdas;
}