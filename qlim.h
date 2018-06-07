
void qLimInicio(int nB, int * nPQ, int * nPV, barra barras [], ligacao ligacoes [], lista * listaPQ)
{
	int k;
	lista * li;
	li = listaPQ->prox;
	while(li != NULL){
		k = li->m;
		if (barras[k].tipo == 4)
		{
			if ((barras[k].v <= barras[k].vEsp && barras[k].qg == barras[k].qgMin) || (barras[k].v >= barras[k].vEsp && barras[k].qg == barras[k].qgMax))
			{
				barras[k].v = barras[k].vEsp;
				barras[k].tipo = 2;
				removerLista(listaPQ,k);
				*nPQ = *nPQ - 1;
				*nPV = *nPV + 1;
			}
		}

		li = li->prox;
	}
}


void qLimFinal(int nB, int * nPQ, int * nPV, barra barras [], ligacao ligacoes [], lista * listaPQ)
{
	int k;
	double tKM;


	for(k = 0; k < nB; k++){
		if (barras[k].tipo == 4)
		{
			if ((barras[k].v <= barras[k].vEsp && barras[k].qg == barras[k].qgMin) || (barras[k].v >= barras[k].vEsp && barras[k].qg == barras[k].qgMax))
			{
				//barras[k].v = barras[k].vEsp;
				barras[k].tipo = 2;
				removerLista(listaPQ,k);
				*nPQ = *nPQ - 1;
				*nPV = *nPV + 1;
			}
		}

		if (barras[k].tipo == 2)
		{
			barras[k].qg = 0;

			int m;
			ligacao * lig = ligacoes[k].prox;
			while(lig != NULL){
				m = lig->j ;

				tKM = barras[k].theta - barras[m].theta;

				if (lig->tBarra == TAP)
				{
					barras[k].qg += -(lig->info->b * (1/pow(lig->info->tap,2)) + lig->info->bsh)*pow(barras[k].v,2) + (1/lig->info->tap) * barras[k].v * barras[m].v * (lig->info->b * cos(tKM - lig->info->phi) - lig->info->g * sin(tKM - lig->info->phi));
				}
				else 
				{
					barras[k].qg += -(lig->info->b)*pow(barras[k].v,2) + (1/lig->info->tap) * barras[k].v * barras[m].v * (lig->info->b * cos(tKM + lig->info->phi) - lig->info->g * sin(tKM + lig->info->phi));
				}

				lig = lig->prox;
			}

			barras[k].qg += -barras[k].bsh*pow(barras[k].v,2) + barras[k].qc;

			if (barras[k].qg < barras[k].qgMin)
			{
				barras[k].qg = barras[k].qgMin;
				barras[k].tipo = 4;
				inserirLista(listaPQ,k);
				*nPQ = *nPQ + 1;
				*nPV = *nPV - 1;
			} else if (barras[k].qg > barras[k].qgMax)
			{
				barras[k].qg = barras[k].qgMax;
				barras[k].tipo = 4;
				inserirLista(listaPQ,k);
				*nPQ = *nPQ + 1;
				*nPV = *nPV - 1;
			}
		}

	}

}