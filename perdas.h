double calcPerdas(int nB, barra barras [], ligacao ligacoes [])
{
	double perdas = 0;
	double tKM;

	int k;

	ligacao * lig;
	for(k = 0; k < nB; k++){
		lig = ligacoes[k].prox;
		while(lig != NULL){
			if (lig->info->i == k + 1)
			{
				int m = lig->j;

				tKM = barras[k].theta - barras[m].theta;

				perdas += lig->info->g*((1/pow(lig->info->tap,2))*pow(barras[k].v,2) + pow(barras[m].v,2) - 2 * (1/lig->info->tap) *barras[k].v *barras[m].v * cos(tKM - lig->info->phi));
			}

			lig = lig->prox;
		}
	}

	return perdas;

}