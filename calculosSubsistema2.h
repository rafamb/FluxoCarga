void calculoPRef(int ref, barra barras [], ligacao ligacoes [])
{
	int k = ref;

	barras[k].pg = 0;

	int m;
	double tKM;
	ligacao * lig = ligacoes[k].prox;
	while(lig != NULL){
		m = lig->j ;

		tKM = barras[k].theta - barras[m].theta;

		if (lig->tBarra == TAP)
		{
			barras[k].pg += (lig->info->g * (1/pow(lig->info->tap,2)))*pow(barras[k].v,2) - (1/lig->info->tap) * barras[k].v * barras[m].v * (lig->info->g * cos(tKM - lig->info->phi) + lig->info->b * sin(tKM - lig->info->phi));
		}
		else 
		{
			barras[k].pg += (lig->info->g)*pow(barras[k].v,2) - (1/lig->info->tap) * barras[k].v * barras[m].v * (lig->info->g * cos(tKM + lig->info->phi) + lig->info->b * sin(tKM + lig->info->phi));
		}

		lig = lig->prox;
	}

	barras[k].pg += -barras[k].gsh*pow(barras[k].v,2) + barras[k].pc;
}

void calculoQRefPV(int nB, barra barras [], ligacao ligacoes [])
{

	int k;
	double tKM;

	for(k = 0; k < nB; k++)
	{
		if (barras[k].tipo == 3 || barras[k].tipo == 2)
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

		}

	}

}


void calculosSubsistema2(int nB, int ref, barra barras [], ligacao ligacoes [])
{

	//
	// CALCULO DE Pg DA BARRA DE REFERENCIA
	//
	calculoPRef(ref,barras,ligacoes);

	//
	// CALCULO DE Qg PARA AS BARRAS DE REFERENCIA E PV
	//
	calculoQRefPV(nB,barras,ligacoes);
	
}