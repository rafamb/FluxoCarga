//
// CONSP: FUNCAO QUE RETORNA O VALOR DE Pk
//
double consP(int k,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	int m;
	ligacao * lig = ligacoes[k].prox;
	while(lig != NULL){
		m = lig->j ;

		tKM = barras[k].theta - barras[m].theta;

		if (lig->tBarra == TAP)
		{
			expr += (lig->info->g * (1/pow(lig->info->tap,2)))* pow(barras[k].v,2) - (1/lig->info->tap)*barras[k].v*barras[m].v*(lig->info->g*cos(tKM - lig->info->phi) + lig->info->b*sin(tKM - lig->info->phi));
		}
		else 
		{
			expr += (lig->info->g)* pow(barras[k].v,2) - (1/lig->info->tap)*barras[k].v*barras[m].v*(lig->info->g*cos(tKM + lig->info->phi) + lig->info->b*sin(tKM + lig->info->phi));
		}
		

		lig = lig->prox;
	}

	expr = barras[k].pg + barras[k].gsh*pow(barras[k].v,2) - barras[k].pc - expr;


	return expr;
}


//
// CONSQ: FUNCAO QUE RETORNA O VALOR DE Qk
//
double consQ(int k,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	int m;
	ligacao * lig = ligacoes[k].prox;
	while(lig != NULL){
		m = lig->j ;

		tKM = barras[k].theta - barras[m].theta;


		if (lig->tBarra == TAP)
		{
			expr += -(lig->info->b * (1/pow(lig->info->tap,2)) + lig->info->bsh)* pow(barras[k].v,2) + (1/lig->info->tap)*barras[k].v*barras[m].v*(lig->info->b*cos(tKM - lig->info->phi) - lig->info->g*sin(tKM - lig->info->phi));
		}
		else 
		{
			expr += -(lig->info->b)* pow(barras[k].v,2) + (1/lig->info->tap)*barras[k].v*barras[m].v*(lig->info->b*cos(tKM + lig->info->phi) - lig->info->g*sin(tKM + lig->info->phi));
		}

		

		lig = lig->prox;
	}

	expr = barras[k].qg + barras[k].bsh*pow(barras[k].v,2) - barras[k].qc - expr;


	return expr;
}