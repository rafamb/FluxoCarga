
double calc_Hkk(int k,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	int m;
	ligacao * lig = ligacoes[k].prox;
	while(lig != NULL){
		m = lig->j ;

		tKM = barras[k].theta - barras[m].theta;

		if (lig->tBarra == TAP)
		{
			expr += (1/lig->info->tap) * barras[k].v * barras[m].v * (lig->info->b * cos(tKM - lig->info->phi) - lig->info->g * sin(tKM - lig->info->phi));
		}
		else 
		{
			expr += (1/lig->info->tap) * barras[k].v * barras[m].v * (lig->info->b * cos(tKM + lig->info->phi) - lig->info->g * sin(tKM + lig->info->phi));
		}

		lig = lig->prox;
	}


	return expr;
}

double calc_Hkm(int k,int m,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	ligacao * lig = ligacoes[k].prox;

	while(lig != NULL){
		if (lig->j == m)
		{
			if (lig->tBarra == TAP)
			{
				expr += -( (1/lig->info->tap) * barras[k].v * barras[m].v * (lig->info->b * cos(lig->info->phi - (barras[k].theta-barras[m].theta)) + lig->info->g * sin(lig->info->phi - (barras[k].theta-barras[m].theta))) );
			}
			else
			{
				expr += -( (1/lig->info->tap) * barras[k].v * barras[m].v * (lig->info->b * cos(lig->info->phi + (barras[k].theta-barras[m].theta)) - lig->info->g * sin(lig->info->phi + (barras[k].theta-barras[m].theta))) );
			}
		}
		lig = lig->prox;

	}


	return expr;
}

double calc_Nkk(int k,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	int m;
	ligacao * lig = ligacoes[k].prox;
	while(lig != NULL){
		m = lig->j ;

		tKM = barras[k].theta - barras[m].theta;

		if (lig->tBarra == TAP)
		{
			expr += -2*lig->info->g*(1/pow(lig->info->tap,2))*barras[k].v + (1/lig->info->tap)*barras[m].v*( lig->info->g*cos(tKM - lig->info->phi)  +  lig->info->b*sin(tKM - lig->info->phi) );
		}
		else 
		{
			expr += -2*lig->info->g*barras[k].v + (1/lig->info->tap)*barras[m].v*( lig->info->g*cos(tKM + lig->info->phi)  +  lig->info->b*sin(tKM + lig->info->phi) );
		}

		lig = lig->prox;
	}


	return expr;
}

double calc_Nkm(int k,int m,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	tKM = barras[k].theta - barras[m].theta;
	ligacao * lig = ligacoes[k].prox;

	while(lig != NULL){
		if (lig->j == m)
		{
			if (lig->tBarra == TAP)
			{
				expr += (1/lig->info->tap)*barras[k].v*( lig->info->g*cos(tKM - lig->info->phi)  +  lig->info->b*sin(tKM - lig->info->phi) );
			}
			else
			{
				expr += (1/lig->info->tap)*barras[k].v*( lig->info->g*cos(tKM + lig->info->phi)  +  lig->info->b*sin(tKM + lig->info->phi) );
			}
		}
		lig = lig->prox;

	}


	return expr;
}

double calc_Mkk(int k,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	int m;
	ligacao * lig = ligacoes[k].prox;
	while(lig != NULL){
		m = lig->j ;

		tKM = barras[k].theta - barras[m].theta;

		if (lig->tBarra == TAP)
		{
			expr += -(1/lig->info->tap) * barras[k].v * barras[m].v * ( -lig->info->b * sin(tKM - lig->info->phi) - lig->info->g * cos(tKM - lig->info->phi));
		}
		else 
		{
			expr += -(1/lig->info->tap) * barras[k].v * barras[m].v * ( -lig->info->b * sin(tKM + lig->info->phi) - lig->info->g * cos(tKM + lig->info->phi));;
		}

		lig = lig->prox;
	}


	return expr;
}

double calc_Mkm(int k,int m,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	tKM = barras[k].theta - barras[m].theta;
	ligacao * lig = ligacoes[k].prox;

	while(lig != NULL){
		if (lig->j == m)
		{
			if (lig->tBarra == TAP)
			{
				expr += -(1/lig->info->tap) * barras[k].v * barras[m].v * ( -lig->info->b * sin(-tKM + lig->info->phi) + lig->info->g * cos(-tKM + lig->info->phi));
			}
			else 
			{
				expr += -(1/lig->info->tap) * barras[k].v * barras[m].v * ( lig->info->b * sin(tKM + lig->info->phi) + lig->info->g * cos(tKM + lig->info->phi));
			}
		}
		lig = lig->prox;

	}


	return expr;
}

double calc_Lkk(int k,barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	int m;
	ligacao * lig = ligacoes[k].prox;
	while(lig != NULL){
		m = lig->j ;

		tKM = barras[k].theta - barras[m].theta;

		if (lig->tBarra == TAP)
		{
			expr += 2*((1/pow(lig->info->tap,2)*lig->info->b ) + lig->info->bsh)* barras[k].v - (1/lig->info->tap)*barras[m].v*(lig->info->b*cos(tKM - lig->info->phi) - lig->info->g*sin(tKM - lig->info->phi));
		}
		else 
		{
			expr += 2*(lig->info->b)* barras[k].v - (1/lig->info->tap)*barras[m].v*(lig->info->b*cos(tKM + lig->info->phi) - lig->info->g*sin(tKM + lig->info->phi));
		}

		lig = lig->prox;
	}


	return expr;
}

double calc_Lkm(int k, int m, barra barras[],ligacao ligacoes []){
	double expr = 0.0;
	double tKM;
	tKM = barras[k].theta - barras[m].theta;
	ligacao * lig = ligacoes[k].prox;
	while(lig != NULL){
		if (lig->j == m)
		{

			if (lig->tBarra == TAP)
			{
				expr += -(1/lig->info->tap)*barras[k].v*(lig->info->b*cos(tKM - lig->info->phi) - lig->info->g*sin(tKM - lig->info->phi));
			}
			else 
			{
				expr += -(1/lig->info->tap)*barras[k].v*(lig->info->b*cos(tKM + lig->info->phi) - lig->info->g*sin(tKM + lig->info->phi));
			}

		}

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

		int k = li->m;

		// PARA CADA BARRA (lj->m => m) PQ OU PV
		lj = listaPQPV.prox;
		while(lj != NULL){
			
			int m = lj->m;

			jac[i][j] = 0.0;

			if (k == m){

				jac[i][j] = calc_Hkk(k,barras,ligacoes);

			} else {

				jac[i][j] = calc_Hkm(k,m,barras,ligacoes);
			}
			

			lj = lj->prox;
			j++;
		}

		// PARA CADA BARRA (lj->m => m) PQ
		lj = listaPQ.prox;
		while(lj != NULL){

			int m = lj->m;

			jac[i][j] = 0.0;

			
			if (k == m){

				jac[i][j] = (2*(barras[k].gsh*barras[k].v) + calc_Nkk(k,barras,ligacoes));

			} else {

				jac[i][j] = calc_Nkm(k,m,barras,ligacoes);
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

			
			if (k == m){

				jac[i][j] = calc_Mkk(k,barras,ligacoes);

			} else {

				jac[i][j] = calc_Mkm(k,m,barras,ligacoes);
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

				jac[i][j] = (2*(barras[k].bsh*barras[k].v) + calc_Lkk(k,barras,ligacoes));

			} else {

				jac[i][j] = calc_Lkm(k,m,barras,ligacoes);
			}

			lj = lj->prox;
			j++;
		}

		li = li->prox;
		i++;
	}

}