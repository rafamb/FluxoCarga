
#define PI 3.14159265

#define MOD_GUILHERME 0
#define MOD 1

#define TAP 1
#define IMPEDANCIA 0


int modo;

FILE * carregarArquivo(char arquivo []){
	FILE *arq;
	arq = fopen(arquivo, "r");
	if (arq == NULL){
		printf("Erro, nao foi possivel abrir o arquivo\n");
		return NULL;
	} else {
		char aux[6];
		int i = 0;

		do {
			aux[i] = fgetc(arq);
			i++;
		}while(i < 5);

		if (strcmp(aux,"TAPE\n"))
		{
			modo = MOD_GUILHERME;
		}
		else
		{
			modo = MOD;
			fseek(arq,0,SEEK_SET);
		}

		return arq;
	}
}

int carregarBaseMVA(FILE * arq){
	fseek(arq,30,SEEK_CUR);
	double bMVA;
	fscanf(arq,"%lf",&bMVA);
	fseek(arq,59,SEEK_CUR);
	return bMVA;
}

int carregarnB(FILE * arq){
	int nB,aux;
	ssize_t volta = 0;
	char * c = NULL;
	size_t len = 0;
	fscanf(arq,"%*d ITEMS\n");

	do{

		fscanf(arq,"%d",&aux);
		if (aux != -999)
		{
			nB = aux;
		}
		else
		{
			break;
		}

		volta += getline(&c,&len,arq) + sizeof(int);

	}while(1);

	fseek(arq,-(volta+4),SEEK_CUR);
	/*nB = volta/(130+sizeof(int)) -1;*/

	/*if (modo == MOD)
	{
		do{
			fscanf(arq,"%d",&nB);
			
			
			fseek(arq,125,SEEK_CUR);
			volta = volta + 125 + sizeof(nB);
		}while(nB != -999);
		
		
		fseek(arq,-volta+2,SEEK_CUR);
		nB = volta/(125+sizeof(int)) -1;
	}
	else 
	{
		do{
			fscanf(arq,"%d",&nB);
			
			fseek(arq,130,SEEK_CUR);
			volta = volta + 130 + sizeof(nB);
		}while(nB != -999);
		
		
		fseek(arq,-volta+2,SEEK_CUR);
		nB = volta/(130+sizeof(int)) -1;
	}*/
	/*getline(&c,&len,arq);
	printf("%s\n", c);*/


	if (c)
	{
		free(c);
	}

	
	
	return nB;
}

void carregarBarras(FILE *arq, barra barras [], double baseMVA, int * nPQ, int *nPV, int *ref, lista * listaPQPV, lista * listaPQ){
	int i;
	int count = 0;
	fscanf(arq,"%d",&i);
	while(i != -999){
		fgetc(arq);
		do {
			barras[i-1].nome[count] = fgetc(arq);
			count++;
		}while(count < 12);
		fscanf(arq, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf",&barras[i-1].area,&barras[i-1].zonaPerdas,&barras[i-1].tipo,&barras[i-1].v,&barras[i-1].theta,&barras[i-1].pc,&barras[i-1].qc,&barras[i-1].pg,&barras[i-1].qg,&barras[i-1].baseKV,&barras[i-1].vgO);
		barras[i-1].vEsp = barras[i-1].v;
		barras[i-1].theta = barras[i-1].theta*PI/180.0;
		barras[i-1].pc = barras[i-1].pc/baseMVA;
		barras[i-1].qc = barras[i-1].qc/baseMVA;
		barras[i-1].pg = barras[i-1].pg/baseMVA;
		barras[i-1].qg = barras[i-1].qg/baseMVA;


		if (barras[i-1].tipo == 3 || barras[i-1].tipo == 2)
		{
			if (barras[i-1].tipo == 2)
			{
				*nPV = *nPV+1;
				inserirLista(listaPQPV,i-1);
			} else {
				*ref = i-1;
			}

			fscanf(arq, "%lf %lf ",&barras[i-1].qgMax,&barras[i-1].qgMin);
			if (barras[i-1].qgMin == 0 && barras[i-1].qgMax == 0)
			{
				barras[i-1].qgMin = 999.999999;
				barras[i-1].qgMax = -999.999999;
			} else {
				barras[i-1].qgMin = barras[i-1].qgMin/baseMVA;
				barras[i-1].qgMax = barras[i-1].qgMax/baseMVA;
			}

			barras[i-1].vMin = 0.9;
			barras[i-1].vMax = 1.1;


			
		}else{ 
			*nPQ = *nPQ +1;
			inserirLista(listaPQPV,i-1);
			inserirLista(listaPQ,i-1);
			if (barras[i-1].tipo == 1){
					barras[i-1].qgMin = 0.0;
					barras[i-1].qgMax = 0.0;
		
					fscanf(arq, "%lf %lf ",&barras[i-1].vMax,&barras[i-1].vMin);
			} else {
		
					barras[i-1].qgMin = 0.0;
					barras[i-1].qgMax = 0.0;
		
					barras[i-1].vMin = 0.9;
					barras[i-1].vMax = 1.1;
					fscanf(arq, "%*f %*f");
		
			}
		}

		fscanf(arq, "%lf %lf %lf",&barras[i-1].gsh,&barras[i-1].bsh,&barras[i-1].ctrlREM);
		

		if (modo == MOD_GUILHERME)
		{
			fscanf(arq, "%*d");
		}
		
		count = 0;
		fscanf(arq,"%d",&i);
	}
	char c;
	do {
		c = fgetc(arq);
	}while(c != '\n');
}

void solucaoInicial(barra barras [], int nB, int ref){
	int i;
	for(i=0; i<nB; i++){
		if (barras[i].tipo == 0 || barras[i].tipo == 1)
		{
			barras[i].v = 1.0;
			barras[i].theta = barras[ref].theta;
		} else if (barras[i].tipo == 2)
		{
			barras[i].theta = barras[ref].theta;
		}
	}
}


void inicializarLigacoes(barra barras[],ligacao ligacoes [], int nB){
	for (int i = 0; i < nB; ++i)
	{
		ligacoes[i].j = i;
		ligacoes[i].info = NULL;
		ligacoes[i].prox = NULL;
	}
}

void carregarLigacoes(FILE *arq,barra barras[], ligacao ligacoes []){
	fscanf(arq,"BRANCH DATA FOLLOWS %*d ITEMS\n");
	int i,j;
	double gkm, bkm;

	fscanf(arq,"%d",&i);
	while(i != -999){
		infoLigacao * novaLigacao = (infoLigacao *)malloc(sizeof(infoLigacao));
		novaLigacao->i = i;
		fscanf(arq,"%d",&j);
		novaLigacao->j = j;
		fscanf(arq,"%d",&novaLigacao->area);
		fscanf(arq,"%d",&novaLigacao->zonaPerdas);
		fscanf(arq,"%d",&novaLigacao->circParalelos);
		fscanf(arq,"%d",&novaLigacao->tipo);
		fscanf(arq,"%lf",&novaLigacao->r);
		fscanf(arq,"%lf",&novaLigacao->x);
		fscanf(arq,"%lf",&novaLigacao->bsh);
		novaLigacao->bsh = novaLigacao->bsh/2.0;
		fscanf(arq,"%d",&novaLigacao->lineRat1);
		fscanf(arq,"%d",&novaLigacao->lineRat2);
		fscanf(arq,"%d",&novaLigacao->lineRat3);
		fscanf(arq,"%d",&novaLigacao->tControl);
		fscanf(arq,"%d",&novaLigacao->side);
		fscanf(arq,"%lf",&novaLigacao->tap);
		fscanf(arq,"%lf",&novaLigacao->phi);
		novaLigacao->phi = novaLigacao->phi*PI/180.0;

		if (novaLigacao->tipo == 0 || novaLigacao->tap == 0)
		{
			novaLigacao->tap = 1.0;
		}

		/*novaLigacao->tap = 1/novaLigacao->tap;*/

		if (novaLigacao->tipo == 2 || novaLigacao->tipo == 3){
			fscanf(arq,"%lf %lf %lf %lf %lf",&novaLigacao->tapMin,&novaLigacao->tapMax,&novaLigacao->passo,&novaLigacao->ctrlMin,&novaLigacao->ctrlMax);
		} else {
			fscanf(arq,"%*f %*f %*f %*f %*f");
			novaLigacao->tapMin = 0.0;
			novaLigacao->tapMax = 0.0;
			novaLigacao->passo = 0.0;
			novaLigacao->ctrlMax = 0.0;
			novaLigacao->ctrlMin = 0.0;
		}

		ligacao * ligI = (ligacao *)malloc(sizeof(ligacao));
		ligacao * ligJ = (ligacao *)malloc(sizeof(ligacao));

		ligI->info = novaLigacao;
		ligI->j = j-1;
		ligJ->j = i-1;
		ligJ->info = novaLigacao;
		novaLigacao->pI = ligI;
		novaLigacao->pJ = ligJ;

		ligI->prox = ligacoes[i-1].prox;
		ligJ->prox = ligacoes[j-1].prox;
		ligacoes[i-1].prox = ligI;
		ligacoes[j-1].prox = ligJ;

		gkm = novaLigacao->r/(pow(novaLigacao->r,2)+pow(novaLigacao->x,2));
		bkm = -novaLigacao->x/(pow(novaLigacao->r,2)+pow(novaLigacao->x,2));
		novaLigacao->g = gkm;
		novaLigacao->b = bkm;

		/*ligI->g = - sin(novaLigacao->phi)*(novaLigacao->tap)*bkm - cos(novaLigacao->phi)*(novaLigacao->tap)*gkm;
		ligI->b = sin(novaLigacao->phi)*(novaLigacao->tap)*gkm - cos(novaLigacao->phi)*(novaLigacao->tap)*bkm;

		ligJ->g = - sin(-novaLigacao->phi)*(novaLigacao->tap)*bkm - cos(-novaLigacao->phi)*(novaLigacao->tap)*gkm;
		ligJ->b = sin(-novaLigacao->phi)*(novaLigacao->tap)*gkm - cos(-novaLigacao->phi)*(novaLigacao->tap)*bkm;

		ligacoes[i-1].g = ligacoes[i-1].g + pow(novaLigacao->tap,2)*gkm;
		ligacoes[i-1].b = ligacoes[i-1].b + pow(novaLigacao->tap,2)*bkm + novaLigacao->bsh;

		ligacoes[j-1].g = ligacoes[j-1].g + pow(novaLigacao->tap,2)*gkm;
		ligacoes[j-1].b = ligacoes[j-1].b + pow(novaLigacao->tap,2)*bkm + novaLigacao->bsh;*/

		if (novaLigacao->tipo == 0)
		{
			ligI->tBarra = TAP;
			ligJ->tBarra = TAP;
		}
		else
		{
			ligI->tBarra = TAP;
			ligJ->tBarra = IMPEDANCIA;
		}


		if (modo == MOD_GUILHERME)
		{
			fscanf(arq, "%*d");
		}


		fscanf(arq,"%d",&i);
	}

}

void liberarMemoriaLigacoes(ligacao ligacoes [], int nB){
	int i;
	ligacao * atual;
	infoLigacao * info;
	for(i = 0; i < nB; i++){
		atual = ligacoes[i].prox;
		while(atual != NULL){
			if (atual->info != NULL)
			{
				info = atual->info;
				info->pI->info = NULL;
				info->pJ->info = NULL;
				free(info);
			}
			ligacoes[i].prox = atual->prox;
			free(atual);
			atual = ligacoes[i].prox;
		}
	}
}

void printSolucao(int nB, barra barras []){
	int i;

	printf("Barra\tV[p.u.]\t\tAng[rad]\tPg[p.u.]\tQg[p.u.]\tPc[p.u.]\tQc[p.u.]\n");

	for(i = 0; i < nB; i++)
	{
		printf("%d:\t", i + 1);
		printf("%lf\t", barras[i].v);
		printf("%lf\t", barras[i].theta);
		printf("%lf\t", barras[i].pg);
		printf("%lf\t", barras[i].qg);
		printf("%lf\t", barras[i].pc);
		printf("%lf\n", barras[i].qc);
	}

}