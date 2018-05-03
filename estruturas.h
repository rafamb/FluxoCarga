//
// AQUI SÃO DEFINIDAS AS ESTRUTURAS BARRA, INFOLIGACAO, LIGACAO, E LISTA
//


//
// BARRA: GUARDA AS INFORMACOES DE UMA BARRA DO SISTEMA
//
typedef struct{
	char nome[20];
	int area;
	int zonaPerdas;
	int tipo;
	double v;
	double theta;
	double pc;
	double qc;
	double pg;
	double qg;
	double baseKV;
	double vgO;
	double qgMax;
	double qgMin;
	double vMax;
	double vMin;
	double gsh;
	double bsh;
	double ctrlREM;
	double p;
	double q;
}barra;


//
// INFOLIGACAO: GUARDA AS INFORMACOES DE UMA LIGACAO DO SISTEMA
//
typedef struct{
	int i;
	int j;
	int area;
	int zonaPerdas;
	int circParalelos;
	int tipo;
	double r;
	double x;
	double bsh;
	int lineRat1;
	int lineRat2;
	int lineRat3;
	int tControl;
	int side;
	double tap;
	double phi;
	double tapMin;
	double tapMax;
	double passo;
	double ctrlMin;
	double ctrlMax;
	double g;
	double b;
	struct lig *pI;
	struct lig *pJ;
}infoLigacao;

//
// LIGACAO: POSSUI UM PONTEIRO PARA AS INFORMACOES DA LIGACAO A QUE SE REFERE,
// ALEM DOS CAMPOS G E B, DA INFORMACAO DE PARA QUAL NO ESSA LIGACAO SE ENCAMINHA (J),
// E UM PONTEIRO PARA A PROXIMA LIGACAO DO NO
//
typedef struct lig{
	infoLigacao * info;
	int j;
	double g;
	double b;
	struct lig *prox;
}ligacao;

//
// LISTA: POSSUI UM PONTEIRO PARA A PROXIMA LIGACAO E PARA A LIGACAO ANTERIOR, ALEM DE GUARDAR
// UM VALOR NUMERICO M
//
typedef struct list{
	int m;
	struct list *prox;
	struct list *ant;
}lista;


//
// INSERIRLISTA: FUNCAO RESPONSAVEL POR ADICIONAR UM ELEMENTO NA LISTA, RECEBE UM NUMERO (NUM),
// E O INSERE AO FIM DA LISTA PL
//
void inserirLista(lista *pL, int num){
	lista * elem = (lista *)malloc(sizeof(lista));
	elem->m = num;
	if (pL->prox == NULL)
	{
		elem->prox = pL->prox;
		elem->ant = pL->ant;
		pL->prox = elem;
		pL->ant = elem;
	} else {
		elem->prox = NULL;
		elem->ant = pL->ant;
		pL->ant->prox = elem;
		pL->ant = elem;
	}
	
}


//
// LIBERARLISTA: FUNCAO RESPONSAVEL POR DESALOCAR A MEMORIA UTILIZADA PELA LISTA PL
//
void liberarLista(lista * pL){
	lista * atual = pL->prox;
	
	while(atual != NULL){
		pL->prox = atual->prox;
		free(atual);
		atual = pL->prox;
	}
	pL->ant = NULL;

}