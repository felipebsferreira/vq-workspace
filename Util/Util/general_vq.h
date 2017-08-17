//
//  general_vq.h
//  Util
//
//  Created by Felipe Ferreira on 29/01/16.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef general_vq_h
#define general_vq_h

#include <cstdlib>
#include <ctime>
#include <cmath>
#include "inout.h"
#include "definicoes.h"
#include "general.h"

using namespace std;

bool validar_dicionarios(const char *imagem, int qtd_dic, int qtd_sim);
void gerar_multiplos_dicionarios(const char *imagem, int N_inicial, int N_final, int qtd_dic, int qtd_sim);
void gerar_dicionarios(const char *nome_dic, const char *nome_treino, int N, int K, int qtd_dic);
double **gerar_dicionario(double **treino, int N, int K, int n);

double calcular_distancia_euclidiana(double *x, double *y, int K);

double CalcularPSNR(char orig[50], char quant[50]);
double CalcularPSNR(double **orig, double **quant, int nTreino, int K, int dim);

double **montar_matriz_probabilidades(int N, double epsilon);

const char *converter_nns_string(NNS nns);
NNS converter_str_nns(const char * nns);
const char *converter_swarm_string(SWARM_TYPE type);
SWARM_TYPE converter_str_swarm(const char * type);

void reconstruir_imagem(const char *nome, double **dic, int *codigo, int K, int nTreino, int bloco, int dim_img);

bool validar_dicionarios(const char *imagem, int qtd_dic, int qtd_sim)
{
    char nome[50];
    double **dic;
    bool igual = true;

    for (int n = 32; n <= 256; n = n*2)
    {
        for (int c = 0; c < qtd_sim; c++)
        {
            sprintf(nome, "dics/%i_%s_%i.dic", n, imagem, qtd_sim);
            dic = ler_arquivo(nome, n, qtd_dic, 16);

            for (int h = 1; h < qtd_dic; h++)
            {
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < 16; j++)
                    {
                        if(dic[i][j] != dic[(h * n) + i][j])
                        {
                            igual = false;
                            break;
                        }
                    }
                }
            }
        }
    }

    return igual;
}

void gerar_multiplos_dicionarios(const char *imagem, int N_inicial, int N_final, int qtd_dic, int qtd_sim)
{
    char dic[50], treino[50];

    for (int n = N_inicial; n <= N_final; n = n*2)
    {
        cout << "N = " << n << endl;

        for (int i = 0; i < qtd_sim; i++)
        {
            cout << "Simulacao = " << i << endl;
            sprintf(dic, "dics/%i_%s_%i.dic", n, imagem, i + 1);
            sprintf(treino, "imagens/%s.dat", imagem);
            gerar_dicionarios(dic, treino, n, 16, qtd_dic);
        }
    }
}

void gerar_dicionarios(const char *nome_dic, const char *nome_treino, int N, int K, int qtd_dic)
{
    double **treino, **dic;
    int nTreino;

    for(int i = 0; i < qtd_dic; i++)
    {
        treino = ler_arquivo_treino(nome_treino, &nTreino, K);
        dic = gerar_dicionario(treino, N, K, nTreino);
        escreve_arquivo(nome_dic, dic, N, K, "a");
    }
}

double **gerar_dicionario(double **treino, int N, int K, int n)
{
    double **dic;
    int indice;

    dic = aloca_matrizd(N, K);

    for(int i = 0; i < N; i++)
    {
        indice = rand()%n;

        for(int j = 0; j < K; j++)
            dic[i][j] = treino[indice][j];
    }

    return dic;
}

// Calcula a distancia euclidiana entre os vetores (dimensao K) pta e ptb
double calcular_distancia_euclidiana(double *x, double *y, int K)
{
    int i;
    double dist;

    dist = 0.0;

    for(i = 0; i < K; i++)
        dist += (x[i] - y[i]) * (x[i] - y[i]);

    dist = sqrt(dist);

    return dist;
}

double CalcularPSNR(char orig[50], char quant[50])
{
    FILE *origp, *quantp;
    int i;
    double x[SEG], y[SEG], snrtotal, varxyseg;

    /* abre arquivos para leitura e escrita de resultados */

    origp = fopen(orig, "r");
    quantp = fopen(quant, "r");

    while(!feof(origp))
    {
        for(i=0; i < SEG; i++)
            fscanf(origp, "%lf\n", &x[i]);

        for(i=0; i < SEG; i++)
            fscanf(quantp, "%lf\n", &y[i]);

        for(i=0; i < SEG; i++)
            varxyseg += (x[i]-y[i]) * (x[i]-y[i]);
    }

    snrtotal = 10*(log10(pow(255,2.0)) - log10(varxyseg/65536));

    fclose(origp);
    fclose(quantp);

    return snrtotal;
}

double CalcularPSNR(double **orig, double **quant, int nTreino, int K, int dim)
{
    int i, j, c = 0;
    double x[SEG], y[SEG], snrtotal, varxyseg;
    double *orig_seg, *quant_seg;

    varxyseg = snrtotal = 0;

    orig_seg = aloca_vetor_double(dim*dim);
    quant_seg = aloca_vetor_double(dim*dim);

    for(i = 0; i < nTreino; i++)
    {
        for(j = 0; j < K; j++)
        {
            orig_seg[c] = orig[i][j];
            quant_seg[c] = quant[i][j];
            c++;
        }
    }

    for(j = 0; j < dim*dim; j += SEG)
    {
        for(i = 0; i < SEG; i++)
            x[i] = orig_seg[j + i];

        for(i=0; i < SEG; i++)
            y[i] = quant_seg[j + i];

        for(i=0; i < SEG; i++)
            varxyseg += (x[i]-y[i]) * (x[i]-y[i]);
    }

    snrtotal = 10*(log10(pow(255,2.0)) - log10(varxyseg/65536));

    free(orig_seg);
    free(quant_seg);

    return snrtotal;
}


double **montar_matriz_probabilidades(int N, double epsilon)
{
    /******************************************************************
     * Esta função monta uma matriz "p" para representar a             *
     * a probabilidade de cruzamento dos índices devido                *
     * ao erro de canal P(j|b(i))                                      *
     ******************************************************************/

    double comp;
    double **matriz_prob;
    int i, j, d, var_1, var_2, qtd_bits, desloc_bit, op_xor;

    matriz_prob = aloca_matrizd(N,N);

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            matriz_prob[i][j] = 1;
        }
    }

    /* Cálculo do número de bits dos índices do dicionário */
    qtd_bits = log2(N);

    var_1 = 0;
    var_2 = 0;

    for(i=var_1;i<N;i++){
        for(j=var_2;j<N;j++){
            if(i==j){
                matriz_prob[i][j] = pow((1 - epsilon), qtd_bits);
            }

            else{
                op_xor = i ^ j;
                desloc_bit = 1;
                for(d=0;d<qtd_bits;d++){
                    comp = op_xor & desloc_bit;
                    if(comp==0){
                        matriz_prob[i][j] = matriz_prob[i][j] * (1 - epsilon);
                        matriz_prob[j][i] = matriz_prob[i][j];
                    }
                    else {
                        matriz_prob[i][j] = matriz_prob[i][j] * epsilon;
                        matriz_prob[j][i] = matriz_prob[i][j];
                    }
                    desloc_bit = desloc_bit << 1;
                }
            }
        }
        var_1++;
        var_2++;
    }

    return matriz_prob;
}

const char *converter_nns_string(NNS nns)
{
    switch (nns)
    {
        case BT:
            return "BT";

        case PDS:
            return "PDS";

        case ENNS:
            return "ENNS";

        case IENNS:
            return "IENNS";

        case EENNS:
            return "EENNS";

        case IEENNS:
            return "IEENNS";

        case EEENNS:
            return "EEENNS";

        case DTA:
            return "DTA";

        case IDTA:
            return "IDTA";

        default:
            return "";
    }
}

NNS converter_str_nns(const char * nns)
{
    if (str_isequal(nns, "PDS")) return PDS;
    if (str_isequal(nns, "ENNS")) return ENNS;
    if (str_isequal(nns, "IENNS")) return IENNS;
    if (str_isequal(nns, "EENNS")) return EENNS;
    if (str_isequal(nns, "IEENNS")) return IEENNS;
    if (str_isequal(nns, "EEENNS")) return EEENNS;
    if (str_isequal(nns, "DTA")) return DTA;
    if (str_isequal(nns, "IDTA")) return IDTA;

    return BT;
}

const char *converter_swarm_string(SWARM_TYPE type)
{
    switch (type)
    {
        case NO_SWARM:
            return "NO_MOEA";

        case FSS:
            return "FSS";

        case PSO:
            return "PSO";

        case FF:
            return "FF";

        case HBMO:
            return "HBMO";

        case FUZZY:
            return "FUZZY";

        case FF_NH:
            return "FF_NH";

        default:
            return "";
    }
}

SWARM_TYPE converter_str_swarm(const char * type)
{
    if (str_isequal(type, "NO_MOEA")) return NO_SWARM;
    else if (str_isequal(type, "FSS")) return FSS;
    else if (str_isequal(type, "PSO")) return PSO;
    else if (str_isequal(type, "FF")) return FF;
    else if (str_isequal(type, "FUZZY")) return FUZZY;
    else if (str_isequal(type, "FF_NH")) return FF_NH;
    else return HBMO;
}

const char *converter_scaling_string(SCALING type)
{
    switch (type)
    {
        case NONE:
            return "NONE";

        case FIXED:
            return "FIXED";

        case VARIABLE:
            return "VAR";

        default:
            return "";
    }
}

SCALING converter_str_scaling(const char * type)
{
    if (str_isequal(type, "NONE")) return NONE;
    else if (str_isequal(type, "FIXED")) return FIXED;
    else return VARIABLE;
}

void reconstruir_imagem(const char *nome, double **dic, int *codigo, int K, int nTreino, int bloco, int dim_img)
{
    FILE *arq;
    int tam, k;
    int x, y, x_aux, y_aux;
    int **matrix;

    tam = dim_img/bloco;
    matrix = aloca_matriz_int(dim_img, dim_img);

    arq = fopen(nome, "w");

    fprintf(arq, "P2 256 256 255 ");

    x_aux = y_aux = 0;

    for (int i = 0; i < nTreino; i++)
    {
        k = 0;

        for (y = 0; y < bloco; y++)
        {
            for (x = 0; x < bloco; x++)
            {
                matrix[y + y_aux * bloco][x + x_aux * bloco] = dic[codigo[i]][k];
                k++;
            }
        }

        x_aux++;

        if (x_aux == 64)
        {
            x_aux = 0;
            y_aux++;
        }
    }

    for (y = 0; y < dim_img; y++)
    {
        for (x = 0; x < dim_img; x++)
        {
            fprintf(arq, "%i ", matrix[y][x]);
        }
    }

    desaloca_matriz_int(matrix, dim_img);

    fclose(arq);
}

#endif /* general_vq_h */
