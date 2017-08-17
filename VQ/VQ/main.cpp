//
//  main.cpp
//  VQ
//
//  Created by Felipe Ferreira on 08/12/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#include <fstream>
#include "vq.hpp"
#include "fss_vq.hpp"
#include "ff_vq.hpp"
#include "pso_vq.hpp"
#include "hbmo_vq.hpp"
#include "swarm_vq.hpp"
#include "fuzzy_k_means.hpp"
#include "../../Util/Util/general.h"
#include "../../Util/Util/general_vq.h"
#include "../../Util/Util/codec.h"
#include "../../Util/Util/bsc.h"

void RodarSimulacao(const char *imagem, int N_init, int N_fim, NNS nns, SWARM_TYPE swarm_type, SCALING scaleType, double scale, double m);
SWARM_VQ *CriarInstancia(SWARM_TYPE swarm_type, int N, int K, int qtd_part, NNS nns, SCALING scaleType, double scale, double m);
double *RodarSwarmVq(SWARM_VQ *swarm_vq, const char *imagem, int simulacao);
void TestarParametros(SWARM_TYPE swarm_type, int N);
void ExecutarEmLote();

char _str_swarm[50];
int _N_init, _N_fim;

int main()
{
    reset_rand();

    sprintf(_str_swarm, "fss");
    RodarSwarmVq(CriarInstancia(NO_SWARM, 32, 16, 10, IEENNS, NONE, 9, 2), "lena", 1);

    return 0;
}

void ExecutarEmLote()
{
    char imagem[50], str_nns[50], str_scaling[50];
    int qtd_nns = 0, qtd_imagens = 0;
    NNS *nns_list = new NNS[8];
    char imagens[10][50];
    SWARM_TYPE swarm_type;
    SCALING scaleType;
    float scale;
    FILE *entrada = fopen("entrada.txt", "rt");

    while(true)
    {
        fscanf(entrada, "%s ", imagem);

        if (!str_isequal(imagem, "eoi"))
        {
            sprintf(imagens[qtd_imagens], "%s", imagem);
            qtd_imagens++;
        }
        else
            break;
    }

    fscanf(entrada, "%i %i %s %s %f ", &_N_init, &_N_fim, _str_swarm, str_scaling, &scale);

    swarm_type = converter_str_swarm(_str_swarm);
    scaleType = converter_str_scaling(str_scaling);

    while(true)
    {
        fscanf(entrada, "%s ", str_nns);

        if (!str_isequal(str_nns, "eof"))
        {
            nns_list[qtd_nns] = converter_str_nns(str_nns);
            qtd_nns++;
        }
        else
            break;
    }

    for(int i = 0; i < qtd_imagens; i++)
    {
        for(int j = 0; j < qtd_nns; j++)
        {
            RodarSimulacao(imagens[i], _N_init, _N_fim, nns_list[j], swarm_type, scaleType, scale, 0);
        }
    }
}

void TestarParametros(SWARM_TYPE swarm_type, int N)
{
    SWARM_VQ *instance;
    double *resultados, psnr;
    ofstream myfile;

    myfile.open("dados.csv");
    myfile << "V1" << SEPARADOR << "V2" << SEPARADOR << "PSNR\n";

    double add = 0.1, add2 = 0.2, add_aux;

    for (double i = 0.1; i < 2.1; i = i + 0.3)
    {
        add_aux = add;
        add = add2;
        add2 = add;

        for (double j = add; j < 2.1; j = j + 0.3)
        {
            psnr = 0;

            for (int n = 0; n < 10; n++)
            {
                cout << i << " / " << j << " - " << n + 1 << endl;

                instance = CriarInstancia(swarm_type, N, 16, 10, IEENNS, NONE, 9, 0);

                switch (swarm_type)
                {
                    case FSS:
                        ((FSS_VQ *)instance)->stepInd = i;
                        ((FSS_VQ *)instance)->stepVol = j;
                        break;

                    default:
                        exit(1);
                }

                resultados = RodarSwarmVq(instance, "clock", n + 1);

                psnr += resultados[3];
                free(resultados);
            }

            myfile << FormatNumber(i) << SEPARADOR << FormatNumber(j) << SEPARADOR  << FormatNumber(psnr/10.0) << "\n";
        }
    }

    myfile.close();
}

void RodarSimulacao(const char *imagem, int N_init, int N_fim, NNS nns, SWARM_TYPE swarm_type, SCALING scaleType, double scale, double m)
{
    double *resultados;
    double iteracoes, duracao, val_psnr, dist;
    int N, N_limite, qtd_part = 10, qtd_simulacoes = 30;
    char nome[50];
    SWARM_VQ *swarm_vq;
    ofstream myfile;

    if (N_init != N_fim)
        N_limite = 0;
    else
        N_limite = N_init;

    sprintf(nome, "resultados/r_%s_%i_%i_%s_%s_%s.csv", imagem, N_init, N_fim, converter_nns_string(nns), converter_swarm_string(swarm_type), converter_scaling_string(scaleType));
    myfile.open(nome);

    myfile << "N;Algoritmo\n";
    myfile << "N;Algoritmo;PSNR;Iter;Tempo;Distorcao\n";

    N = N_init;
    N_limite = N_fim;

    for(; N <= N_limite; N = N*2)
    {
        myfile << N << SEPARADOR << converter_nns_string(nns);

        iteracoes = duracao = val_psnr = dist = 0;

        for(int i = 0; i < qtd_simulacoes; i++)
        {
            cout << "N = " << N << " / " << imagem << " / " << i + 1 << endl << endl;

            swarm_vq = CriarInstancia(swarm_type, N, 16, qtd_part, nns, scaleType, scale, m);
            resultados = RodarSwarmVq(swarm_vq, imagem, i + 1);

            iteracoes += resultados[0];
            duracao += resultados[1];
            dist += resultados[2];
            val_psnr += resultados[3];
        }

        iteracoes = iteracoes/qtd_simulacoes;
        duracao = duracao/qtd_simulacoes;
        val_psnr = val_psnr/qtd_simulacoes;
        dist = dist/qtd_simulacoes;

        myfile << SEPARADOR << val_psnr << SEPARADOR << iteracoes;
        myfile << SEPARADOR << duracao << SEPARADOR << dist;

        myfile << "\n";
    }

    myfile.close();
}

SWARM_VQ *CriarInstancia(SWARM_TYPE swarm_type, int N, int K, int qtd_part, NNS nns, SCALING scaleType, double scale, double m)
{
    SWARM_VQ *swarmVq;

    switch (swarm_type)
    {
        case NO_SWARM:
            swarmVq = new VQ(N, K, nns, scaleType, scale);
            break;

        case FSS:
            swarmVq = new FSS_VQ(N, K, qtd_part, nns);
            break;

        case PSO:
            swarmVq = new PSO_VQ(N, K, qtd_part, nns, scaleType, scale);
            break;

        case FF:
            swarmVq = new FF_VQ(N, K, qtd_part, nns);
            break;

        case HBMO:
            cout << endl << "Tecnica SWARM não implementada!" << endl;
            exit(1);

        case FUZZY:
            swarmVq = new FuzzyKMeans(N, K, nns);
            break;

        default:
            cout << endl << "Tecnica SWARM não selecionada corretamente!" << endl;
            exit(1);
    }

    return swarmVq;
}

double *RodarSwarmVq(SWARM_VQ *swarm, const char *imagem, int simulacao)
{
    double **dic, **treino, **imagem_quantizada;
    double *resultados;
    double psnr = 0;
    char nome[50], nome_recon[50];
    int *p_codificada;
    int nTreino;

    sprintf(nome, "imagens/%s.dat", imagem);
    treino = ler_arquivo_treino(nome, &nTreino, swarm->K);

    sprintf(nome, "dics/%i_%s_%i.dic", swarm->N, imagem, simulacao);
    dic = ler_arquivo(nome, swarm->N * swarm->qtd_part, swarm->K);

    // Executar quantização
    swarm->nTreino = nTreino;
    resultados = swarm->Run(dic, treino);

    // Codificar
    p_codificada = codificar(dic, treino, swarm->N, swarm->K, nTreino);
    imagem_quantizada = decodificar(dic, p_codificada, swarm->K, swarm->N, nTreino);
    psnr = CalcularPSNR(treino, imagem_quantizada, nTreino, swarm->K, 256);

    resultados[3] = psnr;

    cout << endl << "Iteracoes: " << resultados[0] << endl;
    cout << "Duracao: " << resultados[1] << "s" << endl;
    cout << "Distorcao: " << resultados[2] << endl;
    cout << "PSNR = " << resultados[3] << endl << endl;

    sprintf(nome_recon, "%i_%s_%s_%f.pgm", swarm->N, imagem, _str_swarm, psnr);
    reconstruir_imagem(nome_recon, dic, p_codificada, swarm->K, swarm->nTreino, 4, 256);

    // Desalocar memória
    desaloca_matrizd(dic, swarm->N * swarm->qtd_part);
    desaloca_matrizd(treino, nTreino);
    desaloca_matrizd(imagem_quantizada, nTreino);
    free(p_codificada);
    delete swarm;

    return resultados;
}
