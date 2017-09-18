//
//  main.cpp
//  COVQ
//
//  Created by Felipe Ferreira on 08/12/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#include <fstream>
#include "covq.hpp"
#include "pso_covq.hpp"
#include "fss_covq.hpp"
#include "ff_covq.hpp"
#include "fw_covq.hpp"
#include "swarm_covq.hpp"
#include "../../Util/Util/general.h"
#include "../../Util/Util/codec.h"
#include "../../Util/Util/bsc.h"

void RodarSimulacao(const char *imagem, int N_limite, NNS nns, SWARM_TYPE swarm_type, int qtd_simulacoes, SCALING scaleType, double scale);
SWARM_COVQ *CriarInstancia(SWARM_TYPE swarm_type, int N, int K, int qtd_part, double erro, SCALING scaleType, double scale);
double *RodarMoeaCovq(SWARM_COVQ *swarm, const char *imagem, NNS nns, int simulacao);

int main()
{
    reset_rand();

    const char *str_image = "lena";
    int N = 32;
    NNS nns = IEENNS;
    SWARM_TYPE swarm_type = FW;
    int qtd_simulacoes = 10;
    SCALING scale_type = NONE;
    double scale = 0;
    
    /*Rodar conjunto de simulações*/
    /******************************/
    RodarSimulacao(str_image, N, nns, swarm_type, qtd_simulacoes, scale_type, scale);

    /******************************/
    /*Rodar apenas uma simulação  */
    /******************************/
    //RodarMoeaCovq(CriarInstancia(swarm_type, N, 16, 10, 0.005, scale_type, scale), str_image, nns, 1);

    return 0;
}

void RodarSimulacao(const char *imagem, int N_limite, NNS nns, SWARM_TYPE swarm_type, int qtd_simulacoes, SCALING scaleType, double scale)
{
    double *resultados;
    double iteracoes, duracao, val_psnr, dist;
    double erro[4] = {0.005, 0.01, 0.05, 0.1};
    int N, qtd_part = 10;
    char nome[50];
    SWARM_COVQ *swarm_covq;
    ofstream myfile;

    sprintf(nome, "resultados/r_%s_%i_%s_%s_%s.csv", imagem, N_limite, converter_nns_string(nns), converter_swarm_string(swarm_type), converter_scaling_string(scaleType));
    myfile.open(nome);

    myfile << "N;Algoritmo;0.005;0.005;0.005;0.005;0.010;0.010;0.010;0.010;0.050;0.050;0.050;0.050;0.100;0.100;0.100;0.100\n";
    myfile << "N;Algoritmo;PSNR;Iter;Tempo;Distorcao;PSNR;Iter;Tempo;Distorcao;PSNR;Iter;Tempo;Distorcao;PSNR;Iter;Tempo;Distorcao\n";

    if(N_limite == 0)
    {
        N = 32;
        N_limite = 256;
    }
    else
    {
        N = N_limite;
    }

    for(; N <= N_limite; N = N*2)
    {
        myfile << N << SEPARADOR << converter_nns_string(nns);

        for(int j = 0; j < 4; j++)
        {
            iteracoes = duracao = val_psnr = dist = 0;

            for(int i = 0; i < qtd_simulacoes; i++)
            {
                cout << "N = " << N << " / erro = "<< erro[j] << " / " << imagem << " / " << i + 1 << endl << endl;

                swarm_covq = CriarInstancia(swarm_type, N, 16, qtd_part, erro[j], scaleType, scale);
                resultados = RodarMoeaCovq(swarm_covq, imagem, nns, i + 1);

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
        }

        myfile << "\n";
    }

    myfile.close();
}

SWARM_COVQ *CriarInstancia(SWARM_TYPE swarm_type, int N, int K, int qtd_part, double erro, SCALING scaleType, double scale)
{
    SWARM_COVQ *swarmCovq;

    switch (swarm_type)
    {
        case NO_SWARM:
            swarmCovq = new COVQ(N, K, erro, scaleType, scale);
            break;

        case FSS:
            swarmCovq = new FSS_COVQ(N, K, qtd_part, FSS, erro, scaleType, scale);
            break;

        case PSO:
            swarmCovq = new PSO_COVQ(N, K, qtd_part, PSO, erro, scaleType, scale);
            break;

        case FF:
            swarmCovq = new FF_COVQ(N, K, qtd_part, FF, erro, scaleType, scale);
            break;

        case FW:
            swarmCovq = new FW_COVQ(N, K, qtd_part, FW, erro, scaleType, scale);
            break;

        default:
            cout << endl << "Swarm technique not implemented!" << endl;
            exit(1);
    }

    return swarmCovq;
}

double *RodarMoeaCovq(SWARM_COVQ *swarm, const char *imagem, NNS nns, int simulacao)
{
    double **dic, **treino, **imagem_quantizada;
    double *resultados;
    double psnr_medio = 0, psnr;
    char nome[50], nome_recon[50];
    int *p_codificada, *p_transmitida;
    int nTreino, qtd_transmissoes = 20;
    int N = swarm->N, K = swarm->K;

    ofstream myfile;

    sprintf(nome, "resultados/b_%s_%.3f_%s.csv", imagem, swarm->erro, 
        converter_swarm_string(swarm->swarm_type));
    myfile.open(nome, std::ofstream::out | std::ofstream::app);

    if(simulacao == 1 && swarm->N == 32)
    {
        myfile << "N;PSNR;Iteracoes;Tempo;Distorcao\n";
    }

    sprintf(nome, "imagens/%s.dat", imagem);
    treino = ler_arquivo_treino(nome, &nTreino, K);

    sprintf(nome, "dics/%i_%s_%i.dic", N, imagem, simulacao);
    dic = ler_arquivo(nome, N * swarm->qtd_part, K);

    // Executar quantização
    swarm->nTreino = nTreino;
    resultados = swarm->Run(dic, treino, nns);

    // Codificar
    p_codificada = codificar(dic, treino, N, K, nTreino);

    // Transmitir pelo canal e medir o PNSR a cada transmissão
    for (int i = 0; i < qtd_transmissoes; i++)
    {
        cout << "Transmissao " << i + 1 << endl;

        p_transmitida = bsc(p_codificada, nTreino, swarm->erro);
        imagem_quantizada = decodificar(dic, p_transmitida, K, N, nTreino);

        psnr = CalcularPSNR(treino, imagem_quantizada, nTreino, K, 256);
        psnr_medio += psnr;

		sprintf(nome_recon, "reconstrucoes/%i_%.3f_%s_%s_%.2f.pgm", 
            swarm->N, swarm->erro, imagem, converter_swarm_string(swarm->swarm_type), psnr);
		//reconstruir_imagem(nome_recon, dic, p_transmitida, swarm->K, swarm->nTreino, 4, 256);

        free(p_transmitida);
        desaloca_matrizd(imagem_quantizada, nTreino);
    }

    resultados[3] = psnr_medio / (double)qtd_transmissoes;

    cout << endl << "Iteracoes: " << resultados[0] << endl;
    cout << "Duracao: " << resultados[1] << "s" << endl;
    cout << "Distorcao: " << resultados[2] << endl;
    cout << "PSNR = " << resultados[3] << endl << endl;

    myfile << swarm->N;
    myfile << SEPARADOR << resultados[3] << SEPARADOR << resultados[0];
    myfile << SEPARADOR << resultados[1] << SEPARADOR << resultados[2];
    myfile << "\n";

    // Desalocar memória
    desaloca_matrizd(dic, N * swarm->qtd_part);
    desaloca_matrizd(treino, nTreino);
    free(p_codificada);
    delete swarm;

    myfile.close();

    return resultados;
}
