//
//  FF_COVQ_NH.hpp
//  COVQ
//
//  Created by Felipe Ferreira on 04/01/16.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef FF_COVQ_NH_H
#define FF_COVQ_NH_H

#include "swarm_covq.hpp"

using namespace std;

class FF_COVQ_NH : public SWARM_COVQ
{
public:
    double alpha = 0.01;
    double alpha_init = alpha;
    double alpha_final = alpha;
    double beta = 1;
    double beta_min = 0;
    double gama = 0.00001;

    double * Run(double **dic, double **treino, NNS nns);

    FF_COVQ_NH();
    FF_COVQ_NH(int N, int K, int qtd_part, double erro, SCALING scaleType, double scale);

private:

    FireFly * InitiateSwarm(double **dic);
    void DealocateSwarm(FireFly *swarm);

    void FFUpdate(FireFly *swarm, double **treino, int niter);

    void Move(FireFly *swarm, int index);

    double CalculateFireFlyDistance(double **f1, double **f2);
    void MoveFireFly(double **f1, double **f2, double attractiveness);
    void MoveFireFlyRandomly(FireFly f1);
    int FindBestFireFly(FireFly *swarm);

    void UpdateStep(int cycle);
};

FF_COVQ_NH::FF_COVQ_NH() : SWARM_COVQ()
{

}

FF_COVQ_NH::FF_COVQ_NH(int N, int K, int qtd_part, double erro, SCALING scaleType, double scale) : SWARM_COVQ(N, K, qtd_part, erro, scaleType, scale)
{

}

//Retorno:
//    [0] - Quantidade de iterações realizadas
//    [1] - Duração do algoritmo em segundos
double * FF_COVQ_NH::Run(double **dic, double **treino, NNS nns)
{
    double *result;
    int *particao;
    int niter;
    double limiar, dist_atual, dist_anterior, reducao;
    clock_t inicio, fim;

    result = aloca_vetor_double(4);
    particao = aloca_vetor_int(nTreino);
    CalcularParamTreino(treino, nns);

    limiar = 0.001;
    dist_anterior = 1e20;
    niter = 1;

    cout << "Inicio COVQ" << endl;

    inicio = clock();

    while(true)
    {
        Classificar(treino, dic, particao, nns);
        dist_atual = CalcularDistorcaoMedia(treino, dic, particao);

        reducao = ((dist_anterior - dist_atual)/dist_anterior);

        if(dist_atual == 0 || fabs(reducao) < limiar) break;
        dist_anterior = dist_atual;

        AtualizarCentroides(treino, dic, particao, niter);

        niter++;
    }

    free(particao);

    if (nns != BT && nns != PDS)
        DesalocarParams(nns);

    cout << "Fim COVQ" << endl;
    cout << "Inicio FireFly" << endl;

    FireFly *swarm;

    swarm = InitiateSwarm(dic);
    CalcularParamTreino(treino, nns);

    gbest = aloca_matrizd(N, K);
    dist_gbest = 1e20;

    dist_anterior = 1e20;
    niter = 0;

    while(niter < 100)
    {
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, swarm[i].dicionario, swarm[i].particao, nns);
        }

        for(int i = 0; i < qtd_part; i++)
        {
            swarm[i].dist_atual = CalcularDistorcaoMedia(treino, swarm[i].dicionario, swarm[i].particao);

            if(swarm[i].dist_atual < dist_gbest)
            {
                dist_gbest = swarm[i].dist_atual;
                copy(gbest, swarm[i].dicionario, N, K);
            }
        }

        FFUpdate(swarm, treino, niter);

        niter++;
    }

    fim = clock();

    cout << "Fim" << endl;

    result[0] = niter;
    result[1] = double(fim - inicio)/CLOCKS_PER_SEC;
    result[2] = dist_gbest;

    DealocateSwarm(swarm);

    if (nns != BT && nns != PDS)
        DesalocarParams(nns);

    copy(dic, gbest, N, K);

    return result;
}

void FF_COVQ_NH::FFUpdate(FireFly *swarm, double **treino, int niter)
{
    int bestIndex;

    UpdateStep(niter - 1);

    bestIndex = FindBestFireFly(swarm);
    Move(swarm, bestIndex);
}

void FF_COVQ_NH::Move(FireFly *swarm, int best_index)
{
    double distance, attractiveness;
    int index_rand;
    bool busca;

    for (int i = 0; i < qtd_part; i++)
    {
        if (i == best_index || swarm[i].dist_atual == swarm[best_index].dist_atual)
        {
            MoveFireFlyRandomly(swarm[i]);
        }
        else
        {
            busca = false;

            do
            {
                index_rand = rand()%qtd_part;

                if (swarm[i].dist_atual > swarm[index_rand].dist_atual)
                {
                    busca = true;
                }
            } while (busca != true);

            // Calcular distância entre o vagalume i e j
            distance = CalculateFireFlyDistance(swarm[i].dicionario, swarm[index_rand].dicionario);

            // Calcular atratividade
            attractiveness = (beta - beta_min) * exp(-gama * distance * distance) + beta_min;

            // Mover o vagalume i em direção ao vagalume j
            MoveFireFly(swarm[i].dicionario, swarm[index_rand].dicionario, attractiveness);
        }
    }
}

double FF_COVQ_NH::CalculateFireFlyDistance(double **f1, double **f2)
{
    double distance = 0;

    for (int i = 0; i < N; i++)
    {
        distance += CalcularDEQ(f1[i], f2[i]);
    }

    return sqrt(distance/N);
}

// Movimenta o vagalume f1 em direção ao vagalume f2
void FF_COVQ_NH::MoveFireFly(double **f1, double **f2, double attractiveness)
{
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < K; k++)
        {
            f1[i][k] = f1[i][k] + attractiveness * (f2[i][k] - f1[i][k]) + alpha * (unif() - 0.5);

            if(f1[i][k] < 0) f1[i][k] = 0;
            else if(f1[i][k] > 255) f1[i][k] = 255;
        }
    }
}

void FF_COVQ_NH::MoveFireFlyRandomly(FireFly f1)
{
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < K; k++)
        {
            f1.dicionario[i][k] = f1.dicionario[i][k] + alpha * (unif() - 0.5);

            if(f1.dicionario[i][k] < 0) f1.dicionario[i][k] = 0;
            else if(f1.dicionario[i][k] > 255) f1.dicionario[i][k] = 255;
        }
    }
}

int FF_COVQ_NH::FindBestFireFly(FireFly *swarm)
{
    int i, best_index = 0;
    double best_dist = 1e20;

    for (i = 0; i < qtd_part; i++)
    {
        if (swarm[i].dist_atual < best_dist)
        {
            best_index = i;
            best_dist = swarm[i].dist_atual;
        }
    }

    return best_index;
}

FireFly * FF_COVQ_NH::InitiateSwarm(double **dic)
{
    FireFly *swarm = new FireFly[qtd_part];

    for(int i = 0; i < qtd_part; i++)
    {
        swarm[i].dicionario = aloca_matrizd(N, K);
        swarm[i].particao = aloca_vetor_int(nTreino);

        for(int j = 0; j < N; j++)
        {
            for(int k = 0; k < K; k++)
            {
                swarm[i].dicionario[j][k] = dic[(i * N) + j][k];
            }
        }
    }

    return swarm;
}

void FF_COVQ_NH::DealocateSwarm(FireFly *swarm)
{
    for(int i = 0; i < qtd_part; i++)
    {
        desaloca_matrizd(swarm[i].dicionario, N);
        free(swarm[i].particao);
    }

    delete swarm;
}

void FF_COVQ_NH::UpdateStep(int cycle)
{
    if (alpha_init != alpha_final)
    {
        alpha = alpha_final + alpha_init/(alpha_init + cycle);
    }
}

#endif /* FF_COVQ_NH_H */
