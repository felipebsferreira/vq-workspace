//
//  ff_covq.hpp
//  COVQ
//
//  Created by Felipe Ferreira on 04/01/16.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef FF_COVQ_H
#define FF_COVQ_H

#include "swarm_covq.hpp"

using namespace std;

class FF_COVQ : public SWARM_COVQ
{
public:
    double alpha = 0.7;
    double alpha_init = alpha;
    double alpha_final = alpha;
    double beta = 0.4;
    double beta_min = 0;
    double gama = 0.00001;
    
    double * Run(double **dic, double **treino, NNS nns);
    
    FF_COVQ();
    FF_COVQ(int N, int K, int qtd_part, SWARM_TYPE swarm_type, double erro, SCALING scaleType, double scale);
    
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

FF_COVQ::FF_COVQ() : SWARM_COVQ()
{
    
}

FF_COVQ::FF_COVQ(int N, int K, int qtd_part, SWARM_TYPE swarm_type, double erro, SCALING scaleType, double scale) : SWARM_COVQ(N, K, qtd_part, swarm_type, erro, scaleType, scale)
{
    
}

//Retorno:
//    [0] - Quantidade de iterações realizadas
//    [1] - Duração do algoritmo em segundos
double * FF_COVQ::Run(double **dic, double **treino, NNS nns)
{
    int niter, index_gbest = 0;
    double *result;
    double limiar, dist_atual, dist_anterior, reducao;
    clock_t inicio, fim;
    FireFly *swarm;
    
    result = aloca_vetor_double(4);
    swarm = InitiateSwarm(dic);
    CalcularParamTreino(treino, nns);
    
    limiar = 0.001;
    dist_anterior = 1e20;
    niter = 1;
    
    cout << "Inicio FF" << endl;
    
    inicio = clock();

    while(true)
    {
        cout << "Iteracao = " << niter << endl;
        
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
                index_gbest = i + 1;
            }
        }
        
        dist_atual = dist_gbest;
        reducao = ((dist_anterior - dist_atual)/dist_anterior);
        
        cout << "Reducao Percentual = " << reducao << endl;
        cout << "Distorcao = " << dist_gbest << endl;
        
        if(dist_atual == 0 || fabs(reducao) < limiar) break;
        
        dist_anterior = dist_atual;
        
        FFUpdate(swarm, treino, niter);
        
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, swarm[i].dicionario, swarm[i].particao, nns);
        }
        
        for(int i = 0; i < qtd_part; i++)
        {
            AtualizarCentroides(treino, swarm[i].dicionario, swarm[i].particao, niter);
        }

        niter++;
    }
    
    fim = clock();
    
    result[0] = niter;
    result[1] = double(fim - inicio)/CLOCKS_PER_SEC;
    result[2] = dist_gbest;
    
    DealocateSwarm(swarm);
    
    if (nns != BT && nns != PDS)
        DesalocarParams(nns);
    
    copy(dic, gbest, N, K);
    
    return result;
}

void FF_COVQ::FFUpdate(FireFly *swarm, double **treino, int niter)
{
    int bestIndex;
    
    UpdateStep(niter - 1);
    
    bestIndex = FindBestFireFly(swarm);
    Move(swarm, bestIndex);
}

void FF_COVQ::Move(FireFly *swarm, int best_index)
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

double FF_COVQ::CalculateFireFlyDistance(double **f1, double **f2)
{
    double distance = 0;
    
    for (int i = 0; i < N; i++)
    {
        distance += CalcularDEQ(f1[i], f2[i]);
    }
    
    return sqrt(distance/N);
}

// Movimenta o vagalume f1 em direção ao vagalume f2
void FF_COVQ::MoveFireFly(double **f1, double **f2, double attractiveness)
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

void FF_COVQ::MoveFireFlyRandomly(FireFly f1)
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

int FF_COVQ::FindBestFireFly(FireFly *swarm)
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

FireFly * FF_COVQ::InitiateSwarm(double **dic)
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
    
    gbest = aloca_matrizd(N, K);
    dist_gbest = 1e20;
    
    return swarm;
}

void FF_COVQ::DealocateSwarm(FireFly *swarm)
{
    for(int i = 0; i < qtd_part; i++)
    {
        desaloca_matrizd(swarm[i].dicionario, N);
        free(swarm[i].particao);
    }
    
    delete swarm;
}

void FF_COVQ::UpdateStep(int cycle)
{
    if (alpha_init != alpha_final)
    {
        alpha = alpha_final + alpha_init/(alpha_init + cycle);
    }
}

#endif /* FF_COVQ_H */
