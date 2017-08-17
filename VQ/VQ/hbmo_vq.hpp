//
//  hbmo_vq.hpp
//  VQ
//
//  Created by Felipe Ferreira on 6/16/16.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef hbmo_vq_h
#define hbmo_vq_h

#include "swarm_vq.hpp"

using namespace std;

class HBMO_VQ : public SWARM_VQ
{
public:
    double c1 = 0.1; //0.8         // Constante de aceleração cognitiva
    double c2 = 0.2; //0.2         // Constante de aceleração social
    double w = 0.7;                // Fator de inércia
    
    double * Run(double **dic, double **treino);
    
    HBMO_VQ();
    HBMO_VQ(int N, int K, int qtd_part, NNS nns);
    HBMO_VQ(int N, int K, int qtd_part, NNS nns, SCALING scaleType, double scale);
    
private:
    Bee * InitiateSwarm(double **dic);
    void DealocateSwarm(Bee *swarm);
    
    void HBMOUpdate(Bee *swarm, double **treino, int niter);
};

HBMO_VQ::HBMO_VQ() : SWARM_VQ()
{
    
}

HBMO_VQ::HBMO_VQ(int N, int K, int qtd_part, NNS nns) : SWARM_VQ(N, K, qtd_part, nns)
{
    
}

HBMO_VQ::HBMO_VQ(int N, int K, int qtd_part, NNS nns, SCALING scaleType, double scale) : SWARM_VQ(N, K, qtd_part, nns, scaleType, scale)
{
    
}

//Retorno:
//    [0] - Quantidade de iterações realizadas
//    [1] - Duração do algoritmo em segundos
double * HBMO_VQ::Run(double **dic, double **treino)
{
    int niter;
    double *result;
    double limiar, dist_atual, dist_anterior, reducao;
    clock_t inicio, fim;
    Bee *swarm;
    
    result = aloca_vetor_double(4);
    swarm = InitiateSwarm(dic);
    CalcularParamTreino(treino);
    
    limiar = 0.001;
    dist_anterior = 1e20;
    niter = 1;
    
    cout << "Inicio HBMO" << endl;
    
    inicio = clock();
    
    while(true)
    {
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, swarm[i].dicionario, swarm[i].particao);
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
        
        dist_atual = dist_gbest;
        reducao = ((dist_anterior - dist_atual)/dist_anterior);
        
        if(dist_atual == 0 || fabs(reducao) < limiar)
            break;
        
        dist_anterior = dist_atual;
        
        HBMOUpdate(swarm, treino, niter);
        
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, swarm[i].dicionario, swarm[i].particao);
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
        DesalocarParams();
    
    copy(dic, gbest, N, K);
    
    return result;
}

void HBMO_VQ::HBMOUpdate(Bee *swarm, double **treino, int niter)
{

}

Bee * HBMO_VQ::InitiateSwarm(double **dic)
{
    Bee *swarm = new Bee[qtd_part];
    
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

void HBMO_VQ::DealocateSwarm(Bee *swarm)
{
    for(int i = 0; i < qtd_part; i++)
    {
        desaloca_matrizd(swarm[i].dicionario, N);
        free(swarm[i].particao);
    }
    
    delete swarm;
}

#endif /* hbmo_vq_h */
