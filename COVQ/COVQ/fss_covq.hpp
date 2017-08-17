//
//  fss_covq.hpp
//  COVQ
//
//  Versão do FSS acomodado às iterações do algoritmo COVQ.
//
//  Created by Felipe Ferreira on 17/12/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#ifndef _FSS_COVQ_H_
#define _FSS_COVQ_H_

#include <thread>
#include <chrono>
#include "swarm_covq.hpp"

using namespace std;

class FSS_COVQ : public SWARM_COVQ
{
public:
    double stepInd = 0.4;
    double stepIndInitial = stepInd;
    double stepIndFinal = 0.1;
    double stepVol = 1.2;
    double stepVolInitial = stepVol;
    double stepVolFinal = 0.5;
    double weightScale = 5000;
    
    double * Run(double **dic, double **treino, NNS nns);
    
    FSS_COVQ();
    FSS_COVQ(int N, int K, int qtd_part, double erro, SCALING scaleType, double scale);
    
private:
    Fish * InitiateSchool(double **dic);
    void DealocateSchool(Fish *school);
    
    void FSSUpdate(Fish *school, double **treino, int niter);
    
    double ApplyIndividualOperator(Fish *school, double **treino);
    double ApplyFeedingOperator(Fish *school, double maxDeltaD, double *pastWeight);
    void ApplyInstinctiveMovement(Fish *school);
    void ApplyVolitiveMoviment(Fish *school, double **treino, double currentWeight, double pastWeight);
    void ApplyBreeding(Fish *school);
    
    double CalcularDistanciaEuclidiana(double *x, double *y);
    
    void UpdateSteps(int cycle);
};

FSS_COVQ::FSS_COVQ() : SWARM_COVQ()
{
    
}

FSS_COVQ::FSS_COVQ(int N, int K, int qtd_part, double erro, SCALING scaleType, double scale) : SWARM_COVQ(N, K, qtd_part, erro, scaleType, scale)
{
    
}

//Retorno:
//    [0] - Quantidade de iterações realizadas
//    [1] - Duração do algoritmo em segundos
double * FSS_COVQ::Run(double **dic, double **treino, NNS nns)
{
    int niter;
    double *result;
    double limiar, dist_atual, dist_anterior, reducao;
    Fish *school;
    
    result = aloca_vetor_double(4);
    school = InitiateSchool(dic);
    CalcularParamTreino(treino, nns);
    
    limiar = 0.001;
    dist_anterior = 1e20;
    niter = 1;
    
    cout << "Inicio FSS" << endl;
    
    auto t_start = std::chrono::high_resolution_clock::now();
    
    while(true)
    {
        cout << "Iteracao = " << niter << endl;
        
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, school[i].dicionario, school[i].particao, nns);
        }
        
        for(int i = 0; i < qtd_part; i++)
        {
            school[i].dist_atual = CalcularDistorcaoMedia(treino, school[i].dicionario, school[i].particao);
            
            if(school[i].dist_atual < dist_gbest)
            {
                dist_gbest = school[i].dist_atual;
                copy(gbest, school[i].dicionario, N, K);
            }
        }
        
        dist_atual = dist_gbest;
        reducao = ((dist_anterior - dist_atual)/dist_anterior);
        
        cout << "Reducao Percentual = " << reducao << endl;
        cout << "Distorcao = " << dist_gbest << endl;
        
        if(dist_atual == 0 || fabs(reducao) < limiar) break;
        
        dist_anterior = dist_atual;
        
        FSSUpdate(school, treino, niter);
        
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, school[i].dicionario, school[i].particao, nns);
        }
        
        for(int i = 0; i < qtd_part; i++)
        {
            AtualizarCentroides(treino, school[i].dicionario, school[i].particao, niter);
        }

        niter++;
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    
    result[0] = niter;
    result[1] = std::chrono::duration<double, std::milli>(t_end-t_start).count() / 1000.0;
    result[2] = dist_gbest;
    
    DealocateSchool(school);
    
    if (nns != BT && nns != PDS)
        DesalocarParams(nns);
    
    copy(dic, gbest, N, K);
    
    return result;
}

void FSS_COVQ::FSSUpdate(Fish *school, double **treino, int niter)
{
    double pastWeight, currentWeight, maxDeltaD;

    maxDeltaD = ApplyIndividualOperator(school, treino);
    currentWeight = ApplyFeedingOperator(school, maxDeltaD, &pastWeight);
    ApplyInstinctiveMovement(school);
    ApplyVolitiveMoviment(school, treino, currentWeight, pastWeight);
    
    //ApplyBreeding(school);
    
    UpdateSteps(niter - 1);
}

Fish * FSS_COVQ::InitiateSchool(double **dic)
{
    Fish *school = new Fish[qtd_part];
    
    for(int i = 0; i < qtd_part; i++)
    {
        school[i].dicionario = aloca_matrizd(N, K);
        school[i].neighbor = aloca_matrizd(N, K);
        school[i].deltaX = aloca_matrizd(N, K);
        school[i].particao = aloca_vetor_int(nTreino);
        school[i].weight = weightScale/2;
        school[i].deltaD = 0;
        
        for(int j = 0; j < N; j++)
        {
            for(int k = 0; k < K; k++)
            {
                school[i].dicionario[j][k] = dic[(i * N) + j][k];
            }
        }
    }
    
    gbest = aloca_matrizd(N, K);
    dist_gbest = 1e20;
    
    return school;
}

void FSS_COVQ::DealocateSchool(Fish *school)
{
    for(int i = 0; i < qtd_part; i++)
    {
        desaloca_matrizd(school[i].dicionario, N);
        desaloca_matrizd(school[i].neighbor, N);
        desaloca_matrizd(school[i].deltaX, N);
        free(school[i].particao);
    }
    
    delete school;
}

double FSS_COVQ::ApplyIndividualOperator(Fish *school, double **treino)
{
    double prox_dist, maxDeltaD = 0;
    int i, j, k;

    for(i = 0; i < qtd_part; i++)
    {
        for(j = 0; j < N; j++)
        {
            for(k = 0; k < K; k++)
            {
                school[i].deltaX[j][k] = unif_neg()*stepInd;
                school[i].neighbor[j][k] = school[i].dicionario[j][k] + school[i].deltaX[j][k];
                
                if(school[i].neighbor[j][k] < 0)
                {
                    school[i].neighbor[j][k] = 0;
                }
                
                else if(school[i].neighbor[j][k] > 255)
                {
                    school[i].neighbor[j][k] = 255;
                }
            }
        }
        
        prox_dist = CalcularDistorcaoMedia(treino, school[i].neighbor, school[i].particao);
        
        school[i].deltaD = school[i].dist_atual - prox_dist;
        
        if(abs(school[i].deltaD) > maxDeltaD)
        {
            maxDeltaD = abs(school[i].deltaD);
        }
        
        school[i].moved = false;
        
        if(prox_dist < school[i].dist_atual)
        {
            school[i].moved = true;
        }
        
        copy(school[i].dicionario, school[i].neighbor, N, K);
        school[i].dist_atual = prox_dist;
        
        if(school[i].dist_atual < dist_gbest)
        {
            copy(gbest, school[i].dicionario, N, K);
            dist_gbest = school[i].dist_atual;
        }
    }

    return maxDeltaD;
}

double FSS_COVQ::ApplyFeedingOperator(Fish *school, double maxDeltaD, double *pastWeight)
{
    double currentWeight = 0;
    
    *pastWeight = 0;
    
    for(int i = 0; i < qtd_part; i++)
    {
        *pastWeight += school[i].weight;
        
        school[i].weight = school[i].weight + (school[i].deltaD/maxDeltaD);
        
        if(school[i].weight < 1)
        {
            school[i].weight = 1;
        }
        else if(school[i].weight > weightScale)
        {
            school[i].weight = weightScale;
        }
        
        currentWeight += school[i].weight;
    }
    
    return currentWeight;
}

void FSS_COVQ::ApplyInstinctiveMovement(Fish *school)
{
    int i, j, k;
    double **sumDeltaX, **direction;
    double sumDeltaD = 0;
    
    sumDeltaX = aloca_matrizd(N, K);
    direction = aloca_matrizd(N, K);
    
    for(i = 0; i < qtd_part; i++)
    {
        if(school[i].moved)
        {
            sumDeltaD += school[i].deltaD;
            
            for(j = 0; j < N; j++)
            {
                for(k = 0; k < K; k++)
                {
                    sumDeltaX[j][k] += school[i].deltaX[j][k]*school[i].deltaD;
                }
            }
        }
    }
    
    for(j = 0; j < N; j++)
    {
        for(k = 0; k < K; k++)
        {
            if(sumDeltaD != 0)
            {
                direction[j][k] = sumDeltaX[j][k]/sumDeltaD;
            }
            else
            {
                direction[j][k] = 0;
            }
        }
    }
    
    for(i = 0; i < qtd_part; i++)
    {
        for(j = 0; j < N; j++)
        {
            for(k = 0; k < K; k++)
            {
                school[i].neighbor[j][k] = school[i].dicionario[j][k] + direction[j][k];
            }
        }
    }
    
    desaloca_matrizd(sumDeltaX, N);
    desaloca_matrizd(direction, N);
}

void FSS_COVQ::ApplyVolitiveMoviment(Fish *school, double **treino, double currentWeight, double pastWeight)
{
    int i, j, k, direction;
    double distance, prox_dist, sumWeight = 0;
    double **baricenter, **sumPosition;
    
    baricenter = aloca_matrizd(N, K);
    sumPosition = aloca_matrizd(N, K);
    
    for(i = 0; i < qtd_part; i++)
    {
        sumWeight += school[i].weight;
        
        for(j = 0; j < N; j++)
        {
            for(k = 0; k < K; k++)
            {
                sumPosition[j][k] = school[i].neighbor[j][k]*school[i].weight;
            }
        }
    }
    
    for(j = 0; j < N; j++)
    {
        for(k = 0; k < K; k++)
        {
            baricenter[j][k] = sumPosition[j][k]/sumWeight;
        }
    }
    
    if(currentWeight > pastWeight)
    {
        direction = 1;
    }
    else
    {
        direction = -1;
    }
    
    for(i = 0; i < qtd_part; i++)
    {
        for(j = 0; j < N; j++)
        {
            distance = CalcularDistanciaEuclidiana(school[i].neighbor[j], baricenter[j]);
            
            if(distance != 0)
            {
                for(k = 0; k < K; k++)
                {
                    school[i].neighbor[j][k] +=
                        (direction*stepVol*unif()*(school[i].neighbor[j][k] - baricenter[j][k]))/distance;
                    
                    if(school[i].neighbor[j][k] < 0)
                    {
                        school[i].neighbor[j][k] = 0;
                    }
                    
                    else if(school[i].neighbor[j][k] > 255)
                    {
                        school[i].neighbor[j][k] = 255;
                    }
                }
            }
        }
        
        prox_dist = CalcularDistorcaoMedia(treino, school[i].neighbor, school[i].particao);
        
        copy(school[i].dicionario, school[i].neighbor, N, K);
        school[i].dist_atual = prox_dist;
        
        if(school[i].dist_atual < dist_gbest)
        {
            copy(gbest, school[i].dicionario, N, K);
            dist_gbest = school[i].dist_atual;
        }
    }
    
    desaloca_matrizd(baricenter, N);
    desaloca_matrizd(sumPosition, N);
}

void FSS_COVQ::ApplyBreeding(Fish *school)
{
    int j, k;
    int index_rand, index_rand2;
    int worst_index = 0, best_index = 0;
    double worst_dist = 0, best_dist = 1e20;
    bool busca;

    for (j = 0; j < qtd_part; j++)
    {
        if (school[j].dist_atual > worst_dist)
        {
            worst_dist = school[j].dist_atual;
            worst_index = j;
        }

        if (school[j].dist_atual < best_dist)
        {
            best_dist = school[j].dist_atual;
            best_index = j;
        }
    }

    do
    {
        busca = false;
        index_rand = rand()%qtd_part;

        if (index_rand != best_index && index_rand != worst_index)
        {
            busca = true;
        }
    }
    while (busca != true);

    do
    {
        busca = false;
        index_rand2 = rand()%qtd_part;

        if (index_rand2 != best_index && index_rand2 != worst_index && index_rand2 != index_rand)
        {
            busca = true;
        }
    }
    while (busca != true);

    for(j = 0; j < N; j++)
    {
        for(k = 0; k < K; k++)
        {
            school[worst_index].dicionario[j][k] = (school[best_index].dicionario[j][k] + school[worst_index].dicionario[j][k])/2;
            school[index_rand].dicionario[j][k] = (school[index_rand2].dicionario[j][k] + school[index_rand].dicionario[j][k])/2;
            school[index_rand2].dicionario[j][k] = (school[best_index].dicionario[j][k] + school[index_rand2].dicionario[j][k])/2;
        }
    }

    school[worst_index].weight = (school[best_index].weight + school[worst_index].weight)/2;
    school[index_rand].weight = (school[index_rand2].weight + school[index_rand].weight)/2;
    school[index_rand2].weight = (school[best_index].weight + school[index_rand2].weight)/2;
}

double FSS_COVQ::CalcularDistanciaEuclidiana(double *x, double *y)
{
    double distance = 0;
    
    for(int i = 0; i < K; i++)
    {
        distance += (x[i] - y[i]) * (x[i] - y[i]);
    }
    
    return sqrt(distance);
}

void FSS_COVQ::UpdateSteps(int cycle)
{
    stepInd = stepIndFinal + stepIndInitial/(stepIndInitial + cycle);
    stepVol = stepVolFinal + stepVolInitial/(stepVolInitial + cycle);
}

#endif /* _FSS_COVQ_H_ */
