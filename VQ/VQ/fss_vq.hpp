//
//  fss_vq.hpp
//  VQ
//
//  Created by Felipe Ferreira on 19/01/16.
//  Updated by Charles Fonseca on 25/01/16 - Inclusão da função de Breeding.
//  Copyright © 2016 Felipe Ferreira. All rights reserved.
//

#ifndef fss_vq_h
#define fss_vq_h

#include "swarm_vq.hpp"

using namespace std;

class FSS_VQ : public SWARM_VQ
{
public:
    double stepInd = 0.0001;
    double stepIndInitial = stepInd;
    double stepIndFinal = stepInd;
    double stepVol = 0.01;
    double stepVolInitial = stepVol;
    double stepVolFinal = stepVol;
    double weightScale = 5000;

    double * Run(double **dicionario, double **treino);

    FSS_VQ();
    FSS_VQ(int N, int K, int qtd_part, NNS nns);
    FSS_VQ(int N, int K, int qtd_part, NNS nns, SCALING scaleType, double scale);

private:
    double totalInd = 0;
    double bestInd = 0;
    double totalVol = 0;
    double bestVol = 0;

    Fish * InitiateSchool(double **dic);
    void DealocateSchool(Fish *school);

    void FSSUpdate(Fish *school, double **treino, int niter);

    double ApplyIndividualOperator(Fish *school, double **treino);
    double ApplyFeedingOperator(Fish *school, double maxDeltaD, double *pastWeight);
    void ApplyInstinctiveMovement(Fish *school);
    void ApplyVolitiveMoviment(Fish *school, double **treino, double currentWeight, double pastWeight);
    void ApplyBreeding(Fish *school);
    void ApplyOriginalBreeding(Fish *school);
    double CalcularDistanciaEuclidiana(double *x, double *y);

    void UpdateSteps(int cycle);
};

FSS_VQ::FSS_VQ() : SWARM_VQ()
{

}

FSS_VQ::FSS_VQ(int N, int K, int qtd_part, NNS nns) : SWARM_VQ(N, K, qtd_part, nns)
{

}

FSS_VQ::FSS_VQ(int N, int K, int qtd_part, NNS nns, SCALING scaleType, double scale) : SWARM_VQ(N, K, qtd_part, nns, scaleType, scale)
{

}

// Executar algoritmo LBG
double * FSS_VQ::Run(double **dicionario, double **treino)
{
    int *particao;
    int niter;
    double limiar;
    double distAtual, distAnterior, reducao;
    double *resultado;
    Fish *school;
    clock_t inicio, fim;

    particao = aloca_vetor_int(nTreino);
    resultado = aloca_vetor_double(4);
    school = InitiateSchool(dicionario);

    CalcularParamTreino(treino);

    limiar = 0.001;
    distAnterior = 1e20;
    niter = 1;

    cout << "Inicio FSS QV" << endl;

    inicio = clock();

    while(true)
    {
        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, school[i].dicionario, school[i].particao);
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

        distAtual = dist_gbest;
        reducao = ((distAnterior - distAtual)/distAnterior);

        if(distAtual == 0 || fabs(reducao) < limiar) break;
        distAnterior = distAtual;

        FSSUpdate(school, treino, niter);

        for(int i = 0; i < qtd_part; i++)
        {
            Classificar(treino, school[i].dicionario, school[i].particao);
        }

        for(int i = 0; i < qtd_part; i++)
        {
            AtualizarCentroides(treino, school[i].dicionario, school[i].particao, niter);
        }

        niter++;
    }

    fim = clock();

    resultado[0] = niter;
    resultado[1] = double(fim - inicio)/CLOCKS_PER_SEC;
    resultado[2] = distAtual;

    DealocateSchool(school);

    if (nns != BT && nns != PDS)
        DesalocarParams();

    copy(dicionario, gbest, N, K);

    return resultado;
}

void FSS_VQ::FSSUpdate(Fish *school, double **treino, int niter)
{
    double pastWeight, currentWeight, maxDeltaD;

    maxDeltaD = ApplyIndividualOperator(school, treino);
    currentWeight = ApplyFeedingOperator(school, maxDeltaD, &pastWeight);
    ApplyInstinctiveMovement(school);
    ApplyVolitiveMoviment(school, treino, currentWeight, pastWeight);

    ApplyBreeding(school);
    //ApplyOriginalBreeding(school);

    UpdateSteps(niter - 1);
}

Fish * FSS_VQ::InitiateSchool(double **dic)
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

void FSS_VQ::DealocateSchool(Fish *school)
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

double FSS_VQ::ApplyIndividualOperator(Fish *school, double **treino)
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
        totalInd++;
        if(prox_dist < school[i].dist_atual)
        {
            bestInd++;
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

double FSS_VQ::ApplyFeedingOperator(Fish *school, double maxDeltaD, double *pastWeight)
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

void FSS_VQ::ApplyInstinctiveMovement(Fish *school)
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

void FSS_VQ::ApplyVolitiveMoviment(Fish *school, double **treino, double currentWeight, double pastWeight)
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
                    school[i].neighbor[j][k] += (direction*stepVol*unif()*(school[i].neighbor[j][k] - baricenter[j][k]))/distance;

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

        copy(school[i].dicionario, school[i].neighbor, N, K); totalVol++; if(prox_dist < school[i].dist_atual) bestVol++;
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

void FSS_VQ::ApplyOriginalBreeding(Fish *school)
{
    int j, k, i;
    int worst_index = 0, best_index = 0, chosen = 0;
    double worst_dist = 0, best_dist = 1e20;
    double rate = 0, rate_aux, distance;

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

    for (j = 0; j < qtd_part; j++)
    {
        if (j != best_index)
        {
            rate_aux = school[j].weight;
            distance = 0;

            for (i = 0; i < N; i++)
            {
                for (k = 0; k < K; k++)
                {
                    distance += (school[j].dicionario[i][k] - school[best_index].dicionario[i][k])
                        * (school[j].dicionario[i][k] - school[best_index].dicionario[i][k]);
                }
            }

            distance = distance/(N);

            rate_aux = rate_aux / distance;

            if (rate_aux > rate)
            {
                chosen = j;
                rate = rate_aux;
            }
        }
    }

    for(j = 0; j < N; j++)
    {
        for(k = 0; k < K; k++)
        {
            school[worst_index].dicionario[j][k] = (school[best_index].dicionario[j][k] + school[chosen].dicionario[j][k])/2;
        }
    }

    school[worst_index].weight = (school[best_index].weight + school[chosen].weight)/2;
}

void FSS_VQ::ApplyBreeding(Fish *school)
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

double FSS_VQ::CalcularDistanciaEuclidiana(double *x, double *y)
{
    double distance = 0;

    for(int i = 0; i < K; i++)
    {
        distance += (x[i] - y[i]) * (x[i] - y[i]);
    }

    return sqrt(distance);
}

void FSS_VQ::UpdateSteps(int cycle)
{
    if (stepIndInitial != stepIndFinal)
    {
        stepInd = stepIndFinal + stepIndInitial/(stepIndInitial + cycle);
    }

    if (stepVolInitial != stepVolFinal)
    {
        stepVol = stepVolFinal + stepVolInitial/(stepVolInitial + cycle);
    }
}

#endif /* fss_vq_h */
