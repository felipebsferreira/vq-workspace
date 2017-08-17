//
//  covq.hpp
//  COVQ
//
//  Created by Felipe Ferreira on 17/12/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#ifndef _COVQ_H_
#define _COVQ_H_

#include "swarm_covq.hpp"

using namespace std;

class COVQ : public SWARM_COVQ
{
public:
    double * Run(double ** dic, double **treino, NNS nns);

    COVQ();
    COVQ(int N, int K, double erro, SCALING scaleType, double scale);
};

COVQ::COVQ() : SWARM_COVQ()
{

}

COVQ::COVQ(int N, int K, double erro, SCALING scaleType, double scale) : SWARM_COVQ(N, K, 1, erro, scaleType, scale)
{

}

//Retorno:
//    [0] - Quantidade de iterações realizadas
//    [1] - Duração do algoritmo em segundos
double * COVQ::Run(double ** dic, double **treino, NNS nns)
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

//    OrdenarDicionarioPorMedias(dic);

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

    fim = clock();

    result[0] = niter;
    result[1] = double(fim - inicio)/CLOCKS_PER_SEC;
    result[2] = dist_atual;

    free(particao);

    if (nns != BT && nns != PDS)
        DesalocarParams(nns);

    return result;
}

#endif /* _COVQ_H_ */
