//
//  definicoes.h
//  Util
//
//  Created by Felipe Ferreira on 23/12/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#ifndef definicoes_h
#define definicoes_h

#define SEG 128
#define SEPARADOR ';'

enum NNS
{
    BT,
    PDS,
    ENNS,
    IENNS,
    EENNS,
    IEENNS,
    EEENNS,
    DTA,
    IDTA
};

enum SCALING
{
    NONE,
    FIXED,
    VARIABLE
};

enum SWARM_TYPE
{
    NO_SWARM,
    FSS,
    PSO,
    FF,
    HBMO,
    FUZZY,
    FF_NH
};

struct Params
{
    double *mean_x;
    double *sum_xf;
    double *sum_xs;
    double *variance_x;
    double *norm_x;
    double *sum_squared_x;
    double *sum_x;
    double *max_x;
};

struct Fish
{
    double **dicionario;
    double **neighbor;
    double **deltaX;
    int *particao;
    double weight;
    double deltaD;
    double dist_atual;
    bool moved;
};

struct FireFly
{
    double **dicionario;
    int *particao;
    double dist_atual;
};

struct Particle
{
    double **dicionario;
    double **pbest;
    double **velocidade;
    double dist_atual;
    double dist_pbest;
    int *particao;
};

struct Firework
{
    double **position;
    double ***sparks;
    int qty_sparks;
    int *particao;
    double error;
};

struct Bee
{
    double **dicionario;
    double dist_atual;
    int *particao;
    bool isQueen;
};

#endif /* definicoes_h */
