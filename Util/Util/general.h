//
//  general.h
//  Util
//
//  Created by Felipe Ferreira on 15/12/15.
//  Copyright © 2015 Felipe Ferreira. All rights reserved.
//

#ifndef general_h
#define general_h

#include <cstdlib>
#include <ctime>
#include <cmath>
#include "inout.h"
#include <cstring>
#include <sstream>

using namespace std;

//DEFINIÇÕES

bool str_isequal(const char *s1, const char *s2);

void copy(double **receiver, double **sender, int x, int y);
void copy(double *receiver, double *sender, int x);

namespace Util
{
    void Copy(int **receiver, int **sender, int x, int y);
    void Copy(int *receiver, int *sender, int x);
}

template <typename T> string to_string(T value);
string FormatNumber(double number);

//IMPLEMETAÇÕES

bool str_isequal(const char *s1, const char *s2)
{
    if (strlen(s1) != strlen(s2))
    {
        return false;
    }
    else
    {
        int len = int(strlen(s1));
        bool equal = true;

        for (int i = 0; i < len; i++)
        {
            if (s1[i] != s2[i])
            {
                return false;
            }
        }

        return equal;
    }
}

void copy(double **receiver, double **sender, int x, int y)
{
    int i, j;

    for (i = 0; i < x; i++)
    {
        for (j = 0; j < y; j++)
        {
            receiver[i][j] = sender[i][j];
        }
    }
}

void copy(double *receiver, double *sender, int x)
{
    int i;

    for (i = 0; i < x; i++)
    {
        receiver[i] = sender[i];
    }
}

void Util::Copy(int **receiver, int **sender, int x, int y)
{
    int i, j;

    for (i = 0; i < x; i++)
    {
        for (j = 0; j < y; j++)
        {
            receiver[i][j] = sender[i][j];
        }
    }
}

void Util::Copy(int *receiver, int *sender, int x)
{
    int i;

    for (i = 0; i < x; i++)
    {
        receiver[i] = sender[i];
    }
}

template <typename T> string to_string(T value)
{
    ostringstream os;
    os << value;
    return os.str();
}

string FormatNumber(double number)
{
    string result;

    result = to_string(number);

    for (int i = 0; i < result.size(); i++)
    {
        if (result[i] == '.')
        {
            result[i] = ',';
            break;
        }
    }

    return result;
}

#endif /* general_h */
