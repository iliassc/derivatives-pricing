#include "Header.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <list>
#include <algorithm>
#include <stdio.h>
#include <string.h>

using namespace std;

std::default_random_engine generator;

float GaussianLaw (float t)
{
  return erf(t/sqrt(2))/2 - erf(-1000000000)/2;
}

float generateGaussianRandomValue()
{
    std::normal_distribution<float> distribution(0.0, 1.0);
    float randomValue = distribution(generator);

    return randomValue;
}

float generateUniformRandomValue()
{
    std::uniform_real_distribution<float> distribution(0.0, 1.0);
    float randomValue = distribution(generator);

    return randomValue;
}

float uniformLaw(float inf, float sup)
{
    std::uniform_real_distribution<float> distribution(inf, sup);
    float randomValue = distribution(generator);

    return randomValue;
}

vector<float> generateGaussianRandomVector(int N)
{
    vector<float> x(N);

    for (int i=0; i<N; i++)
    {
        x[i]= generateGaussianRandomValue();
    }

    return x;
}

float norme2(std::vector<float> x)
{
    float a=0;

    for (int i=0; i<x.size(); i++)
    {
        a+= pow(x[i],2) ;
    }

    return sqrt(a);
}

float norme(std::vector<float> x, std::vector<float> y)
{//x and y must have same dimension
    float a= x.size();

    return std::max(norme2(x)/sqrt(a),norme2(y)/sqrt(a));
}

// SIMULATION OF GAUSSIAN 2D COPULA
vector<float> gaussian_copula_2d(float rho)
{
    vector<float> U(2);

    float X1= generateGaussianRandomValue();
    float X2= rho*X1 + sqrt(1-rho*rho)*generateGaussianRandomValue();

    U[0]= GaussianLaw(X1);
    U[1]= GaussianLaw(X2);

    return U;
}

// SIMULATION OF CLAYTON 2D COPULA
vector<float> clayton_copula_2d(float delta)
{
    vector<float> U(2);
    float S= 0;
    U[0]= generateUniformRandomValue();
    S= S+pow(U[0],-delta)-1;
    float V= generateUniformRandomValue();
    U[1]= pow((1+S)*pow(V,-delta/(1+delta))-S,-1/delta);

    return U;
}

// SIMULATION OF CLAYTON N-DIM COPULA
vector<float> clayton_copula_N(int N, float delta)
{
    vector<float> U(N);
    float S= 0;
    U[0]= generateUniformRandomValue();

    for (int i=1; i<N; i++)
    {
        S= S+pow(U[i-1],-delta)-1;
        float V= generateUniformRandomValue();
        U[i]= pow((1+S)*pow(V,(-delta)/(1+i*delta))-S,-1/delta);
    }

    return U;
}

// CORREL MATRIX rho "flat"
vector<vector<float>> Correl_Matrix(float rho,int N)
{
    vector<vector<float>> x(N, vector<float>(N,rho));

    for (int i=0; i<N; i++)
    {
        x[i][i]= 1;
    }

    return x;
}

// VECTOR L IN CHOLESKY DECOMPOSITION
vector<vector<float>> Decomposition_Cholesky(vector<vector<float>> A)
{
    int n= A.size();
    float sum1=0;
    float sum2=0;
    float sum3=0;
    vector<vector<float>> l(n, vector<float> (n));

    l[0][0]= sqrt(A[0][0]);

    for (int j=1; j<n; j++)
    {
        l[j][0]= A[j][0]/l[0][0];
    }

    for (int i=1; i<(n-1); i++)
    {
        for (int k=0; k<i; k++)
        {
            sum1 += pow(l[i][k], 2);
            l[i][i]= sqrt(A[i][i]-sum1);
        }

        for (int j=(i+1); j<n; j++)
        {
            for (int k=0; k<i; k++)
            {
                sum2 += l[j][k]*l[i][k];
                l[j][i]= (A[j][i]-sum2)/l[i][i];
            }
        }
    }

    for (int k=0; k<(n-1); k++)
    {
        sum3 += pow(l[n-1][k],2);
        l[n-1][n-1]= sqrt(A[n-1][n-1]-sum3);
    }

    return l;
}

vector<float> Matrix_Vector(vector<vector<float>> a, vector<float> x)
{
    int n= x.size();
    vector<float> y(n);

    for (int i=0; i<n; i++)
    {
        float s= 0;

        for (int j=0; j<n; j++)
        {
            s += a[i][j]*x[j];
        }
        y[i]= s;
    }

    return y;
}

// MONTE-CARLO PRICING OF FTD WITH GAUSSIAN COPULA
float FTD_gaussian(int N, float rho, float T, float lambda, float recovery, float r, int nbpaths, float &lambda_global)
{
    int x= 0;

    float Tdefault, price;

    vector<vector<float>> Correl= Correl_Matrix(rho,N);
    vector<vector<float>> L= Decomposition_Cholesky(Correl);

    for (int k=0; k<nbpaths; k++)
    {
        int a=0;

        //Simulation de copule gaussien en dimension N
        vector<float> Z= generateGaussianRandomVector(N);
        vector<float> Y= Matrix_Vector(L,Z);
        vector<float> U(N);
        for (int i=0; i<N; i++)
        {
            U[i]= GaussianLaw(Y[i]);
        }

        //Calcul des temps de defaut
        for (int j=0; j<N; j++)
        {
            Tdefault= -log(1-U[j])/lambda;

            if (Tdefault<=T)
            {//Un defaut a eu lieu avant maturite
                a= 1;
                break;
            }
        }

        x += a;
    }

    price= x*exp(-r*T)*(1-recovery)/nbpaths;
    lambda_global = -log(price)/T;

    return price;
}

// MONTE-CARLO PRICING OF FTD WITH CLAYTON COPULA
float FTD_clayton(int N, float delta, float T, float lambda, float recovery, float r, int nbpaths, float &lambda_global)
{
    int x= 0;

    float Tdefault, price;

    for (int k=0; k<nbpaths; k++)
    {
        int a=0;

        //Simulation de copule de Clayton en dimension N
        vector<float> U= clayton_copula_N(N,delta);

        //Calcul des temps de defaut
        for (int j=0; j<N; j++)
        {
            Tdefault= -log(1-U[j])/lambda;

            if (Tdefault<=T)
            {//Un defaut a eu lieu avant maturite
                a= 1;
                break;
            }
        }

        x += a;
    }

    price= x*exp(-r*T)*(1-recovery)/nbpaths;
    lambda_global = -log(price)/T;

    return price;
}
