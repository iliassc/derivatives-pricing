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

//Finite difference method V1 (explicit numerical scheme) for solving Black-Scholes EDP (vanilla call)
vector<float> EDP_schema_explicite(float S0, float K, float r, float T, float sigma, float L,int N,int M)
{
    std::vector<float> P(N+1);
    std::vector<float> Q(N+1);
    float h = 2.*L/N;
    float dt = T/M;

    //Useful coefficients
    float alpha= dt * (sigma*sigma/(2.*h)*(1./h + 1./2.) - r/(2.*h));
    float beta= 1 - dt * (r + sigma*sigma/(h*h));
    float gamma= dt * (sigma*sigma/(2.*h)*(1./h-1./2.) + r/(2.*h));

    //Payoff at maturity (initial condition)
    for (int i=0; i<=N; i++)
    {
        P[i]= std::max(exp(-L+i*h)-K,float(0.));
    }

    //Resolution
    for (int j=0; j<=M; j++)
    {
        Q[0]=0; //Border conditions

        for (int i=1; i<N; i++)
        {
            Q[i]= alpha*P[i+1] + beta*P[i] + gamma*P[i-1];
        }

        Q[N]= exp(L) - K*exp(-r*(j+1)*dt); //Border conditions

        P = Q;
    }

    return P;
}

//Finite difference method V2 (explicit numerical scheme) for solving Black-Scholes EDP (vanilla call)
vector<vector<float>> EDP_schema_explicite_2(float S0, float K, float r, float T, float sigma, float L,int N,int M)
{
    vector<vector<float>> P(M+1, vector<float>(N+1,0) );
    float h = 2.*L/N;
    float dt = T/M;

    //Useful coefficients
    float alpha= dt * (sigma*sigma/(2.*h)*(1./h + 1./2.) - r/(2.*h));
    float beta= 1 - dt * (r + sigma*sigma/(h*h));
    float gamma= dt * (sigma*sigma/(2.*h)*(1./h-1./2.) + r/(2.*h));

    //Payoff at maturity
    for (int i=0; i<=N; i++)
    {
        P[i][0]= std::max(exp(-L+i*h)-K,float(0.));
    }

    //Resolution
    for (int j=1; j<=M; j++)
    {
        P[0][j]=0; //Border conditions

        for (int i=1; i<N; i++)
        {
            P[i][j]= alpha*P[i+1][j-1] + beta*P[i][j-1] + gamma*P[i-1][j-1];
        }

        P[N][j]= exp(L) - K*exp(-r*(j+1)*dt); //Border conditions
    }

    return P;
}

//Finite difference method (implicit numerical scheme) for solving Black-Scholes EDP (vanilla call)
vector<float> EDP_schema_implicite(float S0, float K, float r, float T, float sigma, float L,int N,int M)
{
    std::vector<float> U(N+1);
    std::vector<float> V(N+1);
    std::vector<float> c(N+1);
    std::vector<float> d(N+1);
    float h= 2.*L/N;
    float dt= T/M;

    float alpha= -dt*(2*sigma*sigma + 2*h*r - sigma*sigma*h)/(4*h*h);
    float beta= -dt*(-sigma*sigma - r*h*h)/(h*h);
    float gamma= -dt*(2*sigma*sigma - 2*h*r + sigma*sigma*h)/(4*h*h);

    for (int i=0; i<=N; i++)
    {
        U[i]= std::max(exp(-L+i*h)-K,float(0.));
    }

    c[0]= alpha/(1+beta);

    for (int k=1; k<=N; k++)
    {
        c[k]= alpha/(1+beta-gamma*c[k-1]);
    }

    for (int j=0; j<=M; j++)
    {

        d[0]= U[0]/(1+beta);
        for (int k=1; k<=N; k++)
        {
            d[k]= (U[k]-gamma*d[k-1])/(1+beta-gamma*c[k-1]);
        }

        V[0]= 0;
        V[N]= exp(L)-K*exp(-r*(j+1)*dt);
        for (int i=N-1; i>=1; i--)
        {
            V[i]= d[i]-c[i]*V[i+1];
        }
        V[N]= exp(L)-K*exp(-r*(j+1)*dt);

        U= V;
    }

    return U;
}

float BS_call(float S0, float K, float r, float T, float sigma)
{
    float d1= 1/(sigma*sqrt(T))*(log(S0/K)+(r+0.5*sigma*sigma)*T);
    float d2= d1 - sigma*sqrt(T);

    return (S0*GaussianLaw(d1) - K*exp(-r*T)*GaussianLaw(d2));
}

vector<float> BS(float S0, float K, float r, float T, float sigma, float L,int N)
{
    float h = 2.*L/N;
    std::vector<float> P(N+1);

    for (int i=0; i < N+1; i++)
    {
        P[i]= BS_call(exp(-L+i*h),K,r,T,sigma);
    }

    return P;
}
