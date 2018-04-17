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

float BlackScholes(float S, float sigma, float t, float r, float K)
{
    float d1= 1/(sigma*sqrt(t))*(log(S/K)+(r+sigma*sigma/2)*t);
    float d2= 1/(sigma*sqrt(t))*(log(S/K)+(r-sigma*sigma/2)*t);

    return (S*GaussianLaw(d1)-K*exp(-r*t)*GaussianLaw(d2));
}

//European call pricing with Monte-Carlo
float montecarlo(int nbpaths,float S,float sigma,float t,float r,float K)
{

    float x=0;
    float St;
    float payoff;

    for (int i=0; i<nbpaths; i++){

        St= S * exp((r-sigma*sigma/2.0)*t + sigma*sqrt(t)*generateGaussianRandomValue());

        payoff= exp(-r*t) * std::max(St-K,float(0.));

        x += payoff;
    }

    float price= x/nbpaths;

    return price;
}

//European call pricing with Monte-Carlo (with variance computation)
float montecarlo_call(int nbpaths, float S, float sigma, float t, float r, float K, float &sig)
{
    float x=0, x2=0;
    float St;
    float payoff;

    for (int i=0; i<nbpaths; i++){

        St= S * exp((r-sigma*sigma/2.0)*t + sigma*sqrt(t)*generateGaussianRandomValue());

        payoff= exp(-r*t) * std::max(St-K,float(0.));

        x += payoff;
        x2 += payoff*payoff;
    }

    float price = x/nbpaths ;
    sig = x2/nbpaths - price*price;

    return price;
}

//Simulation of stock-path with Monte-Carlo
void simulation_stockpath(int N, float S, float sigma, float t, float r, float M, float &ST, float &maxi, float &avg)
{
    float stock[N+1];
    float x= S;
    stock[0]= S;
    maxi= S;

    for (int i=1; i<=N; i++)
    {
        stock[i]= stock[i-1] * exp((r-sigma*sigma/2.0)*t/N + sigma*sqrt(t/N)*generateGaussianRandomValue());
        x += stock[i];

        if (stock[i]>maxi)
        {
            maxi= stock[i];
        }
    }

    ST= stock[N];
    avg= x/(N+1);
}

//Piricing of asian call option with Monte-Carlo
float MC_asiatique(int nbpaths, int N, float S, float sigma, float t, float r, float K, float &sig)
{
    float x=0, x2=0;
    float payoff;
    float ST, maxi, avg;

    for (int i=0; i<nbpaths; i++)
    {
        simulation_stockpath(N,S,sigma,t,r,0,ST,maxi,avg);

        payoff= exp(-r*t) * std::max(avg-K,float(0.));

        x += payoff;
        x2 += payoff*payoff;
    }

    float price= x/nbpaths;
    sig= x2/nbpaths - price*price;

    return price;
}

//Monte-Carlo for European Call with variance reduction technique (antithetic variables)
float MC_call_antithetique(int nbpaths, float S, float sigma, float t, float r, float K, float &sig)
{
    float x=0, x1=0, x2=0, sum1=0, sum2=0, sum_cov=0;
    float St_1, St_2, payoff1, payoff2;

    for (int i=0; i<nbpaths; i++)
    {
        float Z= generateGaussianRandomValue();

        St_1= S * exp((r-sigma*sigma/2.0)*t + sigma*sqrt(t)*Z);
        payoff1= exp(-r*t) * std::max(St_1-K,float(0.));

        St_2= S * exp((r-sigma*sigma/2.0)*t + sigma*sqrt(t)*(-Z));
        payoff2= exp(-r*t) * std::max(St_2-K,float(0.));

        x += 0.5*(payoff1+payoff2);

        x1 += payoff1;
        sum1 += payoff1*payoff1;

        x2 += payoff2;
        sum2 += payoff2*payoff2;

        sum_cov += payoff1*payoff2;
    }

    float var_1= sum1/nbpaths - (x1/nbpaths)*(x1/nbpaths);
    float var_2= sum2/nbpaths - (x2/nbpaths)*(x2/nbpaths);
    float cov_1_2= sum_cov/nbpaths - (x1/nbpaths)*(x2/nbpaths);

    sig= 0.25*(var_1 + var_2 + 2*cov_1_2);
    float price= x/nbpaths;

    return price;
}

//Stock-path simulation with antithetic variables technique
void simulation_stockpath_antithetique(int N, float S, float sigma, float t, float r, float &avg1, float &avg2)
{

    float stock1[N], stock2[N];
    float sum1= S, sum2= S;
    stock1[0]= S;
    stock2[0]= S;

    for (int i=1; i<N; i++){

        float Z= generateGaussianRandomValue();

        stock1[i] = stock1[i-1] * exp((r-sigma*sigma/2.0)*t/N + sigma*sqrt(t/N)*Z);
        sum1 += stock1[i];

        stock2[i] = stock2[i-1] * exp((r-sigma*sigma/2.0)*t/N + sigma*sqrt(t/N)*(-Z));
        sum2 += stock2[i];
    }

    avg1 = sum1/N;
    avg2= sum2/N;
}

//Asian option pricing with Monte-Carlo with antithetic variables
float MC_asiatique_antithetique(int nbpaths, int N, float S, float sigma, float t, float r, float K, float &sig)
{
    float x=0, x2=0;
    float S1, S2, payoff1, payoff2;

    for (int i=0; i<nbpaths; i++)
    {
        simulation_stockpath_antithetique(N,S,sigma,t,r,S1,S2);

        payoff1= exp(-r*t) * std::max(S1-K,float(0.));
        payoff2= exp(-r*t) * std::max(S2-K,float(0.));

        x += 0.5*(payoff1+payoff2);
        x2 += (0.5*(payoff1+payoff2))*(0.5*(payoff1+payoff2));
    }

    float price= x/nbpaths;
    sig= x2/nbpaths - price*price;

    return price;
}

//Asian call option pricing with Monte-Carlo using control variable Z=ST
float MC_asiatique_controle1(int nbpaths, int N, float S, float sigma, float t, float r, float K, float &sig)
{
    float x=0, x2=0, sum_payoff=0, sum_cov=0;
    float payoff, ST, maxi, avg;

    /* Compute Cov(h(X),Z) by Monte-Carlo, with h(X)=payoff and Z=ST */
    for (int j=0; j<1000; j++)
    {
        simulation_stockpath(N,S,sigma,t,r,0,ST,maxi,avg);

        payoff= exp(-r*t) * std::max(avg-K,float(0.));

        sum_payoff += payoff;
        sum_cov += payoff*ST;
    }

    float cov= (sum_cov/1000) - (S*exp(r*t))*(sum_payoff/1000);

    /* Compute Var(Z) with Z following a log-normal distribution (mu, s2) */
    float mu= log(S) + (r-sigma*sigma/2)*t;
    float s2= sigma*sigma*t;
    float varZ= (exp(s2) - 1) * exp(2*mu + s2);

    float c= -cov/varZ;

    /* Compute price and variance of option */
    for (int i=0; i<nbpaths; i++)
    {
        simulation_stockpath(N,S,sigma,t,r,0,ST,maxi,avg);

        payoff = exp(-r*t) * std::max(avg-K,float(0.));

        x += (payoff + c*(ST-S*exp(r*t)));
        x2 += (payoff + c*(ST-S*exp(r*t)))*(payoff + c*(ST-S*exp(r*t)));
    }

    float price= x/nbpaths;
    sig= x2/nbpaths - price*price;

    return price;
}

//Asian call option pricing with Monte-Carlo using control variable Z=AVERAGE(S)
float MC_asiatique_controle2(int nbpaths, int N, float S, float sigma, float t, float r, float K, float &sig)
{
    float x=0, x2=0, sum_esp=0, sum_payoff=0, sum_avg=0, sum_cov=0, sum_var=0;
    float payoff, ST, maxi, avg;

    /* Compute E(Z) */
    for (int k=1; k<=N; k++)
    {
        sum_esp += exp(r*k*t/N);
    }
    float espZ= S*(1+sum_esp)/(N+1);

    /* Compute cov(h(X),Z) by Monte-Carlo, with h(X)=payoff and Z=avg */
    for (int j=0; j<10000 ; j++)
    {
        simulation_stockpath(N,S,sigma,t,r,0,ST,maxi,avg);

        payoff= exp(-r*t) * std::max(avg-K,float(0.));

        sum_payoff += payoff;
        sum_avg += avg;
        sum_cov += avg*payoff;
        sum_var += avg*avg;
    }

    float cov= sum_cov/10000 - (sum_payoff/10000)*espZ;
    float varZ= sum_var/10000 - espZ*espZ;

    float c= -cov/varZ;

    /* Compute price and variance of option */
    for (int i=0; i<nbpaths; i++)
    {
        simulation_stockpath(N,S,sigma,t,r,0,ST,maxi,avg);

        payoff= exp(-r*t) * std::max(avg-K,float(0.));

        x += (payoff + c*(avg-espZ));
        x2 += (payoff + c*(avg-espZ))*(payoff + c*(avg-espZ));

    }

    float price= x/nbpaths;
    sig= x2/nbpaths - price*price;

    return price;
}

//Asian call option pricing with Monte-Carlo using control variable Z= MAX(ST-K,0)
float MC_asiatique_controle3(int nbpaths, int N, float S, float sigma, float t, float r, float K, float &sig)
{
    float x=0, x2=0, sum_payoff=0, sum_cov=0, sum_var=0;
    float payoff, ST, maxi, avg, Z;

    /* Calcul de cov(h(X),Z) par Monte-Carlo, avec h(X)=payoff et Z=MAX(ST-K,0) */
    for (int j=0; j<100000; j++)
    {
        simulation_stockpath(N,S,sigma,t,r,0,ST,maxi,avg);

        payoff= exp(-r*t) * std::max(avg-K,float(0.));
        Z= exp(-r*t) * std::max(ST-K,float(0.));

        sum_payoff += payoff;
        sum_cov += Z*payoff;
        sum_var += Z*Z;
    }

    float espZ= BlackScholes(S,sigma,t,r,K);
    float varZ= (sum_var/100000) - espZ*espZ;
    float cov= (sum_cov/100000) - espZ*(sum_payoff/100000);

    float c= -cov/varZ;

    /* Calcul du prix et de la variance de l'option asiatique */
    for (int i=0; i<nbpaths; i++)
    {
        simulation_stockpath(N,S,sigma,t,r,0,ST,maxi,avg);

        payoff= exp(-r*t) * std::max(avg-K,float(0.));
        Z= exp(-r*t) * std::max(ST-K,float(0.));

        x += (payoff + c*(Z-espZ));
        x2 += (payoff + c*(Z-espZ))*(payoff + c*(Z-espZ));
    }

    float price= x/nbpaths;
    sig = x2/nbpaths - price*price;

    return price;
}
