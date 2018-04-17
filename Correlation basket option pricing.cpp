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

//MONTE-CARLO FOR BEST-OF CALL
float montecarlo_correlation_max(float S, float K, float r, float t, float sigma, float rho, int nbpaths)
{
    float payoff, St_1, St_2;
    float x= 0;

    for (int i=0; i<nbpaths; i++)
    {
        float X1= generateGaussianRandomValue();
        float X2= generateGaussianRandomValue();

        float B1= X1;
        float B2= rho*X1 + sqrt(1-rho*rho)*X2;

        St_1= S * exp((r-sigma*sigma/2.0)*t + sigma*sqrt(t)*B1);
        St_2= S * exp((r-sigma*sigma/2.0)*t + sigma*sqrt(t)*B2);

        payoff= std::max(std::max(St_1,St_2)-K,float(0.));
        x += payoff;
    }

    float price= exp(-r*t)*x/nbpaths;

    return price;
}

//MONTE-CARLO OPTION CALL ON AVERAGE OF BASKET
float montecarlo_correlation_avg(float S, float K, float r, float t, float sigma, float rho, int nbpaths)
{
    float payoff, St_1, St_2;
    float x= 0;

    for (int i=0; i<nbpaths; i++)
    {
        float X1= generateGaussianRandomValue();
        float X2= generateGaussianRandomValue();

        float B1= X1;
        float B2= rho*X1 + sqrt(1-rho*rho)*X2;

        St_1= S * exp((r-sigma*sigma/2.0)*t + sigma*sqrt(t)*B1);
        St_2= S * exp((r-sigma*sigma/2.0)*t + sigma*sqrt(t)*B2);

        payoff= std::max(float(0.5*(St_1+St_2)-K),float(0.));
        x += payoff;
    }

    float price= exp(-r*t)*x/nbpaths;

    return price;
}

//MONTE-CARLO FOR WORST-OF CALL
float montecarlo_correlation_min(float S, float K, float r, float t, float sigma, float rho, int nbpaths)
{
    float payoff, St_1, St_2;
    float x= 0;

    for (int i=0; i<nbpaths; i++)
    {
        float X1= generateGaussianRandomValue();
        float X2= generateGaussianRandomValue();

        float B1= X1;
        float B2= rho*X1 + sqrt(1-rho*rho)*X2;

        St_1= S * exp((r-sigma*sigma/2.0)*t + sigma*sqrt(t)*B1);
        St_2= S * exp((r-sigma*sigma/2.0)*t + sigma*sqrt(t)*B2);

        payoff= std::max(std::min(St_1,St_2)-K,float(0.));
        x += payoff;
    }

    float price= exp(-r*t)*x/nbpaths;

    return price;
}
