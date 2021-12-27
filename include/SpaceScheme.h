#include<iostream>
#include<cmath>
#include<vector>
#include "Variables.h"
#include "Flux.h"
#include "Limiter.h"
#include "Quadrature.h"



//std::map<int,double> Numeric_flux(std::vector<double>& sol,int j,int c, int flux);


std::vector<double> LH(std::vector<double> sol, int flux);