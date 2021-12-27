#include<iostream>
#include<vector>
#include<cmath>
#include<map>

#include "Variables.h"
#include "Template.h"
#include "InitCondition.h"
#include "Limiter.h"




std::map<int,std::vector<double>> LLF(std::vector<double>& sol,int c);


std::map<int,std::vector<double>> S(std::map<int,std::vector<double>> left, std::map<int,std::vector<double>> wright,int c);

std::vector<double> maxum(std::map<int,std::vector<double>> left, std::map<int,std::vector<double>> wright);

std::map<int,std::vector<double>> Eigen(std::map<int,std::vector<double>> left);


std::map<int,std::vector<double>> Numeric_flux(std::vector<double>& sol,int c,int flux);