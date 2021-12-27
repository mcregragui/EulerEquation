#include<iostream>
#include<cmath>
#include<vector>
#include<map>
#include "Variables.h"
#include "Orthogonal.h"
#include "Template.h"


double u0(int variable,double x);


std::map<int,std::vector<double>> transformp(std::vector<double> sol,double x);

std::map<int,std::vector<double>> FU(std::vector<double> sol,double x);