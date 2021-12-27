#include<iostream>
#include<cmath>
#include<vector>
#include "Variables.h"
#include "Flux.h"
#include "limiter.h"
#include "SpaceScheme.h"
#include "Template.h"


std::vector<double> Update(std::vector<double> sol, int flux);
