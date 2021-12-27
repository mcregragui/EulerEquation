#include<iostream>
#include<vector>
#include<map>

#define RHO 0
#define U  1
#define E 2
#define P 3

//extern std::vector<double> RHO, U, P, E;
extern std::map<int,std::vector<double>> var;
extern double gamma;
extern double dx, dt, L;
extern int ord, r;