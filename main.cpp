#include "include/Variables.h"
#include "include/TimeScheme.h"
#include "include/Quadrature.h"
#include "include/Orthogonal.h"
#include "include/limiter.h"
#include "include/Flux.h"
#include "include/InitCondition.h"
#include<fstream>

double dx(0.1), dt(0.001), L(10),gama(1.4);
int ord(1), r(3);
std::map<int,std::vector<double>> var;

int main()
{

    std::vector<double> sol(3*ord*std::floor(L/dx),0);

    


    sol=Quadrature0();

    


    //Resolution 

    int n=1000;

    for(int i=0;i<n;i++)
    {
        sol=Update(sol,1);
    }
    
    std::ofstream myfile;
    myfile.open ("example.txt");
    std::vector<double> sol0(ord*std::floor(L/dx),0);
    sol0=Quadrature0();
    for(int i=0;i<std::floor(L/dx);i++)
    {
        myfile <<i*dx<<" "<<transformp(sol,i*dx)[RHO].at(i)<<" "<<transformp(sol,i*dx)[U].at(i)<<" "<<transformp(sol,i*dx)[E].at(i)<<" "<<transformp(sol,i*dx)[P].at(i) <<std::endl;
    }
   
    myfile.close();




    return 0;
}