#include "include/Flux.h"



std::map<int,std::vector<double>> LLF(std::vector<double>& sol,int c)
{
    std::map<int,std::vector<double>> Flux, Left, Wright,FLeft, FWright, Minus;
    
    std::vector<double> l=Tleft(sol);
    
    std::vector<double> w=Twright(sol);
    
    Left=transform(l);
   
    Wright=transform(w);
   
    FLeft=transform(l);
   
    FWright=transform(w);
    
    Minus=S(FLeft,FWright, c);
    
    double p=(c==1)*1.0;
    double m=-(c==-1)*1.0;
    
    Flux[RHO].push_back(0.5*(FLeft[RHO].at(0)+FWright[RHO].at(0)+Minus[RHO].at(0)));

    Flux[U].push_back(0.5*(FLeft[U].at(0)+FWright[U].at(0)+Minus[U].at(0)));

    Flux[E].push_back(0.5*(FLeft[E].at(0)+FWright[E].at(0)+Minus[E].at(0)));

    for(int i=1;i<std::floor(L/dx)-1;i++)
    {

        Flux[RHO].push_back(0.5*(FLeft[RHO].at(i+m)+FWright[RHO].at(i+p)+Minus[RHO].at(i)));

        Flux[U].push_back(0.5*(FLeft[U].at(i+m)+FWright[U].at(i+p)+Minus[U].at(i)));

        Flux[E].push_back(0.5*(FLeft[E].at(i+m)+FWright[E].at(i+p)+Minus[E].at(i)));


    }
    Flux[RHO].push_back(0.5*(FLeft[RHO].at(std::floor(L/dx)-1)+FWright[RHO].at(std::floor(L/dx)-1)+Minus[RHO].at(std::floor(L/dx)-1)));

    Flux[U].push_back(0.5*(FLeft[U].at(std::floor(L/dx)-1)+FWright[U].at(std::floor(L/dx)-1)+Minus[U].at(std::floor(L/dx)-1)));

    Flux[E].push_back(0.5*(FLeft[E].at(std::floor(L/dx)-1)+FWright[E].at(std::floor(L/dx)-1)+Minus[E].at(std::floor(L/dx)-1)));



    return Flux;
}


std::map<int,std::vector<double>> S(std::map<int,std::vector<double>> left, std::map<int,std::vector<double>> wright,int c)
{

    std::map<int,std::vector<double>> minus,leigens,weigens;

    std::vector<double> eigens;

    
    eigens=maxum(left,wright);
    

    double a=0.0;

    double p=(c==1)*1.0;
    
    double m=-(c==-1)*1.0;
    

    minus[RHO].push_back(eigens[0]*(left[RHO][0]-wright[RHO][0]));
    minus[U].push_back(eigens[0]*(left[U][0]-wright[U][0]));
    minus[E].push_back(eigens[0]*(left[E][0]-wright[E][0]));
    for(int i=1;i<std::floor(L/dx)-1;i++)
    {
       
        minus[RHO].push_back(eigens[i]*(left[RHO][i+m]-wright[RHO][i+p]));

        minus[U].push_back(eigens[i]*(left[U][i+m]-wright[U][i+p]));

        minus[E].push_back(eigens[i]*(left[E][i+m]-wright[E][i+p]));
    }
    minus[RHO].push_back(eigens[std::floor(L/dx)-1]*(left[RHO][std::floor(L/dx)-1]-wright[RHO][std::floor(L/dx)-1]));
    minus[U].push_back(eigens[std::floor(L/dx)-1]*(left[U][std::floor(L/dx)-1]-wright[U][std::floor(L/dx)-1]));
    minus[E].push_back(eigens[std::floor(L/dx)-1]*(left[E][std::floor(L/dx)-1]-wright[E][std::floor(L/dx)-1]));

    return minus;
}

std::vector<double> maxum(std::map<int,std::vector<double>> left, std::map<int,std::vector<double>> wright)
{
    std::map<int,std::vector<double>> leigens,weigens;

    std::vector<double> eigens(std::floor(L/dx),0.0);

    leigens=Eigen(left);

    weigens=Eigen(wright);

    
    double max=0.0;
    for(int i=0;i<std::floor(L/dx);i++)
    {
     
        max=0.0;
        for(int j=1;j<4;j++)
        {
            
            if(max<fabs(leigens[j].at(i)))
            {
                max=fabs(leigens[j].at(i));
            }

            if(max<fabs(weigens[j].at(i)))
            {
                max=fabs(weigens[j].at(i));
            }
        }
        eigens[i]=max;

    }
 

    return eigens;
}




std::map<int,std::vector<double>> Eigen(std::map<int,std::vector<double>> left)
{

    std::map<int,std::vector<double>> Eigen;
    double a=0.0;
    double c=0.0;
    
    for(int i=0;i<std::floor(L/dx);i++)
    {
        a=left[U].at(i)/left[RHO].at(i);
        c=sqrt(gama*left[P].at(i)/left[RHO].at(i));
        Eigen[1].push_back(a-c);
        Eigen[2].push_back(a);
        Eigen[3].push_back(a+c);
    }

    return Eigen;
}


std::map<int,std::vector<double>> Numeric_flux(std::vector<double>& sol,int c,int flux)
{
    switch (flux)
    {
    case 1:

        return LLF(sol,c);
        break;
    
    default:
        std::cout<<"Such a flux is not yet coded"<<std::endl;
        break;
    }
}