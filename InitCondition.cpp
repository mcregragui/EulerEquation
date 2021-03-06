#include "include/InitCondition.h"






double u0(int variable,double x)
{
    if(x<L/2.0)
    {
        if(variable==RHO)
        {
            return 1.0;
        }
        if(variable==U)
        {
            return 0;
        }
        if(variable==E)
        {
            return 1.0/(gama-1);
        }
        if(variable==P)
        {
            return 1.0;
        }
    }
    else
    {
        if(variable==RHO)
        {
            return 0.125;
        }
        if(variable==U)
        {
            return 0;
        }
        if(variable==E)
        {
            return 0.1/(gama-1);
        }
        if(variable==P)
        {
            return 0.1;
        }
    }
}





std::map<int,std::vector<double>> transformp(std::vector<double> solf,double x)
{

    std::vector<double> sol(std::floor(L/dx),0);
    std::map<int,std::vector<double>> V;
    for(int i=0;i<std::floor(L/dx);i++)
    {
        
        
        if(ord==1)
        {
            V[RHO].push_back(solf[3*ord*i]);
            V[U].push_back(solf[3*ord*i+U]);
            V[E].push_back(solf[3*ord*i+E]);
        }
        if(ord==2)
        {
            //sol[i]=solf[ord*i]+coeff(1)*solf[ord*i+1]*legendre(1,i,x);
            V[RHO].push_back(solf[3*ord*i]+coeff(1)*solf[3*ord*i+1]*legendre(1,i,x));
            V[U].push_back(solf[3*ord*i+2]+coeff(1)*solf[3*ord*i+3]*legendre(1,i,x));
            V[E].push_back(solf[3*ord*i+4]+coeff(1)*solf[3*ord*i+5]*legendre(1,i,x));
        }
        if(ord==3)
        {
            sol[i]=solf[3*ord*i]+coeff(1)*solf[3*ord*i+1]*legendre(1,i,x)+coeff(2)*solf[3*ord*i+2]*legendre(2,i,x);
            V[RHO].push_back(solf[3*ord*i]+coeff(1)*solf[3*ord*i+1]*legendre(1,i,x)+coeff(2)*solf[3*ord*i+2]*legendre(2,i,x));
            V[U].push_back(solf[3*ord*i+3]+coeff(1)*solf[3*ord*i+4]*legendre(1,i,x)+coeff(2)*solf[3*ord*i+5]*legendre(2,i,x));
            V[E].push_back(solf[3*ord*i+6]+coeff(1)*solf[3*ord*i+7]*legendre(1,i,x)+coeff(2)*solf[ord*i+8]*legendre(2,i,x));
        }
    }
    for(int i=0;i<std::floor(L/dx);i++)
    {
        V[U].at(i)=V[U].at(i)/V[RHO].at(i);
    }
    for(int i=0;i<std::floor(L/dx);i++)
    {
        V[P].push_back((gama-1.0)*(V[E].at(i)-0.5*V[RHO].at(i)*V[U].at(i)*V[U].at(i)));
    }
    return V;
}








std::map<int,std::vector<double>> FU(std::vector<double> sol,double x)
{
    std::map<int,std::vector<double>> F;
    std::map<int,std::vector<double>> V;

    V=transformp(sol,x);

    for(int i=0;i<std::floor(L/dx);i++)
    {
        F[RHO].at(i)=V[RHO].at(i)*V[U].at(i);
        F[U].at(i)=V[RHO].at(i)*V[U].at(i)*V[U].at(i)+V[P].at(i);
        F[E].at(i)=V[U].at(i)*(V[E].at(i)+V[P].at(i));
    }

    return V;
}


std::map<int,std::vector<double>> transform(std::vector<double> solf)
{
    std::vector<double> sol(std::floor(L/dx),0);
    std::map<int,std::vector<double>> V;

    for(int i=0;i<std::floor(L/dx);i++)
    {
        V[RHO].push_back(solf[3*i]);

        V[U].push_back(solf[3*i+1]);

        V[E].push_back(solf[3*i+2]);
    }

    for(int i=0;i<std::floor(L/dx);i++)
    {

       V[P].push_back((gama-1)*(V[E].at(i)-0.5*V[U].at(i)*V[U].at(i)/V[RHO].at(i)));

    }


    return V;


}

std::map<int,std::vector<double>> FFU(std::vector<double> sol)
{
    std::map<int,std::vector<double>> F;
    std::map<int,std::vector<double>> V;

    V=transform(sol);

    for(int i=0;i<std::floor(L/dx);i++)
    {
        F[RHO].at(i)=V[U].at(i);
        F[U].at(i)=V[RHO].at(i)*V[U].at(i)*V[U].at(i)/V[RHO].at(i)+V[P].at(i);
        F[E].at(i)=V[U].at(i)*(V[E].at(i)+V[P].at(i))/V[RHO].at(i);
    }

    return V;
}