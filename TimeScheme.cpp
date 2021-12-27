#include "include/TimeScheme.h"



std::vector<double> Update(std::vector<double> sol, int flux)
{

    std::vector<double> U0;
    std::vector<double> U1;
    std::vector<double> U2;

    switch (r)
    {
    case 1:

        U0=sol;
       
        sol=U0+dt*LH(U0,flux);
       
        return sol;

    case 2:

        U0=sol;
        sol=sol+dt*LH(U0,flux);
        U1=sol;
        sol=0.5*U0+0.5*U1+0.5*dt*LH(U1,flux);

        return sol;

    case 3:

        U0=sol;
        sol=sol+dt*LH(U0,flux);
        U1=sol;
        sol=(3.0/4.0)*U0+0.25*U1+0.25*dt*LH(U1,flux);
        U2=sol;
        sol=(1.0/3.0)*U0+(2.0/3.0)*U2+(2.0/3.0)*dt*LH(U2,flux);

        return sol;


    
    default:
        std::cout<<"The order "<< r <<" Runge-Kutta order is not yet implemented "<<std::endl;
        break;
    }
}