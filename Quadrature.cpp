#include "include/Quadrature.h"




std::vector<double> Quadrature0()
{
    std::vector<double> Q(3*ord*std::floor(L/dx),0);
    for(int j=0;j<std::floor(L/dx);j++)
    {
        Q[3*ord*j+RHO*ord]=(0.5)*((5.0/9.0)*u0(RHO,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*u0(RHO,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*u0(RHO,dx*j));


        Q[3*ord*j+U*ord]=(0.5)*((5.0/9.0)*u0(U,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*u0(U,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*u0(U,dx*j));


        Q[3*ord*j+E*ord]=(0.5)*((5.0/9.0)*u0(E,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*u0(E,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*u0(E,dx*j));



        if(ord>1)
        {
            Q[3*ord*j+RHO*ord+1]=(0.5/dx)*((5.0/9.0)*u0(RHO,dx*(j-0.5*sqrt(3.0/5.0)))*legendre(1,j,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*u0(RHO,dx*(j+0.5*sqrt(3.0/5.0)))*legendre(1,j,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*u0(RHO,dx*j)*legendre(1,j,dx*j));

            Q[3*ord*j+U*ord+1]=(0.5/dx)*((5.0/9.0)*u0(U,dx*(j-0.5*sqrt(3.0/5.0)))*legendre(1,j,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*u0(U,dx*(j+0.5*sqrt(3.0/5.0)))*legendre(1,j,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*u0(U,dx*j)*legendre(1,j,dx*j));

            Q[3*ord*j+E*ord+1]=(0.5/dx)*((5.0/9.0)*u0(E,dx*(j-0.5*sqrt(3.0/5.0)))*legendre(1,j,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*u0(E,dx*(j+0.5*sqrt(3.0/5.0)))*legendre(1,j,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*u0(E,dx*j)*legendre(1,j,dx*j));
        }
        if(ord==3)
        {
            Q[3*ord*j+RHO*ord+2]=(0.5/(dx*dx))*((5.0/9.0)*u0(RHO,dx*(j-0.5*sqrt(3.0/5.0)))*legendre(2,j,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*u0(RHO,dx*(j+0.5*sqrt(3.0/5.0)))*legendre(2,j,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*u0(RHO,dx*j)*legendre(2,j,dx*j));

            Q[3*ord*j+U*ord+2]=(0.5/(dx*dx))*((5.0/9.0)*u0(U,dx*(j-0.5*sqrt(3.0/5.0)))*legendre(2,j,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*u0(U,dx*(j+0.5*sqrt(3.0/5.0)))*legendre(2,j,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*u0(U,dx*j)*legendre(2,j,dx*j));

            Q[3*ord*j+E*ord+2]=(0.5/(dx*dx))*((5.0/9.0)*u0(E,dx*(j-0.5*sqrt(3.0/5.0)))*legendre(2,j,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*u0(E,dx*(j+0.5*sqrt(3.0/5.0)))*legendre(2,j,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*u0(E,dx*j)*legendre(2,j,dx*j));
        }
        
    }
    return Q;
}


std::vector<double> Quadrature(std::vector<double>& sol)
{
    std::vector<double> Q(3*ord*std::floor(L/dx),0);
    for(int j=0;j<std::floor(L/dx);j++)
    {
        Q[3*ord*j+RHO*ord]=0.0;


        Q[3*ord*j+U*ord]=0.0;


        Q[3*ord*j+E*ord]=0.0;



        if(ord>1)
        {
            Q[3*ord*j+RHO*ord+1]=(0.5)*((5.0/9.0)*FU(sol,dx*(j-0.5*sqrt(3.0/5.0)))[RHO][j]+(5.0/9.0)*FU(sol,dx*(j+0.5*sqrt(3.0/5.0)))[RHO][j]+(8.0/9.0)*FU(sol,dx*j)[RHO][j]);


            Q[3*ord*j+U*ord+1]=(0.5)*((5.0/9.0)*FU(sol,dx*(j-0.5*sqrt(3.0/5.0)))[U][j]+(5.0/9.0)*FU(sol,dx*(j+0.5*sqrt(3.0/5.0)))[U][j]+(8.0/9.0)*FU(sol,dx*j)[U][j]);


            Q[3*ord*j+E*ord+1]=(0.5)*((5.0/9.0)*FU(sol,dx*(j-0.5*sqrt(3.0/5.0)))[E][j]+(5.0/9.0)*FU(sol,dx*(j+0.5*sqrt(3.0/5.0)))[E][j]+(8.0/9.0)*FU(sol,dx*j)[E][j]);

        }
        if(ord==3)
        {
            Q[3*ord*j+RHO*ord+2]=(0.5/(dx*dx))*((5.0/9.0)*FU(sol,dx*(j-0.5*sqrt(3.0/5.0)))[RHO][j]*legendre(1,j,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*FU(sol,dx*(j+0.5*sqrt(3.0/5.0)))[RHO][j]*legendre(1,j,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*FU(sol,dx*j)[RHO][j]*legendre(1,j,dx*j));

            Q[3*ord*j+U*ord+2]=(0.5/(dx*dx))*((5.0/9.0)*FU(sol,dx*(j-0.5*sqrt(3.0/5.0)))[U][j]*legendre(1,j,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*FU(sol,dx*(j+0.5*sqrt(3.0/5.0)))[U][j]*legendre(1,j,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*FU(sol,dx*j)[U][j]*legendre(1,j,dx*j));

            Q[3*ord*j+E*ord+2]=(0.5/(dx*dx))*((5.0/9.0)*FU(sol,dx*(j-0.5*sqrt(3.0/5.0)))[E][j]*legendre(1,j,dx*(j-0.5*sqrt(3.0/5.0)))+(5.0/9.0)*FU(sol,dx*(j+0.5*sqrt(3.0/5.0)))[E][j]*legendre(1,j,dx*(j+0.5*sqrt(3.0/5.0)))+(8.0/9.0)*FU(sol,dx*j)[E][j]*legendre(1,j,dx*j));
        }
        
    }
    return Q;
}