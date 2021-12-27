#include "include/SpaceScheme.h"




std::map<int,double> Numeric_flux(std::vector<double>& sol,int j,int c, int flux)
{

}


std::vector<double> LH(std::vector<double> sol, int flux)
{
    std::vector<double> output(sol.size());
    std::map<int,double> wright;
    std::map<int,double> left;
    std::vector<double> Q;
    Q=Quadrature(sol);
    for(int j=0;j<std::floor(L/dx);j++)
    {
        wright[RHO]=Numeric_flux(sol,j,1,flux)[RHO];
        left[RHO]=Numeric_flux(sol,j,-1,flux)[RHO];

        wright[U]=Numeric_flux(sol,j,1,flux)[U];
        left[U]=Numeric_flux(sol,j,-1,flux)[U];

        wright[E]=Numeric_flux(sol,j,1,flux)[E];
        left[E]=Numeric_flux(sol,j,-1,flux)[E];


        if(ord==1)
        {
            output[3*ord*j]=-(1/dx)*(wright[RHO]-left[RHO]);  
            output[3*ord*j+1]=-(1/dx)*(wright[U]-left[U]);
            output[3*ord*j+3]=-(1/dx)*(wright[E]-left[E]);
        }
        
        
        if(ord==2)
        {
            output[3*ord*j]=-(1/dx)*(wright[RHO]-left[RHO]);  
            output[3*ord*j+2]=-(1/dx)*(wright[U]-left[U]);
            output[3*ord*j+4]=-(1/dx)*(wright[E]-left[E]);

            output[3*ord*j+1]=-(0.5/dx)*(wright[RHO]-left[RHO])+(1.0/(dx*dx))*Q[3*ord*j+1];  
            output[3*ord*j+3]=-(0.5/dx)*(wright[U]-left[U])+(1.0/(dx*dx))*Q[3*ord*j+3];
            output[3*ord*j+5]=-(0.5/dx)*(wright[E]-left[E])+(1.0/(dx*dx))*Q[3*ord*j+5];
            //output[ord*j+1]=-(0.5/dx)*(wright+left)+sol[ord*j]/dx;
        }
        if(ord==3)
        {
            output[3*ord*j]=-(1/dx)*(wright[RHO]-left[RHO]);  
            output[3*ord*j+3]=-(1/dx)*(wright[U]-left[U]);
            output[3*ord*j+6]=-(1/dx)*(wright[E]-left[E]);

            output[3*ord*j+1]=-(0.5/dx)*(wright[RHO]-left[RHO])+(1.0/(dx*dx))*Q[3*ord*j+1];  
            output[3*ord*j+4]=-(0.5/dx)*(wright[U]-left[U])+(1.0/(dx*dx))*Q[3*ord*j+4];
            output[3*ord*j+7]=-(0.5/dx)*(wright[E]-left[E])+(1.0/(dx*dx))*Q[3*ord*j+7];

            output[3*ord*j+2]=-(1/(6*dx))*(wright[RHO]-left[RHO])+(1.0/(dx*dx*dx))*Q[3*ord*j+2];  
            output[3*ord*j+5]=-(1/(6*dx))*(wright[U]-left[U])+(1.0/(dx*dx*dx))*Q[3*ord*j+5];
            output[3*ord*j+8]=-(1/(6*dx))*(wright[E]-left[E])+(1.0/(dx*dx*dx))*Q[3*ord*j+8];

        
            //output[ord*j+2]=-(1/(6*dx))*(wright-left)+2*sol[ord*j+1]/dx;
        }
    }
    return output;
}