#include "include/Limiter.h"





std::vector<double> left(std::vector<double>& sol)
{
    std::vector<double> l(3*std::floor(L/dx),0.0);
    
    if(ord==2)
    {
        for (int i = 0; i <std::floor(L/dx); i++)
        {
            l[3*i]=6*sol[3*ord*i+1];

            l[3*i+1]=6*sol[3*ord*i+3];

            l[3*i+3]=6*sol[3*ord*i+5];
        }
            
    }
            
    if(ord==3)
    {
        for (int i = 0; i <std::floor(L/dx); i++)
        {
            l[3*i]=6*sol[3*ord*i+1]+30*sol[3*ord*i+2];

            l[3*i+1]=6*sol[3*ord*i+4]+30*sol[3*ord*i+5];

            l[3*i+3]=6*sol[3*ord*i+7]+30*sol[3*ord*i+8];
        }
            
    }
    
    return l;
}


std::vector<double> wright(std::vector<double>& sol)
{
    std::vector<double> l(3*std::floor(L/dx),0.0);
   
    if(ord==2)
    {
         for (int i = 0; i <std::floor(L/dx); i++)
        {
            l[3*i]=6*sol[3*ord*i+1];

            l[3*i+1]=6*sol[3*ord*i+3];

            l[3*i+3]=6*sol[3*ord*i+5];
        }
        

    }
       
    

    if(ord==3)
    {
        for (int i = 0; i <std::floor(L/dx); i++)
        {
            l[3*i]=6*sol[3*ord*i+1]-30*sol[3*ord*i+2];

            l[3*i+1]=6*sol[3*ord*i+4]-30*sol[3*ord*i+5];

            l[3*i+3]=6*sol[3*ord*i+7]-30*sol[3*ord*i+8];
        }     
 
    }
        

    return l;
}


std::vector<double> Tleft(std::vector<double>& sol)
{
    

    std::vector<double> l1(3*std::floor(L/dx),0.0);
    
    l1=left(sol);
    
    std::vector<double> l(3*std::floor(L/dx),0.0);
    
    if(ord==1)
    {
        for (int i = 0; i <std::floor(L/dx); i++)
        {
            l[3*i]=sol[3*ord*i];

            l[3*i+1]=sol[3*ord*i+1];

            l[3*i+3]=sol[3*ord*i+2];
        }
    }    
    if(ord==2)
    {
        for (int i = 0; i <std::floor(L/dx); i++)
        {
            l[3*i]=sol[3*ord*i];

            l[3*i+1]=sol[3*ord*i+2];

            l[3*i+3]=sol[3*ord*i+4];
        }
    
    }
    if(ord==3)
    {
        for (int i = 0; i <std::floor(L/dx); i++)
        {
            l[3*i]=sol[3*ord*i];

            l[3*i+1]=sol[3*ord*i+3];

            l[3*i+3]=sol[3*ord*i+6];
        }

    }

    l=l+l1;

    return l;
    
}

std::vector<double> Twright(std::vector<double>& sol)
{
    

    std::vector<double> l1(3*std::floor(L/dx),0.0);

    l1=wright(sol);

    std::vector<double> l(3*std::floor(L/dx),0.0);

    if(ord==1)
    {
        for (int i = 0; i <std::floor(L/dx); i++)
        {
            l[3*i]=sol[3*ord*i];

            l[3*i+1]=sol[3*ord*i+1];

            l[3*i+3]=sol[3*ord*i+2];
        }

        l=l-l1;

        
    }
        
        
    if(ord==2)
    {
        for (int i = 0; i <std::floor(L/dx); i++)
        {
            l[3*i]=sol[3*ord*i];

            l[3*i+1]=sol[3*ord*i+2];

            l[3*i+3]=sol[3*ord*i+4];
        }
    
        l=l-l1;

        
    }
    if(ord==3)
    {
        for (int i = 0; i <std::floor(L/dx); i++)
        {
            l[3*i]=sol[3*ord*i];

            l[3*i+1]=sol[3*ord*i+3];

            l[3*i+3]=sol[3*ord*i+6];
        }

        l=l-l1;

        
    }

    return l;
    
}

