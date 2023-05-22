#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran2.c"
#define L 30000     ///system size

long sd =-971031158;
double weibull();
double shape_parameter;


int main()
{
    int i,j, n;
    int *store,intact,broken,bundle=10000,store2,*avl_size,sum,count=0;
    double *threshold,min=0,minimum=0,min_1,max=0,r,beta,energy,sq_sum;
    double loadFiber=0,loadBundle=0, load=0,*load_1;
    int intact1=0, broken1=0;
    int avl=0;


    threshold=(double *)malloc(L*sizeof(double));
    store=(int *)malloc(L*sizeof(int));
    avl_size=(int *)malloc(L*sizeof(int));
    load_1= (double *)malloc(L*sizeof(double));


    beta=4.5;
    shape_parameter=1.3;
    FILE*fp;
    fp=fopen("output.txt","w");

    for(j=0;j<bundle;j++)
    {
        min=0;max=0;loadFiber=0;loadBundle=0;load=0;intact1=0;broken1=0;
        for(i=0;i<L;i++){avl_size[i]=0;/*energy_2[n]=0;*/load_1[i]=0;}

    for(i=0;i<L;i++)
    {
       // r=2*ran2(&sd)-1;
         // threshold[i]=(double)i/L;      ///equispaced threshold
         // threshold[i]= ran2(&sd);       ///uniform distribution
        //threshold[i]=pow(10,-1*beta*r);   ///non-uniform distribution
         threshold[i]= weibull();          /// weibull distribution

    }


    max= threshold[0];
    minimum= threshold[0];
    for(i=0;i<L;i++)
    {
         if(max<threshold[i])
            max=threshold[i];
         if(minimum>threshold[i])
           minimum= threshold[i];


    }

    load=minimum;
    intact=L;
    avl=0;
    loadBundle=0;
    while(load<max)
    {
        avl++;
        intact1=0;
        loadFiber=load;
       // loadBundle=(intact*loadFiber);
        loadBundle= loadBundle+100;     ///discrete loading
        energy=0;
       while(2)
         {
            intact=0; broken=0;
            for(i=0;i<L;i++)
             {
               if(threshold[i]>loadFiber)
               {
                intact++;
               }
               else
               {
                 broken++;
                 energy += (threshold[i]*threshold[i])/2;
                 threshold[i]=0;

               }

             }

           if(intact1==intact)
             {
                 store[broken-broken1]++;
                 count++;
                 store2=broken-broken1;
                 sum+=store2;
                 broken1=broken;
                 break;
             }
             intact1=intact;

           loadFiber= (loadBundle/intact);
        }

          avl_size[avl]= store2;
          load_1[avl]=load;
          load=loadFiber;

      /*   min_1=1;
         for(i=0;i<L;i++)
         {
             min=(threshold[i]-loadFiber);
             if(min>0 && min<min_1)
             {
                 min_1=min;
             }
         }
         load=loadFiber+min_1;*/  /// quasistatic loading

    }
    for(i=1;i<avl;i++)
    {
         fprintf(fp,"%d    %d\n",i,avl_size[i]);  /// here avl_size is the avalanche size.
    }
    }
    //printf("\nsum of avalanche size:%d\n",sum);


   free(threshold);
   free(avl_size);
   free(store);

}

double weibull()
{
   double lambda,k,p,r;

   lambda=1; k=shape_parameter;


   r=ran2(&sd);

   p=pow(-1*pow(lambda,k)*log(1-r),1./k);

   return(p);
}
