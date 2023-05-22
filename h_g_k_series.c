#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"ran2.c"

long sd=-937176383;
main()
{
  int i,ii,avc,a,d,dd,L=10000,Max=920,h,h_prev,h_c,h_index,avc_h[Max],step_number,activity,avc_k[Max],k,avc_g[Max],all_inf,u_series[Max],shuffle,ni,nj,tenp,n,NRUN=1,max_shuffle;
  float sum,b,e,c,f,dist[L],series[Max],h_series[Max],delta,delta_prev,k_prev,k_series[Max],temp,cumulative,g,g_series[Max],lorenz[Max],beta,r,count;
  int d1,d2,d3,d4,d5;
  float c1,c2,c3,c4,c5;

  FILE *fp= fopen("input.txt","r"); /// In input file avalanche size should be there.


  for(i=0;i<L;i++) {dist[i]=0;}
  sum=0; avc=0; ii=0;

  for(i=0;i<Max;i++)
  {
    fscanf(fp,"%f\n",&d3);           /// here d3 is the avalanche size

    {series[i]=(float)d3; ii++;}
    sum+=d3;

  }

  Max=ii+1;

 for(n=0;n<NRUN;n++)
 {

    for(i=0;i<Max;i++) {avc_h[i]=0; avc_k[i]=0; h_series[i]=0; k_series[i]=0; g_series[i]=0; avc_g[i]=0; lorenz[i]=0;}

///////////////////////////////////////////////////////////////////////
  h_prev=0;

  for(step_number=0;step_number<Max;step_number++)
  {

/////////////////////////////////////////h-index measurement////////////////////////////////

        for(h=h_prev;h<Max;h++)
        {
           h_c=0;
           for(i=0;i<step_number;i++)
           {
              if(series[i]>=h) {h_c++;}
           }
           if(h_c<h) {h_index=h-1; h=Max;}

        }
        h_series[step_number]+=h_index; avc_h[step_number]++;
        h_prev=h_index;

///////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////k-index measurement////////////////////////////////////////////

   count=0;

   for(k=0;k<step_number;k++) {count+=series[k];}

   activity=1;
   while(activity>0)
   {
     activity=0;
     for(k=0;k<step_number;k++)
     {
        if(series[k]>series[k+1]) {temp=series[k+1];   series[k+1]=series[k]; series[k]=temp; activity++;}
     }

   }
   cumulative=0;
   delta_prev=0;
   k_prev=0;
   for(k=0;k<step_number;k++)
   {
     cumulative+=series[k];
     delta=(float)k/step_number+cumulative/count-1;


     if(delta*delta_prev<=0. && count>0. && k>0)
     {
      if(((k_prev)+(float)k/step_number)/2.<0.5) { k_series[step_number]+=0.5; avc_k[step_number]++;}   //////for plot of k-series, += sign needed in averaging
      else {k_series[step_number]+=((k_prev)+(float)k/step_number)/2.; avc_k[step_number]++;}      //////for plot of k-series, += sign needed in averaging
       k=step_number;

     }
     else {delta_prev=delta; k_prev=(float)k/step_number;}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   }


/////////////////////////////////////////////gini measurement///////////////////////////////////////////////////////////

   cumulative=0; g=0;
   for(k=0;k<step_number;k++)
   {
      cumulative+=series[k];

      g+=((float)k/step_number-cumulative/count)/(float)step_number;
   }

  if(g<0) {g=0;}

   g_series[step_number]+=2*g; avc_g[step_number]++;
}


///////////////////////////////////////////////////////////////////////////////////////////////
  }

  cumulative=0;
  for(i=0;i<=Max;i++)
  {
    cumulative+=series[i];
    lorenz[i]=(float)cumulative/sum;
  }



 for(i=0;i<Max;i++)
 {
  if(avc_h[i]>0 && avc_k[i]>0 && avc_g[i]>0)
  {
      printf("%d %f %f %f\n",i,h_series[i]/avc_h[i],k_series[i]/avc_k[i],g_series[i]/avc_g[i]);
  }

 }
//end of NRUN loop


}
