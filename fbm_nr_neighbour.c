#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"ran2.c"
#define L 100000  ///system size

long sd=-937176386;

int main()
{
  int i,j,n,b,bundle=1,NRUN=1,t,count2,n_l[L],n_r[L],Nr,breaking,time;
  int *stat,count,jl,jr,flag5;
  float r,max_thr ,beta,f_max,sumf,sum4,*thr,*force,del,load,total_load,energy;
  float min,sum,w,sum1,sum2;
  float avalanche_series[L], energy_series[L]/*,load_series[L]*/;
  float load_series[L];
  FILE*fp;
  fp=fopen("output.txt","w");  /// here the avalanche series and energy series are saved in the output file.

  Nr=1;
  beta=1.1;

  stat=(int *)malloc(L*sizeof(int));
  thr=(float *)malloc(L*sizeof(float));
  force=(float *)malloc(L*sizeof(float));

for(b=0;b<bundle;b++)
{

  for(i=0;i<L;i++){avalanche_series[i]=0; energy_series[i]=0;load_series[i]=0;}

   load=0.;time=0;
   for(n=0;n<NRUN;n++)
   {
        load=0.;
        for(i=0;i<L;i++)
        {

         stat[i]=1;
        thr[i]=ran2(&sd);   ///uniform distribution

         //r=2*ran2(&sd)-1;
        // thr[i]=pow(10,-1*beta*r);
         force[i]=0;
         n_l[i]=(i-1+L)%L;
         n_r[i]=(i+1)%L;
         if(max_thr<thr[i])
           {
              max_thr=thr[i];
           }
          // printf("\n%f",thr[i]);
        }
        count=1; count2=1;

	for(w=0.;w<max_thr;)
  	{
        sum4=0;
        count=0;
        for(i=0;i<L;i++)
        {
          count+=stat[i];
        }

        if(count>0)
        {

          del=max_thr;
          for(i=0;i<L;i++)
          {
             if(stat[i]==1)
             {
               if((thr[i]-force[i])<del) {del=thr[i]-force[i];}
             }
          }
          count2=0;
	  for(i=0;i<L;i++)
    	  {
             force[i]=stat[i]*(force[i]+del);
             sum4+=stat[i]*force[i];
             if(stat[i]*force[i]>=thr[i]) {count2++;}
       	  }

        }
        else {w=max_thr+1; /*escape route*/}


        time++;

        breaking=0;energy=0;
        while(count>0 && count2>0)
    	{
          count=0;
          for(i=0;i<L;i++)
          {
            if(force[i]>=thr[i]) {stat[i]=0; n_r[n_l[i]]=n_r[i]; n_l[n_r[i]]=n_l[i]; breaking++; energy+=(thr[i]*thr[i])*0.5;}
            count+=stat[i];
          }
         //fprintf(fp,"\n%d",breaking);

          if(count>0)
          {
            for(i=0;i<L;i++)
            {
               if(stat[i]==0 && force[i]>0)
               {
                  jl=i; jr=i;
                  j=0; flag5=0;
                  while(j<Nr && flag5<L)
                  {
                     if(stat[n_l[jl]]==1) {force[n_l[jl]]+=force[i]/(2*Nr); j++;}
                     jl=n_l[jl];
                     flag5++;
                  }

                  j=0; flag5=0;
                  while(j<Nr && flag5<L)
                  {
                     if(stat[n_r[jr]]==1) {force[n_r[jr]]+=force[i]/(2*Nr); j++;}
                     jr=n_r[jr];
                     flag5++;
                  }

                  force[i]=0;
               }
            }

          }


          count2=0;
         total_load=0;
          for(i=0;i<L;i++)
          {
            if(stat[i]==1 && count>0)
                {total_load+=force[i]; if(force[i]>=thr[i]) {count2++;}}

          }


        }   //////end of while loop/////////////////
      //  fprintf(fp,"\n%d",breaking);
//end of time

        avalanche_series[time]=breaking; energy_series[time]=energy; load_series[time]=total_load/L;

        count=0; total_load=0;
        for(i=0;i<L;i++) {count+=stat[i]; total_load+=force[i];}
        if(count>0) {load+=del*count/L; w=total_load/L;}


         }
//end of force loop(for loop)

  }
//end of realisation loop



  for(i=1;i<time;i++)
  {
       fprintf(fp,"%d %f %f\n",i,avalanche_series[i],energy_series[i]);
   // printf("%d  %f  %f\n",i,avalanche_series[i],load_series[i]);

  }
}


}
