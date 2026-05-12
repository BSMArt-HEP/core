#include "../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include "lib/pmodel.h"
#include <string>

using namespace std;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* MAIN PROGRAM (by M. Goodsell) use from micrOmegas v5.3.41        		    */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int main(int argc, char** argv)
{  
   //nPROCSS=0; /* to switch off multiprocessor calculations */
		int err, i;
	   	char lspname[10], nlspname[10];
		double Omega=-1, Xf=-1;
		double w;
		double cut = 0.01;		// cut-off for channel output								
		int fast = 1;			/* 0 = best accuracy, 1 = "fast option" accuracy ~1% 	     */
 		double Beps = 1.E-5;  		/* Criteqrium for including co-annihilations (1 = no coann.) */
 		VZdecay=0; VWdecay=0; cleanDecayTable();
		ForceUG=1; 		
		err = sortOddParticles(lspname);	
		printMasses(stdout,1);


		
		
		if(CDM1 && CDM2) 
		  {
  
		    Omega= darkOmega2(fast,Beps,&err);
		    printf("Omega h^2=%.2E\n",Omega);
		    printf("Omega_1h^2=%.2E\n", Omega*(1-fracCDM2));
		    printf("Omega_2h^2=%.2E\n", Omega*fracCDM2);
		  }
		else
		  {
		    Omega = darkOmega(&Xf,fast,Beps,&err);
		    printf("Xf=%.2e Omega h^2=%.2E\n",Xf,Omega);
		    printChannels(Xf,cut,Beps,1,stdout);
		  }
		
//   			printChannels(Xf,cut,Beps,1,stdout);
		//	printf("\n");
		
		FILE *omega = fopen("omg.out","w");
			//	fprintf(omega,"%i %6.6lf # relic density \n",1,Omega);
		fprintf(omega,"%i %10.6e # Total relic density Omega h^2 \n",1,Omega);

		
		fprintf(omega,"%i %d # CDM1 \n",2,pNum(CDM1));
		fprintf(omega,"%i %10.6e # CDM1 Mass (GeV) \n",3,pMass(CDM1));
		

			if(CDM2)
			  {
			    fprintf(omega,"%i %10.6e # CDM1 relic density\n",4,Omega*(1-fracCDM2));
			    fprintf(omega,"%i %d # CDM2 \n",5,pNum(CDM2));
			    fprintf(omega,"%i %10.6e # CDM2 Mass (GeV)\n",6,pMass(CDM2));
			    fprintf(omega,"%i %10.6e # CDM2 relic density\n",7,Omega*fracCDM2);
			  }
			else
			  {
			    // Only print channels if there is one CDM candidate
			    w = 1.;
			    i = 0;
			    while (w>cut) 
			      {
				fprintf(omega,"%i %6.6lf # %s %s -> %s %s\n",100+i,omegaCh[i].weight,omegaCh[i].prtcl[0],omegaCh[i].prtcl[1],omegaCh[i].prtcl[2],omegaCh[i].prtcl[3]);
				i++;
				w = omegaCh[i].weight;
			      }
			  }
			




{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;  
  

printf("\n==== Calculation of CDM-nucleon amplitudes  =====\n");   
printf("         TREE LEVEL\n");

    nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
    printf("CDM-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 


printf("         BOX DIAGRAMS\n");  

   
    nucleonAmplitudes(CDM1,  pA0,pA5,nA0,nA5);
    printf("CDM-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0],pA5[0]);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0],nA5[0]); 

  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("CDM-nucleon cross sections[pb]:\n");
    printf(" proton  SI %.3E  SD %.3E\n",SCcoeff*pA0[0]*pA0[0],3*SCcoeff*pA5[0]*pA5[0]);
    printf(" neutron SI %.3E  SD %.3E\n",SCcoeff*nA0[0]*nA0[0],3*SCcoeff*nA5[0]*nA5[0]);
    
    
    fprintf(omega,"201 %10.6e #\n",SCcoeff*pA0[0]*pA0[0]);
    fprintf(omega,"202 %10.6e #\n",3*SCcoeff*pA5[0]*pA5[0]);
    fprintf(omega,"203 %10.6e #\n",SCcoeff*nA0[0]*nA0[0]);
    fprintf(omega,"204 %10.6e # \n",3*SCcoeff*nA5[0]*nA5[0]);
    
    
}


#ifdef NEVENTS
{
  double dNdE[300];
  double nEvents;
  double nEventsCut;

printf("\n======== Direct Detection ========\n");    




  nEvents=nucleusRecoil(Maxwell,73,Z_Ge,J_Ge73,SxxGe73,dNdE);
  printf("73Ge: Total number of events=%.2E /day/kg\n",nEvents);
  nEventsCut=cutRecoilResult(dNdE,10,50);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",nEventsCut);                                   ;
  fprintf(omega,"301 %10.6e # nEvents in Germanium/day/kg\n",nEvents);
                                                                                                         
  nEvents=nucleusRecoil(Maxwell,131,Z_Xe,J_Xe131,SxxXe131,dNdE);
  printf("131Xe: Total number of events=%.2E /day/kg\n",nEvents);
  //  nEventsCut=cutRecoilResult(dNdE,10,50);
  nEventsCut=cutRecoilResult(dNdE,4.9,40.9);
  printf("Number of events in 4.9 - 40.9 KeV region=%.2E /day/kg\n",nEventsCut);
  
  //  fprintf(omega,"302 %10.6e #\n",nEvents);
  fprintf(omega,"302 %10.6e # nEvents in Xenon in 4.9 - 40.9 KeV region /day/kg\n",nEventsCut);
  
  nEvents=nucleusRecoil(Maxwell,23,Z_Na,J_Na23,SxxNa23,dNdE);
  printf("23Na: Total number of events=%.2E /day/kg\n",nEvents);
  nEventsCut=cutRecoilResult(dNdE,10,50);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",nEventsCut);  
  fprintf(omega,"303 %10.6e # nEvents in Sodium/day/kg\n",nEvents);
  
  nEvents=nucleusRecoil(Maxwell,127,Z_I,J_I127,SxxI127,dNdE);
  printf("I127: Total number of events=%.2E /day/kg\n",nEvents);
  nEventsCut=cutRecoilResult(dNdE,10,50);
  printf("Number of events in 10 - 50 KeV region=%.2E /day/kg\n",nEventsCut);  
  fprintf(omega,"304 %10.6e # nEvents in Iodine/day/kg\n",nEvents);


}
#endif




 
{
  


char* expName; 
 printf("\n===== Direct detection exclusion:======\n");
  double pval=DD_pval(AllDDexp, Maxwell, &expName);
         if(pval<0.1 )  printf("Excluded by %s  %.1f%%\n", expName, 100*(1-pval)); 
  else printf("Not excluded by DD experiments  at 90%% level \n");

   fprintf(omega,"401 %8.6f # Direct detection pval \n",pval);
	 

 

}


/* LHC cross-sections, if required */

 
if( (argc > 1)  && (strcmp(argv[1],"--xsections")==0))
{
  //printf("%s\n",argv[1]);
	printf("Computing hadron collider cross-sections\n");
	int SMP[17]={1,2,3,4,5,6, 11,12,13,14,15,16, 21,22,23,24,25};
   int j;
   int nChan=0;
   double PcmMax=6500; // LHC pcm
   double PcmMin=4000;
   double csMinFb=0.001; 
   //if((LHC8|LHC13) == 0) 
   //{ printf("SMODELS: The third parameter has to be either  LHC8 or  LHC13 or  LHC8+LHC13 \n");      
    // return 1;
   //}  
   //if((Run & LHC8) == 0)  PcmMin=6500;
   //if((Run & LHC13) == 0) PcmMax=4000;
   
   //FILE*f=fopen(fileName,"w");
   int np=0;
   int nf=5;
   int wrt=1;
   //char**plist=NULL;
	char* plist[nModelParticles*2];
	for(i=0;i<nModelParticles;i++) if(pMass(ModelPrtcls[i].name) <PcmMax)
   { 
     for(j=0;j<17;j++) if(abs(ModelPrtcls[i].NPDG)==SMP[j]) break; 
     if(j==17 )
     	{ 
        	np++; 
        	//plist=realloc(plist,np*sizeof(char*));
        	plist[np-1]=ModelPrtcls[i].name;
        	if(strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname))
        	{ np++;
          	//plist=realloc(plist,np*sizeof(char*));
          	plist[np-1]=ModelPrtcls[i].aname;
        	}    
        
   		}

   }

int nZprime=0;
	for(i=0; i< np; i++)
     if(plist[i][0] !='~')
       {
           double dcs;
          double Qf=pMass(plist[i]);

	  int spinp1;

       qNumbers(plist[i], &spinp1, NULL,NULL);

       if (spinp1 == 2)
	 {
	   nZprime++;
	   // Write info for ZPEED

	   double coeff[10];

	   double zpeeddata[9]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

	    char *zpexlorerlabels[] ={"PDG","Zp-down-quark vector coupling","Zp-down-quark axial coupling","Zp-up-quark vector coupling","Zp-up-quark axial coupling","Zp-electron vector coupling","Zp-electron axial coupling","Zp-muon vector coupling","Zp-muon axial coupling"};
	   zpeeddata[0]=pNum(plist[i]);
	   int ferm[4]={1,2,13,15};
	   for(j=0;j<4;j++)
   {  char *fe=pdg2name(ferm[j]); if(!fe) continue;
      char *Fe=antiParticle(fe);  if(!Fe) continue; 
      //double mf=pMass(fe); if(mf==0) continue;
      lVert *Zpff=getVertex(Fe,fe,plist[i],NULL); 
      double c[2]={0,0};
      if(Zpff)
      {  getNumCoeff(Zpff,coeff); 
         for(int k=0;k<Zpff->nTerms;k++)
         {

	   //E1      |e1      |Zp      |        |1/2                                     |G(m3)*(1-G5)*v12542+G(m3)*(1+G5)*v12543 
	   printf("%s %12.4E\n",Zpff->SymbVert[k],coeff[k]);
	   if(strcmp(Zpff->SymbVert[k],"G5*G(m3)")==0)    c[1]=-coeff[k]; else   // ZPEED has other way round 
           if(strcmp(Zpff->SymbVert[k],"G(m3)")==0) c[0]=coeff[k];
         }
	 // Z' explorer requires Mz' then sequence of L/R couplings
	 // 

	 // NB in their example model they have
	 // g_V *1/2 * g /cosW, g_A * 1/2*g/cosW in the lagrangian, not the same thing.
	 //
	 // In our CHep model we have 1/2 *(G(m3)*(1-G5)*v12542+G(m3)*(1+G5)*v12543 )
	 // But these get simplified to c[0]* G5*G(m3) + c[1]*G(m3)

	 // Need gL 1/2 G(m3)*(1-G5) + gR 1/2 G(m3)*(1+G5) = 1/2 (gL+gR) G(m3) + 1/2(gR - gL) G(m3) G5
	 // = 1/2 (gL+gR) G(m3) + 1/2(gL - gR) G5 G(m3)

	 // gL = c[0] + c[1], gR = c[0] - c[1]

	 
	 zpeeddata[2*j+1]=c[0];
	 zpeeddata[2*j+2]=c[1];
	 
      }
	 
      }

	   fprintf(omega,"BLOCK ZPEEDIN #\n");
	   fprintf(omega," %d  %d  # %s\n",10*(nZprime-1)+1,int(zpeeddata[0]),zpexlorerlabels[0]);
	   for(int k=1;k<9;k++)
	     {
	       fprintf(omega," %d  %10.4E  # %s\n",10*(nZprime-1)+k+1,zpeeddata[k],zpexlorerlabels[k]);
	     }
	   
	   


	   // Write single production cross-sections for Z' explorer


	   
          for(double Pcm=PcmMin; ; ) 
          {
	    printf("Computing single production of %s\n",plist[i]);
             dcs=hCollider(Pcm,1,nf,Qf,Qf,plist[i],NULL,0,wrt);
             if(dcs>csMinFb*0.001)
             {
	       
               fprintf(omega,"XSECTION  %E   2212  2212  1  %d \n",2*Pcm, pNum(plist[i])); 
/*pb*/         fprintf(omega," 0  0  0  0  0  0 %E micrOMEGAs\n\n", dcs);
               nChan++;
             }
             if(Pcm==PcmMax) break; else Pcm=PcmMax;
     
	  }
	 }
       }


	// End single production cross-sections
	

	for(i=0;i<np;i++) for(j=i;j<np;j++) if(pMass(plist[i])+pMass(plist[j])<PcmMax)
//    if(plist[i][0]=='~' && plist[j][0]=='~')
    {  
		if(plist[i][0]=='~' && plist[j][0]!='~') continue;
		if(plist[i][0]!='~' && plist[j][0]=='~') continue;
		
		int q31,q32,q3,c1,c2;

       qNumbers(plist[i], NULL, &q31,&c1);
       qNumbers(plist[j], NULL, &q32,&c2);
       
       q3=q31+q32;

       if(q3<0) { q3*=-1; if(abs(c1)==3) c1*=-1; if(abs(c2)==3)  c2*=-1;}
       if(c1>c2){ int c=c1; c1=c2;c2=c;}
       if(c1==8) c1=1;
       if(c2==8) c2=1;

       int ok=0;

       switch(q3)
       {  case 0:  if( (c1==-3 && c2==3 ) || (c1==1  && c2==1) ) ok=1;
                   break;
          case 1:  if( (c1==3  && c2==3 ) || (c1==-3 && c2==1) ) ok=1;
                   break;
          case 2:  if( (c1==-3 && c2==-3) || (c1==1  && c2==3) ) ok=1;
                   break;
          case 3:  if( (c1==-3 && c2==-3) || (c1==1  && c2==1) ) ok=1;
                   break;
          case 4:  if( (c1==3  && c2==3) ) ok=1;
                   break;
       }
   
       if(ok==0) continue;
        
       {  double dcs;
          double Qf=0.5*(pMass(plist[i])+pMass(plist[j]));
          for(double Pcm=PcmMin; ; ) 
          { 
             dcs=hCollider(Pcm,1,nf,Qf,Qf,plist[i],plist[j],0,wrt);
             if(dcs>csMinFb*0.001)
             {
               fprintf(omega,"XSECTION  %E   2212  2212  2  %d  %d\n",2*Pcm, pNum(plist[i]),pNum(plist[j])); 
/*pb*/         fprintf(omega," 0  0  0  0  0  0 %E micrOMEGAs\n\n", dcs);
               nChan++;
             }
             if(Pcm==PcmMax) break; else Pcm=PcmMax;
          }  
       }
    }    
 }

       fclose(omega);

  	return 0;
}

