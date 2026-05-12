#include "../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include "lib/pmodel.h"
#include <string>

using namespace std;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* MAIN PROGRAM (by M. Goodsell)			     		    */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

int main(int argc, char** argv)
{  
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
  
		    Omega= darkOmega2(fast,Beps);
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


       fclose(omega);

  	return 0;
}

