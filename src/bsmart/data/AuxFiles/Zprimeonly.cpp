#include "../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include "../include/micromegas_f.h"
#include "lib/pmodel.h"
#include <string>
#include <sys/utsname.h>
#include "../CalcHEP_src/include/rootDir.h" 
#include <unistd.h>
using namespace std;

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* MAIN PROGRAM (by M. Goodsell) use from micrOmegas v5.3.41        		    */
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

		FILE *omega = fopen("omg.out","w");
		
		
 

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
     //printf("Looking at particle %d\n",i);
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
     
     // Compute single production for even spin one particles
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

	   double zpeeddata[9];

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
	       
               fprintf(omega,"XSECTION  %E   2212  2212  2  %d \n",2*Pcm, pNum(plist[i])); 
/*pb*/         fprintf(omega," 0  0  0  0  0  0 %E micrOMEGAs\n\n", dcs);
               nChan++;
             }
             if(Pcm==PcmMax) break; else Pcm=PcmMax;
     
	  }
	 }
       }



	


     
      

       fclose(omega);

  	return 0;
}

