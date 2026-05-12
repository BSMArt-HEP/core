"""
Tool to run HiggsTools


Possible settings:

.. code-block:: text

    "Codes": {
        "HiggsTools": {
            'Run': "True", #
            'Signal Strengths': "False",  # Compute predicted signal strength ratios 
            'Neutral Higgs' : [25, 35, 45, 36, 46],   # list of IDs of neutral Higgs bosons to check
            'Charged Higgs' :  [37],  #  List of IDs of charged Higgs bosons to check
            'Invisible Particles' : [1000022, 1000023]  #  list of BSM invisible particle ids that the Higgs could decay to
            'HiggsBounds Dataset'  : "/path/to/dataset",
            'HiggsSignals Dataset': "/path/to/dataset"        
        }

    }


If 'Signal Strengths' is 'True' then the tool will compute ratios of production of the neutral scalars into different channels.

These channels are as defined in HiggsTools in higgstools-main/include/Higgs/predictions/Channels.hpp 

I select here production channels are `['H','HZ','vbfH','HW','eeHZ']`  (where 'H' means single Higgs production via ggH+ bbH + qqH combined)
and decay channels `['bb','ZZ','WW','gamgam','gg','mumu','tautau']`
More can easily be added

The outputs are stored in the HiggsTools BLOCK, and will typically look like the following (if Signal Strengths is True):

.. code-block:: text

    BLOCK HIGGSTOOLS 
    1  1  # HB result (0/1)
    2  8.06381191e-01  # HB maximum obs ratio
    11  1.53742211e+02  # HS chi2
    12  159  # HS number of observables
    13  6.02821197e-01  # BSMArt Inferred HS pval
    25  8.06381191e-01  # ID: 12045, ref: CMS-PAS-HIG-12-045: LHC8 [vbfH,HW,Htt,H,HZ]>[bb,tautau,WW,ZZ,gamgam]
    1000  9.87350705e-01  # mueff for  H -> 25 -> bb
    1001  9.31152726e-01  # mueff for  H -> 25 -> ZZ
    1002  1.02204352e+00  # mueff for  H -> 25 -> WW
    1003  1.06359389e+00  # mueff for  H -> 25 -> gamgam
    1004  5.47133098e-01  # mueff for  H -> 25 -> gg
    1005  1.06115939e+00  # mueff for  H -> 25 -> mumu
    1006  1.06280238e+00  # mueff for  H -> 25 -> tautau
    1100  9.87338291e-01  # mueff for  HZ -> 25 -> bb
    1101  9.31141018e-01  # mueff for  HZ -> 25 -> ZZ
    1102  1.02203067e+00  # mueff for  HZ -> 25 -> WW
    1103  1.06358052e+00  # mueff for  HZ -> 25 -> gamgam
    1104  5.47126218e-01  # mueff for  HZ -> 25 -> gg
    1105  1.06114605e+00  # mueff for  HZ -> 25 -> mumu
    1106  1.06278902e+00  # mueff for  HZ -> 25 -> tautau
    1200  9.87381710e-01  # mueff for  vbfH -> 25 -> bb
    1201  9.31181966e-01  # mueff for  vbfH -> 25 -> ZZ
    1202  1.02207562e+00  # mueff for  vbfH -> 25 -> WW
    1203  1.06362729e+00  # mueff for  vbfH -> 25 -> gamgam
    1204  5.47150279e-01  # mueff for  vbfH -> 25 -> gg
    1205  1.06119271e+00  # mueff for  vbfH -> 25 -> mumu
    1206  1.06283576e+00  # mueff for  vbfH -> 25 -> tautau
    1300  9.87396661e-01  # mueff for  HW -> 25 -> bb
    1301  9.31196066e-01  # mueff for  HW -> 25 -> ZZ
    1302  1.02209109e+00  # mueff for  HW -> 25 -> WW
    1303  1.06364339e+00  # mueff for  HW -> 25 -> gamgam
    1304  5.47158564e-01  # mueff for  HW -> 25 -> gg
    1305  1.06120878e+00  # mueff for  HW -> 25 -> mumu
    1306  1.06285185e+00  # mueff for  HW -> 25 -> tautau
    1400  9.87340225e-01  # mueff for  eeHZ -> 25 -> bb
    1401  9.31142843e-01  # mueff for  eeHZ -> 25 -> ZZ
    1402  1.02203267e+00  # mueff for  eeHZ -> 25 -> WW
    1403  1.06358260e+00  # mueff for  eeHZ -> 25 -> gamgam
    1404  5.47127290e-01  # mueff for  eeHZ -> 25 -> gg
    1405  1.06114813e+00  # mueff for  eeHZ -> 25 -> mumu
    1406  1.06279110e+00  # mueff for  eeHZ -> 25 -> tautau
    2000  0.00000000e+00  # mueff for  H -> 35 -> bb
    2001  0.00000000e+00  # mueff for  H -> 35 -> ZZ
    2002  0.00000000e+00  # mueff for  H -> 35 -> WW
    2003  0.00000000e+00  # mueff for  H -> 35 -> gamgam
    2004  0.00000000e+00  # mueff for  H -> 35 -> gg
    ...

"""

__meta__ = {
    "name": "HiggsTools",
    "requires": ["scipy", "numpy", "higgstools"],
    "settings": {
        "Invisible Particles": "List of invisible particle PIDs",
        "HiggsBounds Dataset": "Path to HB dataset",
        "HiggsSignals Dataset": "Path to HS dataset",
        "Signal Strengths": "True/False (default False): whether to compute signal strength ratios",
        "Neutral Higgs": "List of neutral Higgs PIDs",
        "Charged Higgs": "List of charged Higgs PIDs"
    }
}

import os
#import shutil
#import subprocess
from bsmart import debug
from bsmart.HEPRun import HepTool, DataPoint
from bsmart import zslha

import Higgs.predictions as HP
import Higgs.bounds as HB
import Higgs.signals as HS
from Higgs.tools.Input import predictionsFromDict #,readHB5SLHA

from collections import OrderedDict, defaultdict
from scipy.stats import chi2

from itertools import product, combinations_with_replacement, chain
import numpy as np

""" Pieces taken/adapted from HiggsTools readHB5SLHA 

Note that in zslha the br are stored with integer keys, while the HBFERMIONS and HBBOSONS are stored as strings.
OTOH the zslha.Value command converts the keys to strings for such blocks, as required
Hence we must provide the pdg numbers as integers


"""
def getCPFromEffC(
    p: str,
    data: dict,
    effCPrefix: str = "effc",
    evenCoupSuffix: str = "s",
    oddCoupSuffix: str = "p",
):
    """
    Use the fermionic effective couplings of a particle with ID p that are
    stored in the data dictionary to determine its CP quantum number.
    """
    coups = {k: v for k, v in data.items() if k.startswith(f"{effCPrefix}_{p}_")}
    evenCoups = [v for k, v in coups.items() if k.endswith(f"_{evenCoupSuffix}")]
    oddCoups = [v for k, v in coups.items() if k.endswith(f"_{oddCoupSuffix}")]
    return np.allclose(oddCoups, 0) - np.allclose(evenCoups, 0)

def readHB5ZSLHA(slha: zslha.SLHA,
    neutralPDGs: list,
    chargedPDGs: list,
    invisibleWidthThreshold: float = 1e-10,### Strange way of deciding invisible. 
    invisiblePDGs: list = [], ## NEW feature to circumvent the above
):
    
    
    #slha = zslha.read(spcfile)

    fermions = {
        **{v + 1: k for v, k in enumerate(["d", "u", "s", "c", "b", "t"])},
        **{11 + v: k for v, k in enumerate(["e", "nu", "mu", "nu", "tau", "nu"])},
    }
    bosons = {21 + v: k for v, k in enumerate(["g", "gam", "Z", "W"])}

    def parseCouplings():
        hbbosons="HIGGSBOUNDSINPUTHIGGSCOUPLINGSBOSONS"
        #if "HIGGSBOUNDSINPUTHIGGSCOUPLINGSBOSONS" in slha.blocks:
            
        if "HIGGSCOUPLINGSBOSONS" in slha.blocks:
            hbbosons="HIGGSCOUPLINGSBOSONS"

        hbfermions="HIGGSBOUNDSINPUTHIGGSCOUPLINGSFERMIONS"
        if "HIGGSCOUPLINGSFERMIONS" in slha.blocks:
            hbfermions="HIGGSCOUPLINGSFERMIONS"
        
        coups = dict()
        for pdg in neutralPDGs:
            coups.update(
                (f"effc_{pdg}_{k}{k}", slha.safeValue(hbbosons,[pdg,v,v])) for v, k in bosons.items() 
            )
            coups.update(
                (f"effc_{pdg}_{k}{k}_s", slha.safeValue2(hbfermions,[pdg,v,v])[0]) for v, k in fermions.items() if "nu" not in k
            )
            coups.update(
                (f"effc_{pdg}_{k}{k}_p", slha.safeValue2(hbfermions,[pdg,v,v])[1]) for v, k in fermions.items() if "nu" not in k
            )
            coups[f"CP_{pdg}"] = getCPFromEffC(str(pdg), coups)
        
        ### Need to handle the HHZ couplings
        
   
        for coup in slha.blocks[hbbosons]:
            if '23' in coup:        
                coupdata=list(map(int,coup.split(',')))
                if len(coupdata) == 3 and coupdata[-1] == 23:
                    higgspresent=sorted(coupdata[:-1])
                    if higgspresent[0] in neutralPDGs and higgspresent[1] in neutralPDGs:
                        coups["paircxn_%d_%d_LEP" %(higgspresent[0],higgspresent[1])] = float(slha.blocks[hbbosons][coup])**2
        ## Fill in the blanks? Don't know if this is necessary:
        for hi, hj in combinations_with_replacement(neutralPDGs, 2):
            tkey=f"paircxn_{hi}_{hj}_LEP"
            if tkey not in coups:
                coups[tkey] = 0.0
        """      
        for coup in coups:
            print('coup '+str(coup)+': '+str(coups[coup]))
        """     
        
        return coups

    def parseMasses():
        ## NB the 'Value' command converts the [pdg] to a string, as required
        return {f"m_{pdg}" : float(slha.Value("MASS",[pdg])) for pdg in chain(neutralPDGs, chargedPDGs)}
        #massVals = dict(slha.blocks["MASS"]["values"])
        #return {f"m_{pdg}": massVals[pdg] for pdg in chain(neutralPDGs, chargedPDGs)}

    def parseMassUncs():
        try:
            if 'DMASS' in slha.blocks:
                
                massUncVals = dict(slha.blocks["DMASS"])
                return { f"dm_{pdg}": massUncVals.get(pdg, 3.0) for pdg in chain(neutralPDGs, chargedPDGs)}
            else:
                return { f"dm_{pdg}": 3.0 for pdg in chain(neutralPDGs, chargedPDGs)}
            #massUncVals = dict(slha["BLOCK"]["DMASS"]["values"])
            
            #return {
            #    f"dm_{pdg}": massUncVals.get(pdg, 0.0)
            #    for pdg in chain(neutralPDGs, chargedPDGs)
            #}
        except KeyError:
            return dict()

    def parseDecays():
        widths = {
            #f"w_{pdg}": slha["DECAY"][str(pdg)]["info"][0]
            #for pdg in chain(neutralPDGs, chargedPDGs)
            f"w_{pdg}": slha.Value("WIDTH",pdg) for pdg in chain(neutralPDGs, chargedPDGs)
        }
        decays = {}
        for p in chain(neutralPDGs, chargedPDGs):
            
            SMparticles = {**fermions, **bosons}
            decayDat = [
                [*sorted(list(k)), val]  #slha.br[p][vals]
                for k,val in zip(slha.br[p].keys(), slha.br[p].values())
                if len(k) == 2
            ]
            
            def toDecay(d1, d2):
                d1 = abs(d1)
                d2 = abs(d2)
                if d1 in SMparticles and d2 in SMparticles:
                    if f"{SMparticles[d1]}{SMparticles[d2]}" in HP.Decay.__members__:
                        return f"br_{p}_{SMparticles[d1]}{SMparticles[d2]}"
                    elif f"{SMparticles[d2]}{SMparticles[d1]}" in HP.Decay.__members__:
                        return f"br_{p}_{SMparticles[d2]}{SMparticles[d1]}"
                    else:
                        return f"unknown_br_{p}_{SMparticles[d1]}_{SMparticles[d2]}"
                elif d1 in SMparticles:
                    return f"br_{p}_{SMparticles[d1]}_{d2}"
                elif d2 in SMparticles:
                    return f"br_{p}_{SMparticles[d2]}_{d1}"
                else:
                    return f"br_{p}_{d1}_{d2}"

            def isInvisible(d1, d2):
                if invisiblePDGs == []:
                    invisibleParticles = [
                        int(k)
                        for k, val in zip(slha.widths.keys(),slha.widths.values())
                        if (abs(int(k)) > 15 and val < invisibleWidthThreshold)
                    ] + [12, 14, 16]
                else:
                    invisibleParticles = invisiblePDGs
                return abs(d1) in invisibleParticles and abs(d2) in invisibleParticles
                
            decays.update(
                (toDecay(d1, d2), v)
                for d1, d2, v in decayDat
                if not isInvisible(d1, d2)
            )

            if p in neutralPDGs:
                decays[f"br_{p}_directInv"] = sum(
                    v for d1, d2, v in decayDat if isInvisible(d1, d2)
                )

        return {**widths, **decays}

    def parseTopDecays():

        """ 
        Get 2 body decays of the top, to catch cases where it decays to a Higgs 
        Note that in zslha the br are stored with integer keys     

        Here we want t -> H+ b. The original HiggsTools algorithm is a bit hokey, it could get confused if
        there was t -> H+ s at the same time that appears last. So I specifically require the b quark.    
        """
        topDecays = {
            #np.max(np.abs(x[-2:])): slha.br[6][x]
            np.max(np.abs(x)): slha.br[6][x]
            for x in slha.br[6]
            if len(x) == 2 and (5 in x or -5 in x) ## corrected: we want two-body decays to a Higgs here, so it assumes the 
        }        
        
        return {f"cxn_{p}_brtHpb_LHC13": topDecays.get(p, 0) for p in chargedPDGs}

    

    return {
        **parseMasses(),
        **parseMassUncs(),
        **parseDecays(),
        **parseCouplings(),
        **parseTopDecays(),
    }


"""

 Now the actual tool!

 
"""



class NewTool(HepTool):
    """ overload the init, name and settings are already given in HepTool """
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)
        #self.fake_uncertainties=True
        if 'Neutral Higgs' in self.settings and 'Charged Higgs' in self.settings:
            #### Override the options line ...
            
            #self.uncertainty_file='MHall_uncertainties.dat'
            self.neutralpdgids=list(self.settings['Neutral Higgs'])
            self.chargedpdgids=list(self.settings['Charged Higgs'])
            self.neutralhiggs=len(self.neutralpdgids)
            self.chargedhiggs=len(self.chargedpdgids)
            
            #if self.options == '':
            #    self.options='latestresults 2 effC %d %d' %(self.neutralhiggs,self.chargedhiggs)
        else:
            raise NameError('Must specify numbers of Charged and Neutral Higgs for HiggsTools')
        
            
        self.neutralIds = [str(i) for i in self.neutralpdgids]
        self.chargedIds = [str(i) for i in self.chargedpdgids]
        if 'Invisible Particles' in settings:
            self.Invisibles=settings['Invisible Particles']
        else:
            self.Invisibles=[]
        if 'HiggsBounds Dataset' not in settings:
            raise NameError("HiggsTools requires the location of the HiggsBounds Dataset")
        
        self.bounds = HB.Bounds(settings['HiggsBounds Dataset'])

        if 'HiggsSignals Dataset' not in settings:
            raise NameError("HiggsTools requires the location of the HiggsSignals Dataset")
        
        
        self.signals = HS.Signals(settings['HiggsSignals Dataset'])

        self.do_signals=False
        if 'Signal Strengths' in settings:
            self.do_signals=eval(settings['Signal Strengths'])


    def run(self, spc_file, temp_dir, log,data_point):
        
        if not os.path.exists(spc_file):
            log.debug('No spc; not running HiggsTools')
            raise
        
        log.debug('Reading SLHA file into HiggsTools')        

        ## This uses the built-in routine from HiggsTools which reads the file with pylha:        
        #all_inputs = readHB5SLHA(spc_file,self.neutralpdgids,self.chargedpdgids)

        if data_point.spc is None:
            data_point.spc= zslha.read(spc_file)

        all_inputs = readHB5ZSLHA(data_point.spc,self.neutralpdgids,self.chargedpdgids,invisiblePDGs=self.Invisibles)

        
        ## Fake uncertainty info if it is not there: do this in our homebrewed routine now
        """
        for neutH in self.neutralIds: ## in case the lightest one is lighter than the SM Higgs ...
            dmstr = 'dm_'+neutH
            if dmstr not in all_inputs:
                all_inputs[dmstr]=3.0
        """
        
        # Run HiggsTools
        log.debug('getting predictions from HiggsTools') 
        pred=predictionsFromDict(all_inputs, self.neutralIds, self.chargedIds, []) #, calcggH=True,calcHgamgam=False)

        if self.do_signals:
            predSM=HP.Predictions()
            hSM=predSM.addParticle(HP.NeutralScalar("hSM", "even"))
            #hSM.setMass(125.09)
            #HP.effectiveCouplingInput(hSM, HP.scaledSMlikeEffCouplings(1.0), reference="SMHiggsEW")

        hb,hs = self.bounds(pred), self.signals(pred)




        nobs=int(self.signals.observableCount())
        mypval=1 - chi2.cdf(hs,nobs)
        
        HBallowed=int(hb.allowed)



        maxobsratio=0
        maxdesc=""
        maxid=0
        maxnh=0

        for h in hb.selectedLimits:
            
            maxobsratio = max(maxobsratio,hb.selectedLimits[h].obsRatio())
            
        if data_point.spc is None:
            data_point.spc = zslha.SLHA()


        entries=OrderedDict()
        comments=OrderedDict()
        #data_point.spc.blocks['HIGGSTOOLS']=OrderedDict()
        #data_point.spc.blockcomments['HIGGSTOOLS']=OrderedDict()

        with open(spc_file,'a') as OF:
            OF.write('BLOCK HIGGSTOOLS \n')
            OF.write('1  %d  # HB result (0/1)\n' %HBallowed)
            OF.write('2  %10.8e  # HB maximum obs ratio\n' %maxobsratio)
            OF.write('11  %10.8e  # HS chi2\n' %hs)
            OF.write('12  %d  # HS number of observables\n' %nobs)
            OF.write('13  %10.8e  # BSMArt Inferred HS pval\n' %mypval)
            entries['1'] = HBallowed
            comments['1'] = 'HB result (0/1)'
            entries['2'] = maxobsratio
            comments['2'] = 'HB maximum obs ratio'
            entries['11'] = hs
            comments['11'] = 'HS chi2'
            entries['12'] = nobs
            comments['12'] = 'HS number of observables'
            entries['13'] = mypval
            comments['13'] = 'BSMArt Inferred HS pval'

            for h in hb.selectedLimits:
                tdesc="ID: %d, ref: %s: %s" %(hb.selectedLimits[h].limit().id(),hb.selectedLimits[h].limit().reference(),hb.selectedLimits[h].limit().processDesc())
                OF.write('%d  %10.8e  # %s\n' %(int(h),hb.selectedLimits[h].obsRatio(),tdesc))
                entries[h]=hb.selectedLimits[h].obsRatio()
                comments[h]=tdesc

            if self.do_signals:
                ## add channels

                ## see include/Higgs/predictions/Channels.hpp
                #prod_channels=['H','HZ','vbfH','HW','eeHZ']
                prod_channels=['H','HZ','vbfH','HW','eeHZ']
                colliders={'H':'LHC13','HZ':'LHC13','vbfH':'LHC13','HW':'LHC13','eeHZ':'LEP'}

                decay_channels=['bb','ZZ','WW','gamgam','gg','mumu','tautau']
                for nhiggs,nh in enumerate(self.neutralpdgids):
                    #chindex=(nhiggs+1)*1000
                    # need to reset mass to the hSM
                    try:
                        #temphmass=float(data_point.spc.Value("MASS",[nh]))
                        #hSM.setMass(temphmass)
                        hSM.setMass(float(data_point.spc.Value("MASS",[nh])))
                        HP.effectiveCouplingInput(hSM, HP.scaledSMlikeEffCouplings(1.0), reference="SMHiggsEW")
                    except Exception as e:
                        log.debug("Failed to reset neutral Higgs mass for "+str(nh)+": "+str(e))
                        continue
                    nhpred=pred.particle(str(nh))
                    for nchan,prodchan in enumerate(prod_channels):
                        #chindex=(nhiggs+1)*1000+nchan*100
                        for ndec,decchan in enumerate(decay_channels):
                            
                            try:
                                #mueff=nhpred.channelRate("LHC13",prodchan,decchan)/hSM.channelRate("LHC13",prodchan,decchan)
                                mueff=nhpred.channelRate(colliders[prodchan],prodchan,decchan)/hSM.channelRate(colliders[prodchan],prodchan,decchan)
                            except: ## presumably channel doesn't exist
                                mueff = 0.0
                                #continue
                            if np.isnan(mueff):
                                log.warning("NaN signal strength for %s -> %d -> %s" %(prodchan,int(nh),decchan))
                                mueff = 0.0
                            chindex=(nhiggs+1)*1000+nchan*100+ndec   # want to guarantee the numbers!!
                            entries[str(chindex)]=mueff
                            comments[str(chindex)]='mueff for  %s -> %d -> %s' %(prodchan,int(nh),decchan)
                            #print(' %d  %.8e  # %s\n' % (int(chindex),mueff,comments[str(chindex)]))  # for debugging ..
                            OF.write(' %d  %.8e  # %s\n' % (int(chindex),mueff,comments[str(chindex)]))
                            #chindex=chindex+1
                        
            data_point.spc.blocks['HIGGSTOOLS']=entries
            data_point.spc.blockcomments['HIGGSTOOLS']=comments
        
