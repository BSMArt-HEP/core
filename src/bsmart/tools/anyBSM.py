"""
Tool to run anyBSM



https://anybsm.gitlab.io/

Provided by Martin Gabelmann -- thanks!!

"""

__meta__ = {
    "name": "anyBSM",
    "requires": ["anyBSM"],
    "settings": {
        "model": "Model name",
        "settings": "Dict of settings",
        "setparameters": "Dict of parameters"
    }
}



from anyBSM import anyBSM as anyBSM_lib
from bsmart.HEPRun import HepTool
from bsmart import zslha

class NewTool(HepTool):
    def __init__(self, name, settings, global_settings=None):
        
        HepTool.__init__(self, name, settings, global_settings)
        self.anybsm = anyBSM_lib(settings['model'],**settings['settings'])
        self.anybsm.warnSSSS = False
        self.anyinput = settings['setparameters']

        ## e.g. 'setparameters' : {} or 'setparameters': 'SPheno.spc.MSSM'

    def run(self, spc_file, temp_dir, log, data_point):
        #try:
        params = {k: eval(v,data_point.var_dict) for k,v in self.anyinput.items()}
        self.anybsm.setparameters(params)
        lam = self.anybsm.lambdahhh()
        lam_1 = lam['total'].real
        lam_0 = lam['treelevel'].real
        kaplam = lam_1/187.2818

        log.debug('AnyBSM parameters: '+str(params))
        if data_point.spc is None:
            data_point.spc = zslha.read(spc_file)
        with open(spc_file,'a') as OF:
            OF.write('BLOCK ANYBSM \n')
            OF.write('1  %10.8e  # lamhhh(0)\n' % lam_0)
            OF.write('2  %10.8e  # lamhhh(1)\n' % lam_1)
            OF.write('3  %10.8e  # kappa_lambda(1)\n' % kaplam)

        
        data_point.spc.blocks['ANYBSM'] = {'1': lam_0, '2': lam_1, '3': kaplam}
        data_point.spc.blockcomments['ANYBSM'] = {'1': 'lamhhh0', '2': 'lamhhh1', '3': 'kappalam'}
        
