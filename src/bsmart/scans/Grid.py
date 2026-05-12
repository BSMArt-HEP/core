"""
Simple grid scan based on 1D ranges of variables specified by the user.


Variables are specified in form of python functions, e.g.

.. code-block:: json

    "Variables": {
        "m0" : "np.geomspace(100,4000,10)",
        "tanb" : "np.linspace(2,50,10)"
    }


etc.

Will support in principle an unlimited number of variables constructing a hypercube, but of course this is not
very efficient!

"""

__meta__ = {
    "name": "Grid",
    "requires": ["numpy","itertools"],
    "settings": {
    }
}


from bsmart.core import Scan as Scan
import itertools
import numpy as np
import math
from bsmart import debug



class NewScan(Scan):
    """Scanner class for Grid Scans"""
 
    def __init__(self, inputs, log):
        Scan.__init__(self, inputs, log)
        self.MPI_mode= False
        if 'MPI' in inputs and inputs['MPI']['Size'] > 0: # MPI running!
            self.MPI_mode=True
        
            

    #def initialise(self):
        ## set the default to store tabbed output
        ## self.runsettings.tabbed_output=True
        
        
    def run(self):
        all_points=self.generate_parameter_points()

        """
        For MPI, I could be clever and have only the master node generating the points, and then
        distribute them among the nodes. Instead I will be stupid and have each node generate the
        whole set (which should be identical ...) but just work out which batch to run.
        This is easier to code but has disadvantages, e.g. at the moment if one node fails it just sits
        spinning on one core while waiting for the rest to finish. 
        """
        if self.MPI_mode:
            from mpi4py import MPI
            mpisize=int(MPI.COMM_WORLD.Get_size())
            mpirank=int(MPI.COMM_WORLD.Get_rank())
            if mpisize > 1: 
                self.log.info('MPI running! Core' + str(mpirank))
            n_input_data=len(all_points)
            chunksize=int(n_input_data/mpisize)
            chunkstart=mpirank*chunksize
            chunkend=min(chunkstart+chunksize,n_input_data)
            batch_points = all_points[chunkstart:chunkend]
            self.RunManager.run_batch(batch_points)
    
            MPI.COMM_WORLD.Barrier()
        else:
            self.RunManager.run_batch(all_points)
        
        # n_points=len(all_points)
        
        # if 'Batch_size' in self.inputs['Setup']:
        #     batch_size=eval(self.inputs['Setup']['Batch_size'])
        # else:
        #     batch_size=n_points

        # number_batches=math.ceil(n_points/batchsize)
        
        # for i in range(number_batches):
            
        
        
    def generate_parameter_points(self):
        vars=self.inputs['Variables']
        all_points = []
        for x in vars:
            all_points.append(eval(vars[x]))
        temp = list(itertools.product(*all_points))
        all_points = [list(xx) for xx in temp]
        return all_points       

    def postprocess(self,Point, observables, data_point,temp_dir,log, lock=None):
        """ No postprocessing result """
        return ''
