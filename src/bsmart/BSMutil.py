"""
Some simple utilities that can be used either in the scans or in external programs

So far this includes routines for scaling up and down variables. 
To invoke them, use

Import BSMutil

in your scan. 

Then 


"""


from bsmart import debug
import math
import scipy

def downscaler_func(vmin,vmax):
    """
        Inverse of uniform prior, maps [vmin,vmax] -> [0,1]
    """
    diff = vmax-vmin
    
    def tfunc(x):
        res = (x-vmin)/diff
        if res < 0:
            return 0.0
        elif res > 1.0:
            return 1.0
        else:
            return res
        
    return tfunc


def upscaler_func(vmin,vmax):
    """
        Uniform prior on [0,1] -> [vmin,vmax]
    """
    diff = vmax-vmin
    def tfunc(x):
        if x < 0:
            return vmin
        elif x > 1:
            return vmax
        else:
            return vmin+x*diff
    return tfunc


def upscaler_capped(vmin,vmax):
    """
        Uniform prior on the interval [0,1] but raises an error if outside
    """
    diff = vmax-vmin
    def tfunc(x):
        if x < 0:
            raise
        elif x > 1:
            raise
        else:
            return vmin+x*diff
    return tfunc

def downscaler_capped(vmin,vmax):
    """
        Uniform prior on the interval [0,1] but raises an error if outside
    """
    diff = vmax-vmin
    
    def tfunc(x):
        res = (x-vmin)/diff
        if res < 0:
            raise
        elif res > 1.0:
            raise
        else:
            return res
        
    return tfunc

def downscaler_open(vmin,vmax):
    """
        Rescale a range onto [0,1] but allow values outside
    """
    diff = vmax-vmin
    
    def tfunc(x):
        res = (x-vmin)/diff
        return res
        
    return tfunc


def upscaler_open(vmin,vmax):
    """
        Rescale a range onto [0,1] but allow values outside
    """
    diff = vmax-vmin
    def tfunc(x):
        return vmin+x*diff
    return tfunc

def downscaler_symlinear(vmin,vmax):
    """
        Inverse of Map [0.5,1] -> [|vmin|,|vmax|] and [0,0.5] -> [-|vmax|,-|vmin|]
        This should be a disjoint interval. 
    """
    avmin = min(abs(vmin),abs(vmax))
    avmax = max(abs(vmin),abs(vmax))

    diff = (avmax-avmin)
    
    def tfunc(x):
        ax = abs(x)
        #_sign = math.copysign(1.0, _box)  # Set _sign equal to the sign of _box
        #res = (abs(x)-avmin)/diff


        if ax < avmin:  # out of lower bound -> map to 0.5
            return 0.5
        elif ax > avmax: # out of upper bound
            return math.copysign(1.0,x)
        else:
            # this gives _box = 2 * r -1
            _z = math.copysign((ax-avmin)/diff,x)
            return 0.5*(_z + 1.0)
        
    return tfunc


def upscaler_symlinear(vmin,vmax):
    """
        Map [0.5,1] -> [|vmin|,|vmax|] and [0,0.5] -> [-|vmax|,-|vmin|]
        This should be a disjoint interval. 
    """
    avmin = min(abs(vmin),abs(vmax))
    avmax = max(abs(vmin),abs(vmax))

    diff = (avmax-avmin)
    def tfunc(x):
        _box = 2*x -1   # in [-1,1]. Take the [0,1] part and do a normal map; the [-1,0] maps to opposite sign
        #_sign = math.copysign(1.0, _box)  # Set _sign equal to the sign of _box
        if x < 0:
            return -avmax
        elif x > 1:
            return avmax
        else:
            return math.copysign(avmin+abs(_box)*diff,_box)
    return tfunc



def upscaler_logprior(vmin,vmax):
    """ Thanks to Luc Darme"""
    if vmin == 0:
        raise 
    if vmin < 0 and vmax > 0:
        raise
    logratio=math.log(vmax/vmin)

    def tfunc(x):
        if x<= 0:
            return vmin
        if x >= 1:
            return vmax
        #return vmin*math.exp(x*math.log(vmax/vmin))
        return vmin*math.exp(x*logratio)
    return tfunc


def downscaler_logprior(vmin,vmax):
    """
        x = log(y/vmin)/log(vmax/vmin) 
    
    """
    if vmin == 0:
        raise 
    if vmin < 0 and vmax > 0:
        raise

    
    logratio=math.log(vmax/vmin)

    def tfunc(x):
        y=x/vmin

        if y<= 0:  # need log(x/vmin) so no good
            return 0.0
        y = math.log(y)/logratio
        if y >= 1:
            return 1.0
        return y
    return tfunc

def upscaler_symlog(vmin,vmax):
    """  

        _box = 2 * box - 1
        _sign = np.sign(_box)
        _exp = low + np.abs(_box) * (high - low)
        return float(_sign * 10**_exp)
  
        Map [0.5,1] -> vmin *exp (x*log(vmax/vmin)) and [0,0.5] -> - vmin *exp (x*log(vmax/vmin))
        This should be a disjoint interval. 
   
    avmin = min(abs(vmin),abs(vmax))
    avmax = max(abs(vmin),abs(vmax))

    diff = avmax-avmin
    """
    #if vmin == 0:
    #    raise 
    if vmin * vmax <= 0 or vmin < 0: # require both positive and non-zero
        raise NameError("Incorrect bounds for symlog: [%f,%f]" %(vmin,vmax))

    # make sure we have correct ordering ...
    avmin = min(vmin,vmax)
    avmax = max(vmin,vmax)
    logratio=math.log(avmax/avmin)

    def tfunc(x):
        
        if x < 0:
            return -avmax
        elif x > 1:
            return avmax
        else:
            _box = 2*x - 1
            return math.copysign(avmin*math.exp(abs(_box)*logratio),_box)
        
    return tfunc


def downscaler_symlog(vmin,vmax):
    """
        Inverse of Map [0.5,1] -> vmin *exp (x*log(vmax/vmin)) and [0,0.5] -> - vmin *exp (x*log(vmax/vmin))
        This should be a disjoint interval. 
        _box = log(y/vmin)/log(vmax/vmin) 
        _box = 2x -1
        x = (_box + 1)/2
    
    """
    if vmin * vmax <= 0 or vmin < 0: # require both positive and non-zero
        raise NameError("Incorrect bounds for symlog: [%f,%f]" %(vmin,vmax))
    
    # make sure we have correct ordering ...
    avmin = min(vmin,vmax)
    avmax = max(vmin,vmax)
    logratio=math.log(avmax/avmin)

    def tfunc(x):
        ay=abs(x/avmin)

        if ay < 1.0:  # below the minimum value, so return 
            return 0.5

        #if y<= 0:  # need log(x/vmin) so no good
        #    return 0.0
        # ay should be in [0,1]
        ay = math.log(ay)/logratio
        if ay >= 1:
            return math.copysign(1.0,x)
        # map to [-1,1] and then to [0,1]
        return 0.5*(1.0+math.copysign(ay,x))
    return tfunc


def upscaler_gaussprior(vmin,vmax,pmean,psigma):
    sqrt2=math.sqrt(2.0)

    def tfunc(x):
        #if x <= 0 or x >= 1:
        #    return 0.0
        
        res=pmean + psigma*sqrt2*scipy.special.erfcinv(2*(1.0-x))
        if res < vmin:
            return vmin
        if res > vmax:
            return vmax
        return res
    return tfunc

def upscaler_lognormalprior(vmin,vmax,pmean,psigma):
    sqrt2=math.sqrt(2.0)
    psigma2=psigma**2

    def tfunc(x):
        
        #if x <= 0 or x >= 1:
        #    return 0.0
        #bracket = psigma2 + psigma*sqrt2*math.erf(2*x)
        bracket = psigma2 + psigma*sqrt2*scipy.special.erfcinv(2*x)
        res=pmean*math.exp(bracket)
        if res < vmin:
            return vmin
        if res > vmax:
            return vmax
        return res 
    return tfunc

def upscaler_expprior(vmin,vmax,psigma):
    """ P(y) ~ e^y. Coincides with MultiNest version when vmin = -infty, vmax=0, psigma > 0 or vmin=infty, vmax=0, psigma < 0 """
    
    ##x = \sigma^{-1}\log[e^{x_1 \sigma} + y(e^{\sigma x_2} - e^{\sigma x_1}) ]
    # sigma = -lambda
    exp1=math.exp(vmin*psigma)
    exp2=math.exp(vmax*psigma)


    def tfunc(x):
        #if x <= 0 or x >= 1:
        #    return 0.0
        
        #res = - math.log(x)/pmean
        res= math.log(exp1 + x*(exp2-exp1))/psigma
        if res < vmin:
            return vmin
        if res > vmax:
            return vmax
        return res 
    return tfunc


def upscaler_sinprior(vmin,vmax):

    #deg2rad=0.017453292
    deg2rad=math.pi/180.0
    cx1=math.cos(vmin*deg2rad)
    cx2=math.cos(vmax*deg2rad)
    def tfunc(x):
        #if x <= 0 or x >= 1:
        #    return 0.0
        res=math.acos(cx1+x*(cx2-cx1))
        return res 
    return tfunc

def upscaler_cauchyprior(vmin,vmax,pmean,psigma):

    
    def tfunc(x):
        #if x <= 0 or x >= 1:
        #    return 0.0
        res=pmean+psigma*math.tan(math.pi*(x-0.5))
        if res < vmin:
            return vmin
        if res > vmax:
            return vmax
        return res 
    return tfunc


def create_scalers(inputs,log=None):
        downscalers=[]
        upscalers=[]
        for varname in inputs['Variables']:
            varmax=math.inf # allow infinite range by default
            varmin=-math.inf
            prior = 'UNIFORM'
            prior_dict={}
            if 'RANGE' in inputs['Variables'][varname]:
                varmax=max(inputs['Variables'][varname]['RANGE'][0],inputs['Variables'][varname]['RANGE'][1])
                varmin=min(inputs['Variables'][varname]['RANGE'][0],inputs['Variables'][varname]['RANGE'][1])
                if varmin > varmax:
                    tt=varmax
                    varmax=varmin
                    varmin=tt

                diff=varmax-varmin
            else:
                prior = 'NONE' # no range -> not compatible with a UNIFORM prior, but it might be ok for a GAUSSIAN one etc (which can have infinite range)
                if log is not None:
                    log.debug('Variable %s has no range: using unscaled inputs!'%(varname))
            
                #downscalers.append(downscaler_func(varmin,varmax))

            if 'PRIOR' in inputs['Variables'][varname]:
                if hasattr(inputs['Variables'][varname]['PRIOR'],"len"):
                    prior_dict=inputs['Variables'][varname]['PRIOR']
                    if 'TYPE' in prior_dict:
                        prior= prior_dict['TYPE'].upper()
                else:
                    prior = inputs['Variables'][varname]['PRIOR'].upper()

                if log is not None:
                    log.info('Creating scalers for %s between %.4e and %.4e with prior %s'%(varname,varmin,varmax,prior))

            try:
                    if prior == 'UNIFORM':
                        upscalers.append(upscaler_func(varmin,varmax))
                        downscalers.append(downscaler_func(varmin,varmax))
                    elif prior == 'SYMLINEAR':
                        upscalers.append(upscaler_symlinear(varmin,varmax))
                        downscalers.append(downscaler_symlinear(varmin,varmax))
                    elif prior == 'INTERVAL':
                        upscalers.append(upscaler_capped(varmin,varmax))
                        downscalers.append(downscaler_capped(varmin,varmax))
                    elif prior == 'OPEN':
                        upscalers.append(upscaler_open(varmin,varmax))
                        downscalers.append(downscaler_open(varmin,varmax))
                    elif prior =='LOG':
                        upscalers.append(upscaler_logprior(varmin,varmax))
                        downscalers.append(downscaler_logprior(varmin,varmax))
                    elif prior == 'SYMLOG':
                        upscalers.append(upscaler_symlog(varmin,varmax))
                        downscalers.append(downscaler_symlog(varmin,varmax))
                    elif prior == 'GAUSS':
                        prior_mean=float(prior_dict['MEAN'])
                        prior_sigma=float(prior_dict['SIGMA'])
                        upscalers.append(upscaler_gaussprior(varmin,varmax,prior_mean,prior_sigma))
                        downscalers.append(downscaler_func(varmin,varmax))  ### TO UPDATE!!
                    elif prior == 'LOGNORMAL':
                        prior_mean=float(prior_dict['MEAN'])
                        prior_sigma=float(prior_dict['SIGMA'])
                        upscalers.append(upscaler_lognormalprior(varmin,varmax,prior_mean,prior_sigma))
                        downscalers.append(downscaler_func(varmin,varmax))  ### TO UPDATE!!
                    elif prior == 'EXP':
                        prior_sigma=float(prior_dict['SIGMA'])
                        upscalers.append(upscaler_expprior(varmin,varmax,prior_sigma))
                        downscalers.append(downscaler_func(varmin,varmax))  ### TO UPDATE!!
                    elif prior == 'SIN':
                        upscalers.append(upscaler_sinprior(varmin,varmax))
                        downscalers.append(downscaler_func(varmin,varmax))  ### TO UPDATE!!
                    elif prior == 'CAUCHY':
                        prior_mean=float(prior_dict['MEAN'])
                        prior_sigma=float(prior_dict['SIGMA'])
                        upscalers.append(upscaler_cauchyprior(varmin,varmax,prior_mean,prior_sigma))
                        downscalers.append(downscaler_func(varmin,varmax))  ### TO UPDATE!!
                    elif prior == 'NONE':  # don't transform!! Might be dangerous?
                        upscalers.append(lambda x: x)
                        downscalers.append(lambda x: x)
                    elif prior == 'FIXED': # fix at a constant value ... merely for debugging
                        
                        mean = float(prior_dict.get('MEAN',0.5*(varmin+varmax)))
                        if math.isnan(mean): # if we haven't set a range, just set the 
                            raise NameError('FIXED prior used for variable %s without a valid RANGE or MEAN' % varname)

                        upscalers.append(lambda x: mean)
                        downscalers.append(lambda x: 0.5)
                    else:
                        upscalers.append(upscaler_func(varmin,varmax))
                        downscalers.append(downscaler_func(varmin,varmax))
            except:
                    if log is not None:
                        log.critical('Problem with definition of priors')
                    raise NameError('Failed to create upscalers')
            
            """
            else: # no range: still need a scaler! So don't!
                if log is not None:
                    log.error('Variable %s has no range: using unscaled inputs!'%(varname))
                upscalers.append(lambda x: x)
                downscalers.append(lambda x: x)
            """
            """
                else:
                    upscalers.append(upscaler_func(varmin,varmax))
                    downscalers.append(downscaler_func(varmin,varmax))
            """
        return downscalers,upscalers




