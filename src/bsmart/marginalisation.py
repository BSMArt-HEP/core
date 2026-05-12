from bsmart import debug

import math
from scipy import optimize
import numpy as np

"""
Scalers for marginal variables where we map from an interval onto (-infty, infty)
"""


"""
def downscaler_func_infty(vmin,vmax):
    diff = vmax-vmin
    a = (vmin+vmax)/2
    b = diff/2

    def tfunc(x):
        return a + b*math.tanh(x*1e6/diff)
        
    return tfunc


def upscaler_func_infty(vmin,vmax):
    diff = vmax-vmin
    a = (vmin+vmax)/2
    b = diff/2
    def tfunc(x):
        return diff*math.atanh((x-a)/b)
    return tfunc
"""     

def downscaler_func_infty(vmin,vmax):
    diff = vmax-vmin
    a = (vmin+vmax)/2
    b = diff/2

    def tfunc(x):
        return a + diff*math.atan(x)/math.pi
        
    return tfunc


def upscaler_func_infty(vmin,vmax):
    diff = vmax-vmin
    a = (vmin+vmax)/2
    b = diff/2
    def tfunc(x):
        return math.pi/2*math.tan((x-a)/b)
    return tfunc


#def chisq(x,mean,err):
#    return (x-mean)**2/(2.*err**2)



def chisq(mean,err):
    def tfunc(x):
        return (x-mean)**2/(err**2)

    return tfunc

def diff_for_brentq(mean,err):
    def tfunc(x):
        return (x-mean)/err
    return tfunc


def create_diff_marginal(inputs,log=None):
    tfs= [ diff_for_brentq(obs['TARGET'],obs['VARIANCE']) for obs in inputs.values()]
    def ttfunc(y):
        return tfs[0](y[0])
    return ttfunc



def create_chisquare_marginal(inputs,log=None):
    #print(inputs)
    
    nfuncs = len(inputs.values())

    #if nfuncs == 1:
    #    tfs= [ diff_for_brentq(obs['MEAN'],obs['VARIANCE']) for obs in inputs.values()]
    #    def ttfunc(y):
    #        return tfs[0](y[0])
    #else:
    myfuncs = [ chisq(obs['TARGET'],obs['VARIANCE']) for obs in inputs.values()]
    def ttfunc(y):
        return sum( [f(x) for f,x in zip(myfuncs,y)])/nfuncs

    return ttfunc





def create_scalers_marginal(inputs,log=None):
        downscalers=[]
        upscalers=[]
        ranges=[]
        for varname in inputs['Variables']:
           
            if 'RANGE' in inputs['Variables'][varname]:
                varmax=max(inputs['Variables'][varname]['RANGE'][0],inputs['Variables'][varname]['RANGE'][1])
                varmin=min(inputs['Variables'][varname]['RANGE'][0],inputs['Variables'][varname]['RANGE'][1])
                if varmin > varmax:
                    tt=varmax
                    varmax=varmin
                    varmin=tt
                ranges.append([varmin,varmax])
                diff=varmax-varmin
                if log is not None:
                    log.info('Creating scalers for %s between %.4e and %.4e' %(varname,varmin,varmax))
            
                downscalers.append(downscaler_func_infty(varmin,varmax))
                upscalers.append(upscaler_func_infty(varmin,varmax))

        return downscalers,upscalers,ranges




def find_a_sign_change(mindiff,max_safe_val,func,sols,start,end,ffa,ffb,minx=None,miny=None):
    """ Function is supposed to return a large value if unsafe, so we use max_safe_val """
    
    if ffa is None:
        fa=func(start)
    else:
        fa=ffa
    if ffb is None:
        fb=func(end)
    else:
        fb=ffb

    #print('start: %f, end: %f, fa: %f, fb: %f' %(start,end,fa,fb))
    if miny is None:
        
        if fa < max_safe_val:
            minx=start
            miny = ffa
        elif fb < max_safe_val:
            minx=end
            miny = ffb

    if fa < max_safe_val:
        if abs(fa) < abs(miny):
            minx=start
            miny = ffa

    if fb < max_safe_val:
        if abs(fb) < abs(miny):
            minx=end
            miny = ffb
        
    if fa < max_safe_val and fb < max_safe_val and np.sign(fa) != np.sign(fb):
        sols.append([start,end])
        return sols,True,minx,miny

    diff=end-start
    
    diff2=diff/2
    if diff2 < mindiff:
        return sols,False,minx,miny
    fc = func(start+diff2)
    
    sols,found1,minx,miny = find_a_sign_change(mindiff,max_safe_val,func,sols,start,start+diff2,fa,fc,minx,miny)
    found2 = False
    if not found1:
        sols,found2,minx,miny = find_a_sign_change(mindiff,max_safe_val,func,sols,start+diff2,end,fc,fb,minx,miny)
    if found1 or found2:
        return sols,True,minx,miny
    else:
        return sols,False,minx,miny


#def oned_usebrentq(func,ranges,log,tolerance=0.1,max_safe_val=1e30,intervals=10):
def oned_usebrentq(func,ranges,log,options={}):
    """ Takes a function with argument a 1D list for compatibility reasons """
    if 'ftol' in options:
        tolerance=options['ftol']
    else:
        tolerance=0.1

    if 'eps' in options:
        eps=options['eps']
    else:
        eps=0.05
    if 'max_safe' in options:
        max_safe_val=options['max_safe']
    else:
        max_safe_val=1e30
    
    if 'disp' in options:
        verbose=options['disp']
    else:
        verbose=False

    start=ranges[0][0]
    end=ranges[0][1]
    fstart=func([start])
    fend=func([end])
    
    def ffunc(x):
        return func([x])

    if fstart < max_safe_val and fend < max_safe_val and np.sign(fstart) != np.sign(fend):
        if verbose:
            print('Searching in range [%f,%f]' %(start,end))
        try:
            sol=optimize.brentq(ffunc,start,end,rtol=tolerance)
            if verbose:
                print('Found solution '+str(sol))
            return sol
        except Exception as e:
            raise NameError("Failed to find solution "+str(e))

    """ Search the range  """
    sols=[]
    #mindiff=abs(end-start)/intervals
    mindiff=eps*abs(end-start)
    sols,findarange,minx,miny=find_a_sign_change(mindiff,max_safe_val,ffunc,sols,start,end,fstart,fend)

    if not findarange:
        if verbose:
            print('Could not find a range! Using the minimum values %f, %f' %(minx,miny))
        log.debug('Could not find a range! Using the minimum values %f, %f' %(minx,miny))
        if minx is None:
            raise NameError("Could not find a single valid value in range")
        return minx
    else:
        log.debug('Found solutions: '+str(sols))
        if verbose:
            print('found solutions! '+str(sols))
    newstart=sols[0][0]
    newend=sols[0][1]
    if verbose:
        print('Newstart %f, fa: %f newend %f, fb: %f ' %(newstart,ffunc(newstart),newend,ffunc(newend)))
    #print(xopt)
    #print(newmin)
    
    
    try:
        sol=optimize.brentq(ffunc,newstart,newend,rtol=tolerance)
        log.debug('Found solution through brentq: '+str(sol))
        if verbose:
            print('Found solution! '+str(sol))
        return sol
    except Exception as e:
        raise NameError("Failed to find solution "+str(e))

    """

    if np.sign(fstart) != np.sign(fend) and fstart < 1e30 and fend < 1e30:
        try:
            sol=optimize.brentq(ffunc,start,end)
            return sol
        except Exception as e:
            raise NameError("Failed to find solution "+str(e))

    # evaluate on a series of points to find a good range

    

    try:
        xmin, fmin, ierr, numfunc=optimize.fminbound(ffunc,start,end,full_output=1)
        #xmin=optimize.fminbound(ffunc,start,end)
    except Exception as e:
        raise NameError("Failed to find minimum "+str(e))
   
    #print(xopt)
    #print(newmin)
    
    if np.sign(fstart) != np.sign(fmin):
        if choose_lower:
            newstart=start
            newend=xmin
        else:
            newstart=xmin
            newend=end
        try:
            sol=optimize.brentq(ffunc,newstart,newend)
            return sol
        except Exception as e:
            raise NameError("Failed to find solution "+str(e))
                
    else: ## can't find solution on this interval, at least not with this method
        ## maybe the problem is that we found a local minimum and not a true one
        ## either way it's time to give up
        raise NameError("Failed to find a solution")

    """
