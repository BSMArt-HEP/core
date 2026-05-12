"""

Basic functions related to likelihoods, given here for convenience for scans that want them

"""

import math


"""
def C_function(obs, obs_lower, obs_upper):
                        C_max = np.log(1 + np.finfo(np.float64).max) + 1
                        C_value = C_max
                        if isinstance(obs, str):
                                try:
                                        obs = float(obs)
                                except:
                                        raise Exception(f"Observable value {obs} not a string or float!")
                        if isinstance(obs, float):
                                C_value = np.log10(max(0.0, -obs + obs_lower, obs - obs_upper) + 1)
                        return C_value  
"""
# Helpers for interval-based CFUNCTION behaviour


def _normalize_bound(val):
    """
    Normalize user-specified bounds to floats or None.

    Accepts:
      - None                     -> None
      - float/int (incl. +/-inf) -> float(val)
      - strings like 'inf', '-inf', 'infty', 'infinity' -> +/-inf
    """
    if val is None:
        return None

    # Already numeric
    if isinstance(val, (int, float)):
        return float(val)

    # Try to interpret strings
    if isinstance(val, str):
        s = val.strip().lower()
        if s in ("inf", "+inf", "infty", "infinity", "+infty", "+infinity"):
            return float("inf")
        if s in ("-inf", "-infty", "-infinity"):
            return float("-inf")
        # fall back to numeric parse
        return float(s)

    # Last resort: try casting
    return float(val)


def LogC_interval(val, lower, upper):
    """
    Log C-function based on an interval [lower, upper] where:

      - If both bounds finite: closed interval [lower, upper]
      - If only lower is given: [lower, +inf[
      - If only upper is given: ]-inf, upper]
      - If neither is given: always inside (zero penalty)

    Returns 0.0 inside the allowed region and a negative log10(1 + distance)
    penalty outside.
    """
    try:
        d = 0.0
        if (lower is not None) and (upper is not None):
            # finite closed interval
            d = max(0.0, lower - val, val - upper)
        elif lower is not None:
            # [lower, +inf[
            d = max(0.0, lower - val)
        elif upper is not None:
            # ]-inf, upper]
            d = max(0.0, val - upper)
        else:
            # no bounds at all -> always inside
            d = 0.0

        return -math.log10(1.0 + d)
    except Exception:
        return -1e100

def LogClog_interval(val, lower, upper):
    """
    Like LogC_interval, but with log-scaled values.
    WARNING: val, lower, upper must be strictly positive numbers.
    """
    try:
        C_value = LogC_interval(math.log(val), math.log(lower), math.log(upper))
        return C_value
    except:
        return -1e100

def C_interval(val, lower, upper):
    """
    Non-log version corresponding to LogC_interval.
    """
    try:
        C_value = math.pow(10.0, LogC_interval(val, lower, upper))
        return C_value
    except:
        return -1e100

def Clog_interval(val, lower, upper):
    """
    Like C_interval, but with log-scaled values.
    WARNING: val, lower, upper must be strictly positive numbers.
    """
    try:
        C_value = C_interval(math.log(val), math.log(lower), math.log(upper))
        return C_value
    except:
        return -1e100

def LogC_function(val, m, s):
    """
    m is mean, s is variance, so
    lower = m-s, upper = m+s

    Gives a value between -1 * max log and 0;
    looks very much like a LL

    """
    lower = m - s
    upper = m + s
    try:
        C_value = LogC_interval(val, lower, upper)
        return C_value
    except:
        return -1e100

def LogClog_function(val, m, s):
    """
    Like LogC_function, but with log-scaled values.
    WARNING: val, m - s must be strictly positive numbers.
    """
    lower = m - s
    upper = m + s
    try:
        C_value = LogClog_interval(val, lower, upper)
        return C_value
    except:
        return -1e100


""" For MultiNest/Diver"""


# bias function to enhance likelihoods
def LogBias(val, m, s):
    return s * math.log(val / m)


# gauss likelihood
def LogGauss(val, m, s):
    # print('%f %f %f' %(val,m,s))
    return -((val - m) ** 2) / (2.0 * (s**2))


## Sigmoid for upper or lower bounds. The variance 's' controls how sharply it should turn off
def LogSigmoid(val, s):
    ## want - math.log(1.+math.exp(-val/s)) but safe for cases of large or small val/s
    rr = val / s
    if rr < -25:
        return rr + math.exp(rr)

    return -math.log(1.0 + math.exp(-val / s))


def SafeLog(val):
    try:
        return math.log(val)
    except:
        return -1e100  ## stupid I know


def SafeExp(val):
    try:
        return math.exp(val)
    except:
        return 1e100  ## could be modified




def MakeLikelihoodFuncLog(scaling, mean, var, lower_bound=None, upper_bound=None, log_scaling=False):
    tmean = mean
    tvar = var
    lb = lower_bound
    ub = upper_bound

    if scaling == "OFF":
        return lambda x: 0.0  # dummy function that returns 1
    # elif scaling == 'LOG':
    # 	return lambda x : logL(x, tmean, tvar)
    elif scaling == "UPPER":
        return lambda x: LogSigmoid(tmean - x, tvar)
    elif scaling == "LOWER":
        return lambda x: LogSigmoid(x - tmean, tvar)
    elif scaling == "BIAS":
        return lambda x: LogBias(x, tmean, tvar)
    elif (
        scaling == "USER"
    ):  #### User supplied likelihood/pval, e.g. from MicrOmegas, but we need to give a log
        return lambda x: SafeLog(x)
    elif (
        scaling == "SIGUSER"
    ):  ## Sigmoid of user-supplied function (gives value between 0 and 1, ideal for likelihood)
        return lambda x: LogSigmoid(x, 1.0)
    elif (
        scaling == "EXPUSER"
    ):  #### User supplied value needs exp to give L, but we need to give a log of that -> x
        return lambda x: x
    elif (
        scaling == "MINUSEXPUSER"
    ):  #### User supplied value needs exp to give L, but we need to give a log of that -> -x
        return lambda x: -x
    elif scaling == "LOGUSER":  ## This is confusing in this context
        return lambda x: -x
    elif scaling == "CFUNCTION":  # see above
        # For CFUNCTION, prefer explicit bounds if provided; otherwise fall back to m±s.
        if (lb is not None) or (ub is not None):
            # Normalize bounds once at creation time, not every call
            norm_lb = _normalize_bound(lb)
            norm_ub = _normalize_bound(ub)
            # log_scaling is a boolean flag that indicates whether the values should be log-scaled
            if log_scaling:
                return lambda x: LogClog_interval(x, norm_lb, norm_ub)
            else:
                return lambda x: LogC_interval(x, norm_lb, norm_ub)
        # log_scaling is a boolean flag that indicates whether the values should be log-scaled
        if log_scaling:
            return lambda x: LogClog_function(x, tmean, tvar)
        else:
            return lambda x: LogC_function(x, tmean, tvar)
    else:
        return lambda x: LogGauss(x, tmean, tvar)


# bias function to enhance likelihoods
def bias(val, m, s):
    return math.pow(val / m, s)


# gauss likelihood
def gauss(val, m, s):
    # print('%f %f %f' %(val,m,s))
    return math.exp(-((val - m) ** 2) / (2.0 * s**2))


# Log likelihood
def logL(val, m, s):
    return 1.0 / (1.0 + (val - m) ** 2 / (2.0 * s**2))


## Sigmoid for upper or lower bounds. The variance 's' controls how sharply it should turn off
def Sigmoid(val, s):
    return 1.0 / (1.0 + math.exp(-val / s))


def MakeLikelihoodFunc(scaling, mean, var, lower_bound=None, upper_bound=None, log_scaling=False):
    tmean = mean
    tvar = var
    lb = lower_bound
    ub = upper_bound

    if scaling == "OFF":
        return lambda x: 1.0  # dummy function that returns 1
    elif scaling == "LOG":
        return lambda x: logL(x, tmean, tvar)
    elif scaling == "UPPER":
        return lambda x: Sigmoid(tmean - x, tvar)
    elif scaling == "LOWER":
        return lambda x: Sigmoid(x - tmean, tvar)
    elif scaling == "BIAS":
        return lambda x: bias(x, tmean, tvar)
    elif scaling == "USER":  #### User supplied likelihood/pval, e.g. from MicrOmegas
        return lambda x: x
    elif (
        scaling == "SIGUSER"
    ):  ## Sigmoid of user-supplied function (gives value between 0 and 1, ideal for likelihood)
        return lambda x: Sigmoid(x, 1.0)
    elif scaling == "LOGUSER":  ## CONFUSING NOW!!
        return lambda x: SafeLog(x)
    elif (
        scaling == "MINUSEXPUSER"
    ):  # the provided function is an NLL, so to get a likelihood we need exp(-x)
        return lambda x: SafeExp(-x)
    elif (
        scaling == "EXPUSER"
    ):  # the provided function is an LL, so to get a likelihood we need exp(x)
        return lambda x: SafeExp(x)
    elif scaling == "CFUNCTION":  # see above
        if (lb is not None) or (ub is not None):
            # Normalize bounds once at creation time, not every call
            norm_lb = _normalize_bound(lb)
            norm_ub = _normalize_bound(ub)
            if log_scaling:
                return lambda x: Clog_interval(x, norm_lb, norm_ub)
            else:
                return lambda x: C_interval(x, norm_lb, norm_ub)
        if log_scaling:
            return lambda x: math.pow(10, LogClog_function(x, tmean, tvar))
        else:
            return lambda x: math.pow(10, LogC_function(x, tmean, tvar))
    else:
        return lambda x: gauss(x, tmean, tvar)


def MakeLikelihoods(observables_dict, loglike: bool = False):
    """
    returns a list of likelihood functions 	and masks (for variables that should be ignored by the likelihood computation)
    """
    out_likelihoods = []
    masks = []

    for obs in observables_dict.values():
        # print('adding observable ' +str(obs))
        if (
            "SCALING" not in obs
        ):  # will try to be Gaussian by default, but need some checks
            if "MEAN" not in obs or "VARIANCE" not in obs:
                scaling = (
                    "OFF"  ## just ignore the observable instead of raising an error
                )
                # raise NameError('MEAN or VARIANCE not declared for an observable')
            else:
                scaling = "GAUSS"
        else:
            scaling = obs["SCALING"]

        if scaling == "OFF":
            masks.append(False)
            continue
        else:
            masks.append(True)

        tmean=0.0
        tvar=0.0 # as a default

        lower_bound = None
        upper_bound = None
        log_scaling = str(obs.get("LOG_SCALING", False)).lower() == "true"

        if 'MEAN' in obs:
            tmean = obs['MEAN']
            if 'VARIANCE' in obs:
                tvar = obs['VARIANCE']
        else:
            tmean=0.0
            if 'RANGE' in obs:
                tmean=0.5*obs['RANGE'][0] + 0.5*obs['RANGE'][1]
                tvar=abs(0.5*obs['RANGE'][1] - 0.5*obs['RANGE'][0])


        # CFUNCTION-specific bounds logic layered on top
        if scaling == "CFUNCTION":
            lower_bound = obs.get("LOWER_BOUND", None)
            upper_bound = obs.get("UPPER_BOUND", None)
            if (lower_bound is not None) or (upper_bound is not None):
                # Explicit bounds take precedence; mean/var ignored by factories
                tmean = None
                tvar = None
        # if 'VARIANCE' in obs:
        # 	tvar = obs['VARIANCE']
        # else:
        # 	tvar = 0.0
        # print(f"Creating likelihood function with {scaling},{tmean},{tvar},{lower_bound},{upper_bound}")
        if loglike:
            out_likelihoods.append(
                MakeLikelihoodFuncLog(scaling, tmean, tvar, lower_bound, upper_bound, log_scaling)
            )
        else:
            out_likelihoods.append(
                MakeLikelihoodFunc(scaling, tmean, tvar, lower_bound, upper_bound, log_scaling)
            )

    return out_likelihoods, masks


def safe_float(x):
    """Convert x to float if possible, else return NaN."""
    try:
        val = float(x)
        # normalize strings like 'NaN' -> NaN
        return val if not math.isnan(val) and not math.isinf(val) else math.nan
    except (ValueError, TypeError):
        return math.nan


def smooth_cap_loss(x,maxloss):
    """
    Caps the loss by applying a sigmoid.
    This is useful for losses that are unbounded.                   

    """
    return -math.copysign(maxloss,x)*math.expm1(-abs(x)/maxloss) 



def MakeGlobalLikelihood(
    observables_dict: dict, 
    user_loglike: bool = False, 
    return_type: str = None,
    max_loss: float = None,
    smooth_losses: bool = True
):
    """
    Makes a function that takes a list of observables and outputs a global likelihood
    Also returns a 'mask' list corresponding to whether a given observable is ignored or not;
    useful for machine-learning tasks.

    This seems to work ok. I wonder whether it would be safer to wrap it into a class that stores the functions?

    user_loglike is a legacy feature. If True it will use log likelihoods; if False it will use regular likelihoods.

    The idea will be that the user tells the code what to provide to give a *likelihood*.
    So if what is provided is a log likelihood, they should provide EXPUSER; if it is an NLL it should be MINUSEXPUSER.
    This may then be undone by the requirements of the scan.


    return_type should be:
            "Likelihood" (A likelihood), by default
            "LL" (Log Likelihood)
            "NLL" (Negative Log Likelihood)

  

    smooth_losses is True by default. If true, then the function will *always* return a number (it should not raise an error if some observables are infinite or NaN),
    and in the case of log likelihoods the value is capped at max_loss.

    max_loss is optional and only relevant for LL and NLL methods and when smooth_losses is True. It will cap the loss at this value (if provided), smoothing using a sigmoid. 
    Otherwise a maximum value corresponding to  710.782712893384  ( = math.log(1+sys.float_info.max) + 1 on a standard system ) is used. 

    """
    out_likelihoods = []
    # masks=[]
    loglike = user_loglike
    multiplier = 1.0
    if return_type is not None:
        if return_type.upper() == "LIKELIHOOD":
            loglike = False
        elif return_type.upper() == "LL":
            loglike = True
        elif return_type.upper() == "NLL":
            loglike = True
            multiplier = -1.0

    out_likelihoods, masks = MakeLikelihoods(observables_dict, loglike)

	# for obs in observables_dict.values():
	# 	#print('adding observable ' +str(obs))
	# 	if 'SCALING' not in obs: # will try to be Gaussian by default, but need some checks
	# 		if 'MEAN' not in obs or 'VARIANCE' not in obs:
	# 			scaling = 'OFF' ## just ignore the observable instead of raising an error
	# 			#raise NameError('MEAN or VARIANCE not declared for an observable')
	# 		else:
	# 			scaling = 'GAUSS'
	# 	else:
	# 		scaling=obs['SCALING']

	# 	if scaling == 'OFF':
	# 		masks.append(False)
	# 		continue
	# 	else:
	# 		masks.append(True)


	# 	tmean=0.0
	# 	tvar=0.0
	# 	if 'MEAN' in obs:
	# 		tmean = obs['MEAN']
	# 		if 'VARIANCE' in obs:
	# 			tvar = obs['VARIANCE']
	# 	else:
	# 		#tmean=0.0
	# 		if 'RANGE' in obs:
	# 			tmean=0.5*obs['RANGE'][0] + 0.5*obs['RANGE'][1]
	# 			tvar=abs(0.5*obs['RANGE'][1] - 0.5*obs['RANGE'][0])
		
		
	# 	if loglike:
	# 		out_likelihoods.append(MakeLikelihoodFuncLog(scaling,tmean,tvar))
	# 	else:
	# 		out_likelihoods.append(MakeLikelihoodFunc(scaling,tmean,tvar)) 

    if loglike:
        if smooth_losses:
            if max_loss is None:
                maxloss = 710.782712893384 # math.log(1+sys.float_info.max) + 1 or np.log(1 + np.finfo(np.float64).max) + 1
            else:
                maxloss = max_loss

            def globallh(observables):
                if len(observables) < len(masks):  # case where we only have provided 'live' observables
                    return multiplier * sum(
                        f(float(v)) for v, f in zip(observables, out_likelihoods)
                    )
                else:
                    likeit = iter(out_likelihoods)
                    return multiplier*sum([smooth_cap_loss((next(likeit))(val),maxloss) if mask and not math.isnan(val := safe_float(v)) else float((next(likeit) and False) or -maxloss) 
                                            for v, mask in zip(observables, masks) if mask])
        else: # no smoothing, no safety checks
            def globallh(observables):
                lh = 0.0
                if len(observables) < len(masks):  # case where we only have provided 'live' observables
                    return multiplier * sum(
                        f(float(v)) for v, f in zip(observables, out_likelihoods)
                    )
                    # for  v, f in zip(observables,out_likelihoods):
                    # 	lh = lh + f(float(v))
                else:  # case where we have all observables
                    likeit = iter(out_likelihoods)
                    return multiplier * sum(
                        next(likeit)(float(v))
                        for v, mask in zip(observables, masks)
                        if mask
                    )
        # def globallh(observables):
        #     lh = 0.0
        #     if len(observables) < len(
        #         masks
        #     ):  # case where we only have provided 'live' observables
        #         return multiplier * sum(
        #             f(float(v)) for v, f in zip(observables, out_likelihoods)
        #         )
        #         # for  v, f in zip(observables,out_likelihoods):
        #         # 	lh = lh + f(float(v))
        #     else:  # case where we have all observables
        #         likeit = iter(out_likelihoods)
        #         return multiplier * sum(
        #             next(likeit)(float(v))
        #             for v, mask in zip(observables, masks)
        #             if mask
        #         )


				# ##return sum((f := next(likeit))(val) for v, mask in zip(observables, masks) if mask and not math.isnan(val := safe_float(v)))
				# return sum([
    			# 	(f := next(likeit))(val) if mask and not math.isnan(val := safe_float(v)) else (next(likeit) and 0 if mask else 0)
    			# 	for v, mask in zip(observables, masks)
				# 	])


				# # equivalent to:

				# total = 0
				# for v, mask in zip(observables, masks):
    			# 	if mask:
				# 		f=next(likeit)
				# 		val = safe_float(v)
				# 		if not math.isnan(val):
				# 			total += f(val)
						
				# return total

				# return sum(
    			# 	next(likeit)(val) if mask and not math.isnan(val) else 0
   				# 	for v, mask in zip(observables, masks)
    			# 	for val in [safe_float(v)]
				# 	)


    else:  # no loglike, no multiplier
        # smoothing is irrelevant here because we just return zero if we have an inf. 
        # BUT we will have two versions: one which will raise an error in case of inf (no smoothing) and one that returns zero likelihood

        
        def globallh(observables):
                lh = 1.0
                if len(observables) < len(masks):  # case where we only have provided 'live' observables
                    for v, f in zip(observables, out_likelihoods):
                        lh *= f(float(v))
                else:  # case where we have all observables
                    likeit = iter(out_likelihoods)
                    for v, mask in zip(observables, masks):
                        if mask:
                            #lh *= next(likeit)(float(v))
                            #lh *= next(likeit)(val) if not math.isnan(val := safe_float(v)) else 0
                            val = safe_float(v)
                            if math.isnan(val):
                                if smooth_losses:
                                    return 0.0
                                else:
                                    raise
                            lh *= next(likeit)(val)
                        
                    return lh

    # return lambda x : globallh(x), masks
    return globallh, masks
