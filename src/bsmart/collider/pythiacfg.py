"""
Script for editing a pythia config file based on a template

The idea is that each entry is of the form 

<string> = <value>

The code will attempt to evaluate the value based on the inputs given, which should be a function of the input variables of the code.


"""

import math
import os
import sys

class PYTHIACFG():
    def __init__(self):
        self.entries=[]

    def read(self,template_in):
        try:
            IF=open(template_in,'r')
        except Exception as e:
            raise NameError("Failed to open %s, %s" %(template_in,str(e)))


        for line in IF:
            linestrip=line.strip()
            if linestrip.startswith('!'):
                continue
            if '=' not in linestrip:
                continue
            if linestrip.startswith('SLHA:file'):  ## needs to be given by the code in LO mode
                continue
            # remove everything after a comment
            if '!' in linestrip:
                linestrip=linestrip.split('!')[0]

            ## need to be able to handle 23:onIfAny = 11 13 15 etc. 

            data=linestrip.split('=')

            
            #firstpart=data[0].replace(' ','')
            #secondpart=data[1].replace(' ','')
            firstpart=data[0].strip()
            secondpart=data[1].strip().replace('\t',' ')
            try:
                ff=float(secondpart)
                thirdpart=0
            except:
                thirdpart=1
            
            # if secondpart.isdigit():
            #     thirdpart='int'
            #     #secondpart=int(secondpart)
            # else:
            #     try:
            #         ff=float(st)
            #         thirdpart='float'
            #         #secondpart=ff
            #     except:
            #         thirdpart='string'
                    
            
            self.entries.append([firstpart,secondpart,thirdpart])

        IF.close()

    def write_dict(self,filename,var_dict={}):           
        with open(filename, 'w+') as f:
            for ll in self.entries:
                if ll[2] ==0:
                    f.write(ll[0]+' = '+ll[1]+'\n')
                else:
                    """
                        RHS doesn't evaluate to a float. So is either a string, or an expression, 
                        or a list of terms separated by whitespace.
                        new algorithm to see if any terms contain one of our variables

                    """
                    rhstr=''
                    terms=ll[1].split() 
                    for term in terms:
                        if any(x in term for x in var_dict):
                            try:
                                tval='%10.6e' % eval(term,var_dict)
                            except:
                                tval=term
                        else:
                            tval=term
                        rhstr+=' '+tval
                    f.write(ll[0]+' = '+rhstr+'\n')
                    """
                        try:
                            tval='%10.6e' % eval(ll[1],var_dict)
                        except:
                            tval=ll[1]
                        f.write(ll[0]+' = '+tval+'\n')
                    """
                        


    def write(self,filename,var_vals, var_names = None):
        if var_names is not None:
            var_dict = { var_name : var_val for var_name,var_val in zip(var_names,var_vals)}
            var_dict['VARIABLE'] = var_vals  ## Added for backward compatibility, maybe remove?
        else:
            var_dict={'VARIABLE':var_vals}   ## Added for backward compatibility, maybe remove?
        var_dict['math'] = math ## to allow maths functions such as sqrt or cos/sin/tan!
           
        with open(filename, 'w+') as f:
            for ll in self.entries:
                if ll[2] ==0:
                    f.write(ll[0]+' = '+ll[1]+'\n')
                else:
                    try:
                        tval='%10.6e' % eval(ll[1],var_dict)
                    except:
                        tval=ll[1]
                    f.write(ll[0]+' = '+tval+'\n')
