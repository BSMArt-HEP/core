"""
Routines to automatically generate simple plots using matplotlib

Function make_plots takes the scan as argument

In the settings the plots to be produced should just be pairs of variables and observables, or the 'result'. But the result is not actually stored in the point in the RunManager ... only in the csv/spectrum file etc.

We could have three strategies:
1) Take the points directly from the store in the runner in the_scan.RunManager.all_points or the_scan.RunManager.valid_points (not necessarily always there for every scan)
2) Read from csv file using pandas (easiest)
3) Read from the output_file which stores all spectrum files if the appropriate option is set

For other things (e.g. using only the selected points in an MCMC scan) this is scan dependent and we can't be sure how or where the results are stored.

It would be useful to write a script that generates the 


"""



from bsmart import core
from bsmart import HEPRun

import os
from bsmart import debug

import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

def say_hello(the_scan):
    print('Hello from BSMplots!')
    for var in the_scan.inputs['Variables']:
        print(var)

def make_plot(xlabel,ylabel,x,y,name,xlim,ylim,plotdict={}):
    fig=plt.figure(figsize=(7,6),constrained_layout='True')
    plt.grid('on')
    #print(x)
    
    plt.scatter(x,y)
    plt.xlabel(xlabel,fontsize=13)
    plt.ylabel(ylabel,fontsize=13)

    if 'Scaling' in plotdict:
        xscale=plotdict['Scaling'][0]
        yscale=plotdict['Scaling'][1]
        plt.xscale(xscale)
        plt.yscale(yscale)

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.savefig(name)
    plt.close()

     
def make_auto_corner_csv(namestub,csvfile,shorts,plotdict):
    #print(','.join(list(map(str,shorts))))
    with open(namestub+'.py','w') as OF:
        OF.write('import matplotlib.pyplot as plt\n')
        OF.write('import matplotlib as mpl\n')
        OF.write('import seaborn as sns\n')
        OF.write('import numpy as np\n')
        OF.write('import pandas as pd\n')
        OF.write('import os\n\n')
        nvars=len(shorts)
        vars_string="\""+"\",\"".join(list(map(str,shorts)))+"\""
        if 'Labels' in plotdict:
            labels_string="\""+"\",\"".join(list(map(str,plotdict["Labels"])))+"\""
        else:
            labels_string=vars_string
        
        #OF.write("fig=plt.figure(figsize=(7,6),constrained_layout='True')\n")
        OF.write("fig, axes = plt.subplots(nrows=%d, ncols=%d, figsize=(7,6))\n"%(nvars,nvars))
        OF.write("vars=[%s]\n" % labels_string)
        OF.write("input_data=pd.read_csv('"+csvfile+"')\n")
        OF.write("data=pd.DataFrame(input_data,columns=["+vars_string+"]).to_numpy()\n")
        
        OF.write("for i in range(%d):\n" %nvars)
        #OF.write("    axes[0,i].set_xlabel(\"%s\" % vars[i])\n") 
        OF.write("    for j in range(%d):\n"%nvars)
        OF.write("        axes[j,i].set_xlabel(\"\")\n") 
        OF.write("        if i > 0:\n")
        OF.write("            axes[j,i].set_yticklabels([])\n") 
        OF.write("        if j < %d:\n" %(nvars-1))
        OF.write("            axes[j,i].set_xticklabels([])\n")
        #OF.write("        axes[j,0].set_ylabel(\"%s\" % vars[j])\n") 
        OF.write("        axes[j,i].grid()\n")
        OF.write("        if i == j:\n")
        OF.write("            sns.kdeplot(data[:, i], ax=axes[j,i])\n")
        OF.write("            if i > 0:\n")
        OF.write("                axes[j,i].set_ylabel(\"\")\n")
        #OF.write("            axes[j,i].set_xlabel(\"%s\" % vars[i])\n")                       
        OF.write("        elif i < j:\n")
        OF.write("            sns.scatterplot(x=data[:, i], y=data[:, j], ax=axes[j,i], alpha=0.4)\n")
        OF.write("            if i == 0:\n")
        OF.write("                axes[j,i].set_ylabel(\"%s\" % vars[j])\n")
        OF.write("            else:\n")
        OF.write("                axes[j,i].set_ylabel(\"\")\n")
        #OF.write("            axes[j,i].set_xlabel(\"%s\" % vars[i])\n")
        #OF.write("            axes[j,i].set_ylabel(\"%s\" % vars[j])\n")                            
        OF.write("        else:\n")
        OF.write("            axes[j,i].set_visible(False)\n")
        OF.write("    axes[j,i].set_xlabel(\"%s\" % vars[i])\n")
        OF.write("plt.tight_layout()\n" )
        OF.write("plt.savefig(\"%s.png\")\n" %namestub)

def make_fancy_auto_corner_csv(namestub,csvfile,shorts,plotdict):
    
    with open(namestub+'.py','w') as OF:
        OF.write('try:\n')
        OF.write('    import corner\n')
        OF.write('except ImportError:\n')
        OF.write('    exit("Please install corner")\n')
        OF.write('import numpy as np\n')
        OF.write('import pandas as pd\n')
        nvars=len(shorts)
        vars_string="\""+"\",\"".join(list(map(str,shorts)))+"\""
        if 'Labels' in plotdict:
            labels_string="\""+"\",\"".join(list(map(str,plotdict["Labels"])))+"\""
        else:
            labels_string=vars_string
        
        #OF.write("fig=plt.figure(figsize=(7,6),constrained_layout='True')\n")
        
        OF.write("vars=[%s]\n" % labels_string)
        OF.write("input_data=pd.read_csv('"+csvfile+"')\n")
        OF.write("data=pd.DataFrame(input_data,columns=["+vars_string+"]).to_numpy()\n")
        OF.write("fig= corner.corner(data,labels=vars,\nplot_datapoints=False,\nfill_contours=True,\nplot_density=False,\ncolor='black',\nhist_kwargs={\"lw\": 1, \"fill\": False},\ncontourf_kwargs={\"cmap\": \"Purples\", \"colors\": None},\ncontour_kwargs={\"colors\": \"#4a007a\", \"linewidths\": 1.5})\n")
        OF.write("fig.savefig('"+namestub+".pdf')\n")



def make_script_csv(xshort,yshorts,xlabel,ylabel,csvfile,namestub,xlim,ylim,plotdict={}):
    #print("Making script "+namestub+'.py')
    with open(namestub+'.py','w') as OF:
        OF.write('import matplotlib.pyplot as plt\n')
        OF.write('import matplotlib as mpl\n')
        OF.write('import numpy as np\n')
        OF.write('import pandas as pd\n')
        OF.write('import os\n\n')
        OF.write("fig=plt.figure(figsize=(7,6),constrained_layout='True')\n")
        OF.write("plt.grid('on')\n")

        if 'Scaling' in plotdict:
            xscale=plotdict['Scaling'][0]
            yscale=plotdict['Scaling'][1]
            OF.write('plt.xscale(\"'+xscale+'\")\n')
            OF.write('plt.yscale(\"'+yscale+'\")\n')
        OF.write("input_data=pd.read_csv('"+csvfile+"')\n")
        OF.write("data=pd.DataFrame(input_data,columns=['"+xshort+"'")
        for yshort in yshorts:
            OF.write(",'"+yshort+"'")
        OF.write("]).to_numpy()\n")


        OF.write("plt.xlabel(r'"+xlabel+"',fontsize=13)\n")
        OF.write("plt.ylabel(r'"+ylabel+"',fontsize=13)\n")
        if xlim is not None:
            OF.write("plt.xlim("+str(xlim)+")\n")  
            OF.write('xmin=%s\n' %str(xlim[0]))
            OF.write('xmax=%s\n' %str(xlim[1]))
            
        
            
        if ylim is not None:
            OF.write("plt.ylim("+str(ylim)+")\n")  
            OF.write('ymin=%s\n' %str(ylim[0]))
            OF.write('ymax=%s\n' %str(ylim[1]))

        if 'Contours' in plotdict and len(yshorts) == 2:
            # contour plot
            if xlim is None:
                OF.write('xmin=min(data[:,0])\n')
                OF.write('xmax=max(data[:,0])\n')
            if ylim is None:
                OF.write('ymin=min(data[:,1])\n')
                OF.write('ymax=max(data[:,1])\n')
            
            OF.write('from scipy.interpolate import LinearNDInterpolator\n')
            OF.write('XY=np.column_stack((data[:,0],data[:,1]))\n')
            OF.write('interp=LinearNDInterpolator(XY,data[:,2])\n')
            OF.write('xI = np.linspace(xmin,xmax,1001)\n')
            OF.write('yI = np.linspace(ymin,ymax,1001)\n')
            OF.write('XI, YI = np.meshgrid(xI,yI)\n')
            OF.write('ZI =interp(XI,YI)\n')
            OF.write('plt.contour(XI,YI,ZI, levels='+str(plotdict["Contours"])+')\n')
            OF.write("plt.scatter(data[:,0],data[:,1])\n")

        else:  ## scatter plot

            for yshort in range(len(yshorts)):
                OF.write("plt.scatter(data[:,0],data[:,%d])\n" % (yshort+1))
        


        
    
        
        OF.write("plt.savefig('"+namestub+".pdf')\n")
        OF.write("#plt.show()\n") ## after savefig, otherwise empty figure
        OF.write("plt.close()\n")


#def make_script_zslha(xlabel,ylabel,outfile,namestub,xlim,ylim,the_scan):
def make_script_zslha(plotstub,plot,the_scan):
    #print("Making script "+namestub+'.py')
    bsmartpath=os.path.join(the_scan.MAINpath,'bsmart')
    outfile=the_scan.RunManager.output_filename
    
    ###labels, values, ranges and shortcuts (for read_small, optional)

    labels=plot['Labels']
    xlabel=labels[0]
    ylabel=labels[1]
    #plotstub=os.path.join(the_scan.plotdir,xlabel+'_'+ylabel)
    valuex=plot['Values'][0]
    valuey=plot['Values'][1]
    if 'Ranges' in plot:
        ranges=plot['Ranges']
        xlim=ranges[0]
        ylim=ranges[1]
    else:
        xlim=None
        ylim=None
          
    

    ## Need to get the info about each point


    with open(plotstub+'.py','w') as OF:
        OF.write('import matplotlib.pyplot as plt\n')
        OF.write('import matplotlib as mpl\n')
        OF.write('import numpy as np\n')
        OF.write('import sys\n')
        OF.write('import os\n\n')
        OF.write('sys.path.insert(0, "'+bsmartpath+'")\n')

        OF.write('from bsmart import zslha\n\n')
        
        if 'Shortcuts' in plot:
            OF.write("all_spcs=zslha.read_small('"+outfile+"',"+plot['Shortcuts']+",sep='ENDOFPARAMETERPOINT')\n")
        else:   
            OF.write("all_spcs=zslha.read('"+outfile+"','ENDOFPARAMETERPOINT')\n")
        #OF.write("data=pd.DataFrame(input_data,columns=['"+xlabel+"','"+ylabel+"']).to_numpy()\n")

        OF.write("datax = [float(spc.Value('"+valuex[0]+"',"+str(valuex[1])+")) for spc in all_spcs]\n")
        OF.write("datay = [float(spc.Value('"+valuey[0]+"',"+str(valuey[1])+")) for spc in all_spcs]\n")

        OF.write("fig=plt.figure(figsize=(7,6),constrained_layout='True')\n")
        OF.write("plt.grid('on')\n")

        if 'Scaling' in plot:
            xscale=plot['Scaling'][0]
            yscale=plot['Scaling'][1]
            OF.write('plt.xscale(\"'+xscale+'\")\n')
            OF.write('plt.yscale(\"'+yscale+'\")\n')

        OF.write("plt.scatter(datax,datay)\n")
        OF.write("plt.xlabel(r'"+xlabel+"',fontsize=13)\n")
        OF.write("plt.ylabel(r'"+ylabel+"',fontsize=13)\n")
        if xlim is not None:
            OF.write("plt.xlim("+str(xlim)+")\n")  
        if ylim is not None:
            OF.write("plt.ylim("+str(ylim)+")\n")  

        
        OF.write("plt.savefig('"+plotstub+".pdf')\n")
        OF.write("#plt.show()\n") ## after savefig, otherwise empty figure
        OF.write("plt.close()\n")
    ## Now run it ... 
    debug.command_line_log("python3 "+plotstub+".py",the_scan.log)
    """
    import subprocess
    command_line = subprocess.Popen("python3 "+plotstub+".py", shell=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
    out, err = command_line.communicate()
    """

def getfrominput(no):
    tvar = no
    return lambda x: x.inputs[tvar]

def getfromoutput(no):
    tvar = no
    return lambda x: x.outputs[tvar]

    

def make_plots(the_scan):

    ## Choose strategy
    strategy='Valid'
    
    if 'Strategy' in the_scan.inputs['Plotting']:
        strategy=the_scan.inputs['Plotting']['Strategy']
        
    if strategy == 'csv':

        """
        Format as dicts:
        "Plots": {
            "myplot1": {"Labels": ["$m_0$ (GeV)","$m_h$ (GeV)"], "Values": ["m0","mh"]}
        }
        """

        try:
            import pandas as pd
        except:
            raise NameError('Failed to import pandas')

        csvfile=the_scan.RunManager.csv_filename

        if not os.path.isfile(csvfile):
            the_scan.log.info('No csv file produced from scan, presumably no valid points')
            return
        #print(csvfile)
        #input_data=pd.read_csv(csvfile)
        allshorts=the_scan.inputs['Variables'].keys()
        label_dict = {k: None for k in allshorts}
        for vv in allshorts:
            if 'Label' in the_scan.inputs['Variables'][vv]:
                label_dict[vv] =  the_scan.inputs['Variables'][vv]["Label"]
        #for automatic corner plot, need the list of variables and labels for them


        for plotname,plotdict in the_scan.plotting['Plots'].items(): 

            
            
            x =plotdict['Values'][0]
            y =plotdict['Values'][1]

            ys = plotdict['Values'][1:]
            if 'Ranges' in plotdict:
                xlim=plotdict['Ranges'][0]
                ylim=plotdict['Ranges'][1]
            else:
                xlim=None
                ylim=None
            if 'Labels' in plotdict:
                xlabel=plotdict['Labels'][0]
                ylabel=plotdict['Labels'][1]
                if x in allshorts and label_dict[x] is None:
                    label_dict[x] = xlabel
                if y in allshorts and label_dict[y] is None:
                    label_dict[y] = ylabel

            else:
                xlabel=x
                ylabel=y
             
            
            plotstub=os.path.join(the_scan.plotdir,plotname)
            the_scan.log.info('Writing script for plot '+plotstub+'.py')
            make_script_csv(x,ys,xlabel,ylabel,csvfile,plotstub,xlim,ylim,plotdict)

            the_scan.log.info('Making plot '+plotstub+'.pdf')
            
            debug.command_line_log('python3 '+plotstub+'.py',the_scan.log)
            #os.system('python3 '+plotstub+'.py')
            #data=pd.DataFrame(input_data,columns=[x,y]).to_numpy()
            #make_plot(xlabel,ylabel,data[:,0],data[:,1],plotstub+'.pdf',xlim,ylim,plotdict)
        # now make automatic corner plot
        for x in label_dict:
            if label_dict[x] is None:
                label_dict[x] = x
        #make_auto_corner_csv(namestub,csvfile,shorts,plotdict)
        plotstub=os.path.join(the_scan.plotdir,'auto_corner')
        make_auto_corner_csv(plotstub,csvfile,allshorts,label_dict)
        the_scan.log.info('Making plot auto_corner.png')
            
        debug.command_line_log('python3 '+plotstub+'.py',the_scan.log)
        plotstub=os.path.join(the_scan.plotdir,'fancy_auto_corner')
        make_fancy_auto_corner_csv(plotstub,csvfile,allshorts,label_dict)
        the_scan.log.info('Making plot fancy_auto_corner.pdf')
            
        debug.command_line_log('python3 '+plotstub+'.py',the_scan.log)


    elif strategy == 'SLHA': ### needs canny use of zslha
        #outfile=the_scan.RunManager.output_filename
        
        """
        ## Here we will expect the plots to be listed as a dict, with labels, values, ranges and shortcuts (for read_small)
        Format, e.g. 
        "Plots":
            "myplot1" : { "Labels": ["m0","mh"], "Values": [["MINPAR",[1]],["MASS",[25]]], "Ranges":[[100,3000],[118,125]]}},
            "myplot2" : ...
        ## It would be more efficient to make all plots in one file, but we want to give flexibility for the user to modify the script
        """
        if not os.path.isfile(the_scan.RunManager.output_filename):
            the_scan.log.info('No output file produced from scan, presumably no valid points')
            return
        for plot in the_scan.plotting['Plots']:
            plotstub=os.path.join(the_scan.plotdir,plot)
            make_script_zslha(plotstub,the_scan.plotting['Plots'][plot],the_scan)
      
    #elif strategy == 'Results':
    #    return
    else: ## take only valid points


        for plotname,plotdict in the_scan.plotting['Plots'].items(): 
            x =plotdict['Values'][0]
            y =plotdict['Values'][1]
            if 'Ranges' in plotdict:
                xlim=plotdict['Ranges'][0]
                ylim=plotdict['Ranges'][1]
            else:
                xlim=None
                ylim=None
            if 'Labels' in plotdict:
                xlabel=plotdict['Labels'][0]
                ylabel=plotdict['Labels'][1]
            else:
                xlabel=x
                ylabel=y


            datax=[]
            datay=[]
            
            if x in the_scan.inputs['Variables']:
                for n,lab in enumerate(the_scan.inputs['Variables']):
                    if lab==x:
                        
                        def getx(point):
                            return float(point.inputs[n])
                        break
            elif x in the_scan.inputs['Observables']: # and (strategy != 'All'):
                
                for n,lab in enumerate(the_scan.inputs['Observables']):
                    if lab==x:
                        def getx(point):
                            return float(point.outputs[n])
                        break
            if y in the_scan.inputs['Variables']:
                for n,lab in enumerate(the_scan.inputs['Variables']):
                    if lab==y:
                        #print('%d, %s found in Variables: ' % (n,lab))
                        def gety(point):
                            return float(point.inputs[n])
                        break
            elif y in the_scan.inputs['Observables']:
                for n,lab in enumerate(the_scan.inputs['Observables']):
                    if lab==y:
                        #print('%d, %s found in Observables ' % (n,lab))
                        #gety =getfromoutput(nnn)
                        def gety(point):
                            return float(point.observables[n])
                        break      

            if strategy=='All':
                for point in the_scan.RunManager.all_points:
                    #print(point.inputs)
                    datax.append(getx(point))
                    datay.append(gety(point))
            else:
                for point in the_scan.RunManager.valid_points:
                    datax.append(getx(point))
                    datay.append(gety(point))

            make_plot(xlabel,ylabel,datax,datay,os.path.join(the_scan.plotdir,plotname+'.pdf'),xlim,ylim,plotdict)
