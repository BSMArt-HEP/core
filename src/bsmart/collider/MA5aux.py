"""

Auxiliary routines for MadAnalysis

"""


import math


def ReadCutFlow(cffile):
    try:
        IF=open(cffile,'r')
    except:
        raise NameError('Could not open cutflow file '+cffile)
    myN0=-1
    myNf=-1
    myNN=-1
    efficiencies=[]
    alleffs=[]
    cutname=""
    IsInitial=False
    IsCounter=False
    for line in IF:
        #print(line)
        if "<InitialCounter>" in line:
            IsInitial = True
            continue
        elif "</InitialCounter>" in line:
            IsInitial = False
            continue
        elif "<Counter>" in line:
            IsCounter = True
            continue
        elif "</Counter>" in line:
            cutname=""
            myNf=-1
            IsCounter = False
            continue
        if IsInitial and "sum of weights" in line and not '^2' in line:
            #print('Found N0')
            myN0 = float(line.split()[0])+float(line.split()[1])
        if IsCounter and line.startswith('\"'):
            #print(line.strip().find('#'))
            cutname=line[:line.strip().find('#')]
        if IsCounter and "sum of weights" in line and not '^2' in line:
            myNf = float(line.split()[0])+float(line.split()[1])
            if myN0 <=0.0 :
                raise NameError('Found cutflow but no initial events')
            eff=myNf/myN0
            if(eff > 0):
                err=math.sqrt((1-eff)/(eff*myNN))
            else:
                err=0.0
            #efficiencies.append([cutname,eff,err])
            efficiencies=[cutname,eff,err] ### return just the last line
            alleffs.append(efficiencies)
        if IsInitial and "nentries" in line:
            myNN=int(line.split()[0])

    #print(efficiencies)
    #close(IF)
    IF.close()
    #if efficiencies[1]>0.0001 and efficiencies[2] > 0.15:
    #    print('Uncertainty %.1f %% on efficiency %.2e for cutflow %s!' % (efficiencies[2]*100.0,efficiencies[1],cffile))
    #return {"efficiency": efficiencies[1], "uncertainty": efficiencies[2], "Total Weights": myN0, "Events": myNN}
    lasteff={"efficiency": efficiencies[1], "uncertainty": efficiencies[2], "Total Weights": myN0, "Events": myNN}
    return alleffs,lasteff



if __name__ == "__main__":
    import sys
    import os

    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        print("Specify a cutflow filename!")
        raise SystemExit

    try:
        alleff,leff=ReadCutFlow(sys.argv[1])
    except Exception as e:
        print("Failed: "+str(e))

    print("=====================================")
    print("File: %s" % sys.argv[1])
    print("Events: %d" % leff["Events"])
    for cut in alleff:
        print("%s: %.4e (=/- %.4f rel. error)" %(cut[0],cut[1],cut[2]))