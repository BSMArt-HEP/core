"""
Simple routine to split a LHE file into pieces for parallel handling in HackAnalysis

"""

import gzip
import os


def splitlhe(n_cores,infile,split_lhe_name='Split/split'):

    splitdir=os.path.dirname(split_lhe_name)
    
    if splitdir != '':
        #if not os.path.exists('Split'):
        #    os.makedirs('Split')
        if not os.path.exists(splitdir):
            os.makedirs(splitdir)
    fhs=[]

    if infile[-3:] == '.gz':
        ending='.gz'
    else:
        ending=''

    try:
        for x in range(n_cores):
            #os.remove('Split/split_'+str(x)+'.lhe.gz')
            #os.remove(split_lhe_name+'_'+str(x)+'.lhe.gz')
            tfile = split_lhe_name+'_'+str(x)+'.lhe'+ending
            if os.path.exists(tfile):
                os.remove(tfile)
    except:
        pass

    try:
        for x in range(n_cores):
            #tfh=gzip.open('Split/split_'+str(x)+'.lhe.gz','wb')
            if ending == '.gz':
                tfh=gzip.open(split_lhe_name+'_'+str(x)+'.lhe.gz','wb')
            else:
                tfh=open(split_lhe_name+'_'+str(x)+'.lhe','w')
            fhs.append(tfh)
    except Exception as e:
        raise  NameError("Could not open split output files "+str(e))

    try:
        if ending=='.gz':

            with gzip.open(infile, 'rb') as f:
                donebanner=False
                done=False
                eventhandle=0
                for line in f:
                    sline=line.decode('ascii')
                    if not donebanner:
                        if '<event' in sline:
                            donebanner=True
                            fhs[eventhandle].write(line)
                        else:
                            for fh in fhs: ## now write the banner to all of the output files
                                fh.write(line)
                    else:
                        if not done:
                            if '</event>' in sline:
                                # finish old event and advance the counter
                                fhs[eventhandle].write(line)
                                eventhandle=(eventhandle+1) % n_cores
                            else:
                                if '</LesHouchesEvent' in sline:
                                    done=True
                                    for fh in fhs:
                                        fh.write(line)
                                        fh.close()
                                else:
                                    # regular line
                                    fhs[eventhandle].write(line)
        else:
            with open(infile, 'r') as f:
                donebanner=False
                done=False
                eventhandle=0
                for line in f:
                    sline=line
                    if not donebanner:
                        if '<event' in sline:
                            donebanner=True
                            fhs[eventhandle].write(line)
                        else:
                            for fh in fhs: ## now write the banner to all of the output files
                                fh.write(line)
                    else:
                        if not done:
                            if '</event>' in sline:
                                # finish old event and advance the counter
                                fhs[eventhandle].write(line)
                                eventhandle=(eventhandle+1) % n_cores
                            else:
                                if '</LesHouchesEvent' in sline:
                                    done=True
                                    for fh in fhs:
                                        fh.write(line)
                                        fh.close()
                                else:
                                    # regular line
                                    fhs[eventhandle].write(line)           
    except Exception as e:
        raise NameError("Could not perform lhe file split, "+str(e))


if __name__=='__main__': ## allow this script to be called as standalone, with arguments
    import argparse
    
    parser = argparse.ArgumentParser(description='Please give the name of the input file.')
    parser.add_argument("--ncores", help="Number of cores", type=int,
                        action="store")
    parser.add_argument('inputfile',
                        metavar='File', type=str, help='Input file name')
    parser.add_argument('--Output', type=str, help='Output file name',action="store",default='Split/split')
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        #print("Please give an input file")
        raise SystemExit

    if args.ncores:
        n_cores=args.ncores
    else:
        import multiprocessing as mp
        n_cores=mp.cpu_count()
    
    infile=args.inputfile
    

    splitlhe(n_cores,infile,split_lhe_name=args.Output)