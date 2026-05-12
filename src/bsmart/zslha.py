import subprocess
import os
from six import string_types
from collections import OrderedDict 

import math

import copy # for copying spc
from typing import Optional # for the declaration of __init__ below
"""
----------------------------------------------------------
SLHA parser for use with BSMArt
By Mark Goodsell (goodsell@lpthe.jussieu.fr)

Initially based on xslha by Florian Staub (florian.staub@gmail.com)
----------------------------------------------------------

NB widths and branching ratios are stored with integer keys.
Each branching ratio is specified by a tuple of integers

XSECTION blocks are specified as a float, then series of integer tuples

The keys of all other blocks are strings.
"""




class SLHA():
    def __init__(self,spc: Optional['SLHA'] = None):
        """
            Initialise object
            If given an existing spectrum file, will deepcopy that
        """

        if spc is None:
            self.blocks = OrderedDict() #### Necessary to preseve ordering of the blocks
            self.blockcomments = OrderedDict() ### store comments after block entries for writing input files
            self.br = {}
            self.widths = {}
            self.br1L = {}
            self.widths1L = {}
            self.xsections = {}
            self.analyses = {}
            
        else:
            self.blocks=copy.deepcopy(spc.blocks)
            self.blockcomments=copy.deepcopy(spc.blockcomments)
            self.br = copy.deepcopy(spc.br)
            self.widths = copy.deepcopy(spc.widths)
            self.br1L = copy.deepcopy(spc.br1L)
            self.widths1L = copy.deepcopy(spc.widths1L)
            self.xsections = copy.deepcopy(spc.xsections)
            self.analyses = copy.deepcopy(spc.analyses)
        
        # Objects necessary for processing the reading of a file
        self.block_name = None
        self.entries = OrderedDict() #### Necessary to preseve ordering of the entries
        self.comments = OrderedDict() ### store comments after block entries for writing input files
        
        
        
        self.reading_block = False
        self.reading_decay = False
        self.reading_xsection = False
        self.reading_hb_fermion = False
        self.reading_hb_boson = False
        self.decay1L = False
        self.decay_part = 0



    def add(self,extra_spc):
        self.blocks.update(extra_spc.blocks)
        self.br.update(extra_spc.br)
        self.br1L.update(extra_spc.br1L)
        self.widths.update(extra_spc.widths)
        self.widths1L.update(extra_spc.widths1L)
        self.blockcomments.update(extra_spc.blockcomments)
        self.xsections.update(extra_spc.xsections)
        self.analyses.update(extra_spc.analyses)

    # return wdith and BR
    def BR(self, init, final):
        # frozenset: make sure that the final states are order-less
        return self.br[init][tuple(sorted(final))]

    def Width(self, pdg): 
        return self.widths[pdg]

    def Value(self, block, number):
        '''return value of a parameter defined by block and entry
            or the width or an BR'''
        if block == 'WIDTH':
            return self.widths[number]
        elif block == 'BR':
            return self.br[number[0]][tuple(sorted(number[1]))]
        elif block == 'WIDTH1L':
            return self.widths1L[number]
        elif block == 'BR1L':
            return self.br1L[number[0]][tuple(sorted(number[1]))]
        elif block == 'XSECTION':
            # Want to return an actual value rather than a list -> just take the first (who computes more than one anyway??!)
            key=tuple([float(number[0]),tuple([int(number[1]),int(number[2])]),tuple(map(int,number[3:]))])
            xsecs=self.xsections[key]
            xs=next(iter(xsecs.values()))
            return float(xs)
            #xs = self.xsections[tuple(number)]
            #return [[x, xs[x]] for x in xs.keys()]
        elif block == 'PROCESS':
            analysis=number[0].upper()
            return self.analyses[analysis]['PROCESS'][number[1]]
        elif block == 'EFFICIENCIES':
            analysis=number[0].upper()
            return self.analyses[analysis]['EFFICIENCIES'][number[1]][int(number[2])]
        elif (block == 'HIGGSCOUPLINGSFERMIONS' or block == 'HIGGSCOUPLINGSBOSONS') and len(number) ==2:
            res=self.blocks[block.upper()][str(number[0])[1:-1].replace(" ", "")]
            return res[int(number[1])]
        else:
            return self.blocks[block.upper()][
                str(number)[1:-1].replace(" ", "")]

    def safeValue(self, block, number,returnval=0.0):
        try:
            return self.Value(block,number)
        except:
            return returnval
    def safeValue2(self, block, number,returnval=[0.0,0.0]):
        try:
            return self.Value(block,number)
        except:
            return returnval
            
    def start_decay(self, li):
        parsed = list(filter(None, li.split(' ')))
        
        self.decay1L = li.upper().startswith("DECAY1L")
        self.decay_part = int(parsed[1])

        if len(parsed) > 2:
            if self.decay1L:
                try:
                    self.widths1L[self.decay_part] = float(parsed[2])
                except: ## should only happen if we have an expression here
                    self.widths1L[self.decay_part] = parsed[2]
                
            else:
                try:
                    self.widths[self.decay_part] = float(parsed[2])
                except:
                    self.widths[self.decay_part] = parsed[2]
        
        self.entries = OrderedDict()
        self.reading_block, self.reading_decay, self.reading_xsection \
            = False, True, False

    def start_block(self, li):
        splitline=list(filter(None, li.split(' ')))
        self.block_name = splitline[1].upper()
        self.entries = OrderedDict()
        self.comments = OrderedDict()
        self.reading_block, self.reading_decay, self.reading_xsection \
            = True, False, False
        self.reading_hb_boson = \
            self.block_name in ["HIGGSBOUNDSINPUTHIGGSCOUPLINGSBOSONS",
                                "HIGGSCOUPLINGSBOSONS"]
        self.reading_hb_fermion = \
            self.block_name in ["HIGGSBOUNDSINPUTHIGGSCOUPLINGSFERMIONS",
                                "HIGGSCOUPLINGSFERMIONS"]
        if self.block_name == 'PROCESS':
            
            self.reading_process = True
            self.current_analysis=splitline[2].upper()
            if not self.current_analysis in self.analyses:
                self.analyses[self.current_analysis]={}
            self.analyses[self.current_analysis]['PROCESS']=OrderedDict()
        else:
            self.reading_process = False
        if self.block_name == 'EFFICIENCIES':
            self.reading_efficiency = True
            self.current_analysis=splitline[2].upper()
            if not self.current_analysis in self.analyses:
                self.analyses[self.current_analysis]={}
            self.analyses[self.current_analysis]['EFFICIENCIES']=OrderedDict()
        else:
            self.reading_efficiency = False
            

    def start_xsection(self, li):
        parsed = list(filter(None, li.split(' ')))
        if "#" in parsed:
            parsed = parsed[:parsed.index("#")]  # remove comments
        ## format is SQRTS PDGIN1 PDGIN2 NFINAL PDGFINAL1 ...
        self.xs_head = tuple(
            [float(parsed[1]),
             tuple([int(parsed[2]), int(parsed[3])]),
             tuple(map(int,parsed[5:]))
             ])

            
        # self.xs_head = tuple(
        #     [float(parsed[1]),
        #      tuple([int(parsed[2]), int(parsed[3])]),
        #      tuple([int(parsed[-2]), int(parsed[-1])])
        #      ])
        self.entries = OrderedDict() ## don't care about ordering here 
        self.reading_block, self.reading_decay, self.reading_xsection \
            = False, False, True

    def flush(self):
        '''store the information once a block is completely parsed'''
        if len(self.entries) > 0:
            if self.reading_block:
                self.blocks[self.block_name] = self.entries
                self.blockcomments[self.block_name] = self.comments
                self.comments=OrderedDict()
            if self.reading_decay:
                if self.decay1L:
                    self.br1L[self.decay_part] = self.entries
                else:
                    self.br[self.decay_part] = self.entries
            if self.reading_xsection:
                self.xsections[self.xs_head] = self.entries
                
    def write(self,
              filename: str,
              append: bool = False
              ):
        """
            Write all data to an SLHA file. Here we do not attempt to evaluate strings
            Takes a filename as argument
            By default will overwrite
            Optional argument to append to the file instead
        """
        
        if append:
            fstr='a'
        else:
            fstr='w'

        with open(filename, fstr) as f:
                
            # First do blocks
            tmaxlen=max(map(len,self.blocks.keys()))+3 # get the maximum blockname length to have uniform sizes
            for bname,bdata in self.blocks.items():
                f.write("BLOCK "+bname.upper()+" #\n")
                comments = self.blockcomments.get(bname,{})
                formatter = ' %%-%ds %%14.8e # %%s\n' % tmaxlen
                for v in bdata:
                    f.write(formatter % (v.replace(","," "), float(bdata[v]), comments.get(v,'')))
                
                
            # write decay info
            for d in self.widths:
                if isinstance(self.widths[d], string_types):
                        twidth=str(self.widths[d])
                else:
                    twidth='%14.8E' % self.widths[d]
                f.write("DECAY %d  %s  # \n" %(d,twidth))    
                
                
                if d in self.br:
                    for decay_mode in self.br[d]:
                        formatter='     %14.8E   %d  '% (self.br[d][decay_mode],len(decay_mode))
                        formatter=formatter + ('%d  '*len(decay_mode))
                        outline=formatter % decay_mode
                        outline=outline + '# \n'
                        f.write(outline)
            # write xsection info
            for xs_head in self.xsections:
                #sqrts=xs_head[0]
                #pdgin1,pdgin2=xs_head[1],xs_head[2]

                line1='XSECTION %.4e %d %d ' %(float(xs_head[0]),int(xs_head[1][0]),int(xs_head[1][0]))
                products=list(xs_head[2])
                line1=line1+'%d ' % len(products)
                for prod in products:
                    line1=line1+'%d ' % prod
                line1=line1+' # \n'
                f.write(line1)
                for xs_line,val in self.xsections[xs_head].items():
                    #print(xs_line)
                    outline=' '+' '.join(list(map(str,xs_line)))
                    if isinstance(val, string_types):
                        outline=outline+' %s ' % val
                    else:
                        outline=outline+' %.8e ' %val
                    #print(outline)
                    f.write(outline+' # \n')


        '''
         with open(filename, 'w+') as f:
             for b in self.blocks:
                 write_block_head(b, f)
                 write_block_entries(self.blocks[b], f)           
         '''       

    def to_dict(self):
        return {
            'blocks' : self.blocks,
            'blockcomments' : self.blockcomments,
            'br' : self.br,
            'widths' : self.widths,
            'br1L' : self.br1L,
            'widths1L' : self.widths1L,
            'xsections' : self.xsections,
            'analyses' : self.analyses
        } # since python 3.7 we preserve ordering so this should be fine

    def from_dict(self,indict):
        self.blocks= indict.get('blocks',OrderedDict()),
        self.blockcomments = indict.get(['blockcomments'],OrderedDict())
        self.br = indict.get('br', {})
        self.widths = indict.get('widths', {})
        self.br1L = indict.get('br1L', {})
        self.widths1L = indict.get('widths1L', {})
        self.xsections = indict.get('xsections', {})
        self.analyses = indict.get('analyses', {})


    #def copy(self): # returns a deep copy of the spectrum file
        
# ----------------------------------------------------------
# Reading
# ----------------------------------------------------------


# now the main function to read the SLHA file
def read(file, separator=None, verbose=False,firstblock=None): 
    """
    file: SLHA file
    separator: e.g. ENDOFPARAMETERPOINT, in order to have multiple spectrum files concatenated into one
    verbose: for extra info
    firstblock: NEW, for reading only everything after a certain input, when e.g. it has been 
                appended to a file (e.g. HiggsBounds, Vevacious)
    """
    spc = SLHA()
    if separator is not None:
        all_files = []
        count = 1
    found_firstblock=False
    with open(file) as infile:
        for line in infile:
            li = line.strip()
            liu = li.upper()

            if li.startswith("#") or len(li) < 1:
                continue
            if firstblock is not None and not found_firstblock:
                if liu.startswith(firstblock):
                    found_firstblock=True
                else:   
                    continue

            if separator is not None:
                if liu.startswith(separator):
                    found_firstblock=False
                    spc.flush()
                    if max(len(spc.blocks.keys()),len(spc.widths.keys())) > 0:
                        all_files.append(spc)
                        # start next point
                        spc = SLHA()
                        count = count + 1
                        if verbose:
                            print("Read spc file:", count)
                        continue

            # New block started
            
            if liu.startswith("BLOCK"):
                spc.flush()  # store information which was read
                spc.start_block(li)
            elif liu.startswith("DECAY"):
                spc.flush()  # store information which was read
                spc.start_decay(li)
            elif liu.startswith("XSECTION"):
                spc.flush()  # store information which was read
                spc.start_xsection(li)
            

            # Reading and parsing values
            else:
                li = li.replace('\t',' ')
                # remove comments first
                if "#" in li: # we know that it is not at the beginning, because we checked
                    firstcomment=li.find('#')

                    comment = li[firstcomment+1:]
                    li = li[:firstcomment]
                else:
                    comment = ''

                parsed = list(filter(None, li.split(' ')))

                

                """
                parsed = list(filter(None, li.split(' ')))
                if "#" in parsed:
#                    comment = ' '.join(parsed[parsed.index("#"):]) # store comments
                    try:
                        comment=line[line.find('#')+1:].strip()
                    except:
                        comment=''
                    parsed = parsed[:parsed.index("#")]  # remove comments, NB this doesn't work if the line contains "##" etc
                else:
                    comment = ''
#                    comment = '#'    ## ensure that each line does have a comment ...
                """
                if spc.reading_block:
                    
                    
                    if spc.reading_hb_fermion:
                        spc.entries[",".join(parsed[3:])] = \
                            [float(parsed[0]), float(parsed[1])]
                    elif spc.reading_hb_boson:
                        spc.entries[",".join(parsed[2:])] = \
                            float(parsed[0])
                    elif spc.reading_process:
                        spc.analyses[spc.current_analysis]['PROCESS'][parsed[0]] = float(parsed[1])   ### For HackAnalysis
                    elif spc.reading_efficiency:
                        spc.analyses[spc.current_analysis]['EFFICIENCIES'][parsed[0]] = [float(parsed[1]),float(parsed[2])] ### For HackAnalysis
                    else:
                        # Value might be a string like in SPINFO block. But these should just be interpreted like any other
                        #value=parsed[-1]
                        if spc.block_name.upper() == "SPINFO":
                            spc.entries[parsed[0]] = ' '.join(parsed[1:])
                        else:
                            value=parsed[-1]

                                
                            spc.entries[",".join(parsed[0:-1])] = value
                            if comment:
                                spc.comments[",".join(parsed[0:-1])] = comment

                        
                if spc.reading_decay:
                    try:
                        spc.entries[
                            tuple(sorted(eval("[" + ",".join(parsed[2:]) + "]")))
                        ] = float(parsed[0]) ##  NB could change this to allow user to specify expression in the branching ratio
                    except:
                        spc.entries[
                            tuple(sorted(eval("[" + ",".join(parsed[2:]) + "]")))
                        ] = parsed[0]

                if spc.reading_xsection:
                    """ Format is SCALE_SCHEME QCD_ORDER EW_ORDER KAPPA_F KAPPA_R  PDF_ID VALUE CODE VERSION """
                    try:
                        spc.entries[
                            tuple(eval("[" + ",".join(parsed[0:6]) + "]"))
                        ] = float(parsed[6])  ##  NB could change this to allow user to specify expression in the branching ratio
                    except:
                        spc.entries[
                            tuple(eval("[" + ",".join(parsed[0:6]) + "]"))
                        ] = parsed[6] 
                    #spc.entries[
                    #    tuple(eval("[" + ",".join(parsed[0:-2]) + "]"))
                    #] = float(parsed[-2])

    spc.flush()  # save the very last block in the file

    if verbose:
        print("Read %i blocks and %i decays" % (len(spc.blocks), len(spc.br)))
    if separator is None:
        return spc
    else:
        if len(spc.entries) > 0:
            all_files.append(spc)
        return all_files


# wrapper for faster read-in of multiple files
# squeeze the file (just keeping the necessary entries) to make the reading more efficient
# example: read_small_spc(filename,["# m0","# m12","# relic"],separator="ENDOFPARAMETERPOINT")
def read_small(file, entries, sep):
    if entries is None:
        out = read(file, separator=sep)
    else:
        string = "--regexp=\"" + sep + "\" --regexp=\"Block\" "
        for i in entries:
            string = string + "--regexp=\"" + i + "\" "
        if os.path.isfile("temp.spc"):
            subprocess.call("rm temp.spc", shell=True)
        subprocess.call("cat " + file + " | grep -i " + string
                        + " > temp_read_small.spc", shell=True)
        out = read("temp_read_small.spc", separator=sep)
        subprocess.call("rm temp_read_small.spc", shell=True)

    return out


def read_dir(dir, entries=None):
    if os.path.isfile("temp_read_dir.spc"):
        subprocess.call("rm temp_read_dir.spc", shell=True)
#    subprocess.check_call("cat "+dir+"/* > temp_read_dir.spc",shell=True)
    subprocess.check_call("tail -n+1 " + dir + "/* > temp_read_dir.spc",
                          shell=True)
    out = read_small("temp_read_dir.spc", entries, "==>")
    subprocess.call("rm temp_read_dir.spc", shell=True)

    return out


# ----------------------------------------------------------
# Writing
# ----------------------------------------------------------

###### Point can be specified in terms of values and their names as two lists
# So you can write the values in the file as e.g m0, m12 etc and they will be converted.
# or you can give just the values, in which case they should be present in the slha file as VARIABLE[0], VARIABLE[1], ...

def write_lh(spc, var_vals, filename, var_names = None):
    if var_names is not None:
        var_dict = { var_name : var_val for var_name,var_val in zip(var_names,var_vals)}
        var_dict['VARIABLE'] = var_vals  ## Added for backward compatibility, maybe remove?
    else:
        var_dict={'VARIABLE':var_vals}   ## Added for backward compatibility, maybe remove?
    var_dict['math'] = math ## to allow maths functions such as sqrt or cos/sin/tan!

    write_lh_dict(spc,var_dict, filename)

def write_lh_dict(spc, var_dict, filename):
    
     
    with open(filename, 'w+') as f:
        for b in spc.blocks:
            f.write("BLOCK "+b.upper()+" #\n")
            #write_block_numbers(b, spc.blocks[b], var_vals, f, comments = spc.blockcomments[b])
            write_block_numbers(b, spc.blocks[b], var_dict, f, comments = spc.blockcomments[b])
        # write decay info
        for d in spc.widths:
            if isinstance(spc.widths[d], string_types):
                try:
                    #twidth='%10.6E' % eval(spc.widths[d],{'VARIABLE':variables})
                    twidth='%14.8E' % eval(spc.widths[d],var_dict)
                
                except:
                    twidth=str(spc.widths[d])
            else:
                twidth='%14.8E' % spc.widths[d]
            f.write("DECAY %d  %s  # \n" %(d,twidth))    
            
            
            if d in spc.br:
                for decay_mode in spc.br[d]:
                    formatter='     %14.8E   %d  '% (spc.br[d][decay_mode],len(decay_mode))
                    formatter=formatter + ('%d  '*len(decay_mode))
                    outline=formatter % decay_mode
                    outline=outline + '# \n'
                    f.write(outline)
        # write xsection info, for the few cases we might need it. For now I don't assume the user will change things much
        for xs_head in spc.xsections:
            #sqrts=xs_head[0]
            #pdgin1,pdgin2=xs_head[1],xs_head[2]

            line1='XSECTION %.4e %d %d ' %(float(xs_head[0]),int(xs_head[1][0]),int(xs_head[1][0]))
            products=list(xs_head[2])
            line1=line1+'%d ' % len(products)
            for prod in products:
                line1=line1+'%d ' % prod
            line1=line1+' # \n'
            f.write(line1)
            for xs_line,val in spc.xsections[xs_head].items():
                #print(xs_line)
                outline=' '+' '.join(list(map(str,xs_line)))
                if isinstance(val, string_types):
                    outline=outline+' %.8e ' %eval(val,var_dict)
                else:
                    outline=outline+' %.8e ' %val
                #print(outline)
                f.write(outline+' # \n')
def write(blocks, filename):
    with open(filename, 'w+') as f:
        for b in blocks:
            write_block_head(b, f)
            write_block_entries(blocks[b], f)


def write_block_entries(values, file):
    tmaxlen=max(map(len,values.keys()))+3
    formatter = ' %%-%ds %%14.8e # \n' % tmaxlen
    for v in values.keys():
        file.write(formatter % (v.replace(","," "), float(values[v])))


def write_les_houches(block, values, point, file):
    write_block_head(block, file)
    write_block_numbers(block, values, point, file)


def write_block_head(name, file):
    file.write("Block " + name.upper() + " # \n")


#def write_block_numbers(name, values, variables, file, comments = OrderedDict()):
def write_block_numbers(name, values, var_dict, file, comments = OrderedDict()):
    tmaxlen=max(map(len,values.keys()))+3
    for v in values.keys():
        try:
            comment=comments[v]
        except:
            comment=name.upper() + "[" + str(v) + "]"

        if isinstance(values[v], string_types):  # to be 2 and 3 compatible
            
            try:
                #teststr=str(eval(values[v],{'VARIABLE':variables}))
                teststr=str(eval(values[v],var_dict))
            except:
                teststr=values[v]
                
            if teststr==values[v]:
                
                formatter = ' %%-%ds %%s # %%s \n' % tmaxlen                
                file.write(formatter
                            % (v.replace(","," "), values[v],
                               comment))

                
            else:
                formatter = ' %%-%ds %% 14.8E # %%s \n' % tmaxlen     ##### space to allow spaces for minus sign
                #file.write(formatter % (v.replace(","," "), float(eval(values[v],{'VARIABLE':variables})),comment))
                file.write(formatter % (v.replace(","," "), float(eval(values[v],var_dict)),comment))
#                file.write(formatter % (v.replace(","," "), float(eval(values[v],{'VARIABLE':variables})),name.upper() + "[" + str(v) + "]"))

           ## now we will assume strings have been caught above, so now compare the string versions of the integer form : nb strings would throw an error here
#        elif isinstance(values[v], int):
        elif len(str(int(values[v]))) == len(str(values[v])):
            formatter = ' %%-%ds %% -14i # %%s \n' % tmaxlen   ##### space to allow spaces for minus sign, padded
            file.write(formatter
                       % (v.replace(","," "), values[v], comment))
#                       % (v.replace(","," "), values[v], name.upper() + "[" + str(v) + "]"))
        else:
            formatter = ' %%-%ds %% 14.8E # %%s \n' % tmaxlen
            file.write(formatter
                       % (v.replace(","," "), float(values[v]),
                          comment))




