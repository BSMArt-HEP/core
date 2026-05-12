"""
Tool to run SPheno


"""

__meta__ = {
    "name": "SPheno",
    "requires": [],
    "settings": {
        "Command": "Path to SPheno",
        "InputFile": "Input SLHA file"
    }
}


import os
import shutil
from bsmart import debug
from bsmart import zslha
from bsmart.HEPRun import HEPRun, RunSettings, CoreRunner, HepTool, DataPoint


references=r'''@article{Porod:2003um,
    author = "Porod, Werner",
    title = "{SPheno, a program for calculating supersymmetric spectra, SUSY particle decays and SUSY particle production at e+ e- colliders}",
    eprint = "hep-ph/0301101",
    archivePrefix = "arXiv",
    reportNumber = "ZU-TH-01-03",
    doi = "10.1016/S0010-4655(03)00222-4",
    journal = "Comput. Phys. Commun.",
    volume = "153",
    pages = "275--315",
    year = "2003"
}

@article{Porod:2011nf,
    author = "Porod, W. and Staub, F.",
    title = "{SPheno 3.1: Extensions including flavour, CP-phases and models beyond the MSSM}",
    eprint = "1104.1573",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.cpc.2012.05.021",
    journal = "Comput. Phys. Commun.",
    volume = "183",
    pages = "2458--2469",
    year = "2012"
}

@article{Staub:2009bi,
    author = "Staub, Florian",
    title = "{From Superpotential to Model Files for FeynArts and CalcHep/CompHep}",
    eprint = "0909.2863",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.cpc.2010.01.011",
    journal = "Comput. Phys. Commun.",
    volume = "181",
    pages = "1077--1086",
    year = "2010"
}


@article{Staub:2010jh,
    author = "Staub, Florian",
    title = "{Automatic Calculation of supersymmetric Renormalization Group Equations and Self Energies}",
    eprint = "1002.0840",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1016/j.cpc.2010.11.030",
    journal = "Comput. Phys. Commun.",
    volume = "182",
    pages = "808--833",
    year = "2011"
}

@article{Staub:2012pb,
    author = "Staub, Florian",
    title = "{SARAH 3.2: Dirac Gauginos, UFO output, and more}",
    eprint = "1207.0906",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "BONN-TH-2012-17",
    doi = "10.1016/j.cpc.2013.02.019",
    journal = "Comput. Phys. Commun.",
    volume = "184",
    pages = "1792--1809",
    year = "2013"
}

@article{Staub:2013tta,
    author = "Staub, Florian",
    title = "{SARAH 4 : A tool for (not only SUSY) model builders}",
    eprint = "1309.7223",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "BONN-TH-2013-17",
    doi = "10.1016/j.cpc.2014.02.018",
    journal = "Comput. Phys. Commun.",
    volume = "185",
    pages = "1773--1790",
    year = "2014"
}

@article{Staub:2015kfa,
    author = "Staub, Florian",
    title = "{Exploring new models in all detail with SARAH}",
    eprint = "1503.04200",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "CERN-PH-TH-2015-051",
    doi = "10.1155/2015/840780",
    journal = "Adv. High Energy Phys.",
    volume = "2015",
    pages = "840780",
    year = "2015"
}

@article{Goodsell:2014bna,
    author = "Goodsell, Mark D. and Nickel, Kilian and Staub, Florian",
    title = "{Two-Loop Higgs mass calculations in supersymmetric models beyond the MSSM with SARAH and SPheno}",
    eprint = "1411.0675",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "BONN-TH-2014-14, CERN-PH-TH-2014-213",
    doi = "10.1140/epjc/s10052-014-3247-y",
    journal = "Eur. Phys. J. C",
    volume = "75",
    number = "1",
    pages = "32",
    year = "2015"
}


@article{Goodsell:2015ira,
    author = "Goodsell, M. and Nickel, K. and Staub, F.",
    title = "{Generic two-loop Higgs mass calculation from a diagrammatic approach}",
    eprint = "1503.03098",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "BONN-TH-2015-04, CERN-PH-TH-2015-044",
    doi = "10.1140/epjc/s10052-015-3494-6",
    journal = "Eur. Phys. J. C",
    volume = "75",
    number = "6",
    pages = "290",
    year = "2015"
}

@article{Goodsell:2017pdq,
    author = "Goodsell, Mark D. and Liebler, Stefan and Staub, Florian",
    title = "{Generic calculation of two-body partial decay widths at the full one-loop level}",
    eprint = "1703.09237",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "DESY-17-042, KA-TP-11-2017",
    doi = "10.1140/epjc/s10052-017-5259-x",
    journal = "Eur. Phys. J. C",
    volume = "77",
    number = "11",
    pages = "758",
    year = "2017"
}

@article{Braathen:2017izn,
    author = "Braathen, Johannes and Goodsell, Mark D. and Staub, Florian",
    title = "{Supersymmetric and non-supersymmetric models without catastrophic Goldstone bosons}",
    eprint = "1706.05372",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "KA-TP-23-2017",
    doi = "10.1140/epjc/s10052-017-5303-x",
    journal = "Eur. Phys. J. C",
    volume = "77",
    number = "11",
    pages = "757",
    year = "2017"
}

@article{Goodsell:2018tti,
    author = "Goodsell, Mark D. and Staub, Florian",
    title = "{Unitarity constraints on general scalar couplings with SARAH}",
    eprint = "1805.07306",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "KA-TP-10-2018",
    doi = "10.1140/epjc/s10052-018-6127-z",
    journal = "Eur. Phys. J. C",
    volume = "78",
    number = "8",
    pages = "649",
    year = "2018"
}
'''


class NewTool(HepTool):
    def __init__(self, name, settings,global_settings=None):
        HepTool.__init__(self, name, settings,global_settings)
        self.citations=references
        """ Check that we can find the command, otherwise not much point carrying on"""
        self.command=os.path.expanduser(self.settings['Command'])
        if not os.path.isfile(self.command):
            #fdir=os.path.dirname
            raise NameError("Cannot identify SPheno executable. Please give it (without additional command line parameters) under SPheno->settings->Command")
        
    def run(self, spc_file, temp_dir, log,data_point):

        ## Remove messages.out and SPheno.out if they exist, otherwise SPheno appends and they can become large
        if os.path.exists(os.path.join(temp_dir,'Messages.out')):
            os.remove(os.path.join(temp_dir,'Messages.out'))
        if os.path.exists(os.path.join(temp_dir,'SPheno.out')):
            os.remove(os.path.join(temp_dir,'SPheno.out'))
        ## default have a timeout for SPheno of 600s.
        debug.command_line_log(self.command+' '+self.settings['InputFile']+' '+spc_file,log,ctimeout=600)
        #debug.command_line_log(self.settings['Command']+' '+self.settings['InputFile']+' '+spc_file,log,ctimeout=60)
        if os.path.exists(spc_file):
            data_point.spc = zslha.read(spc_file)
        
