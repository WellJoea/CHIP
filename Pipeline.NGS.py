#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
***********************************************************
* Author  : Zhou Wei                                      *
* Date    : 2020/09/09 10:42:02                           *
* E-mail  : welljoea@gmail.com                            *
* Version : --                                            *
* You are using the program scripted by Zhou Wei.         *
* If you find some bugs, please                           *
* Please let me know and acknowledge in your publication. *
* Thank you!                                              *
* Best wishes!                                            *
***********************************************************
'''
import pandas as pd
import numpy as  np
import re
import os
import sys
import time
import traceback
from Logger  import Logger
from ArgsPipe import Args
'''
pd.set_option('display.max_rows', 100000)
pd.set_option('display.max_columns', 100000)
pd.set_option('display.width', 100000)
'''
from joblib import Parallel, delayed
class MakeFlow:
    def __init__(self, arg, log,  *array, **dicts):
        self.arg = arg
        self.log = log
        self.array  = array
        self.dicts  = dicts
        self.wdir = self.arg.outdir + '/WorkShell'
        os.makedirs(self.wdir, exist_ok=True)
        self.Pipeline={
            'atac': '/share/home/share/Pipeline/13SCATAC/ATACPipe',
            'ss2': '/share/home/share/Pipeline/12SCRNA/SmartSeq2Pipe',
            'plasmid' : '/share/home/share/Pipeline/01NGSDNA/PlasmidPipe',
        }

    def mkinfo(self, d ):
        infors= '''
SID={SID}
SIDs={Sampleid}
UIDs={Uniqueid}
Lanes={Lane}
Reps={Rep}
R1s={R1}
R2s={R2}
OUTs={Outdir}
Speciess={Species}
Modules={Module}
WORKFLOW_DIRs={WORKFLOW_DIR}'''.format(**d).strip()
        f=open('%s/%s.input'%(self.wdir, d['SID']), 'w')
        f.write(infors)
        f.close()

    def mkshell(self, d):
        pass

    def mkflow(self):
        infodf = pd.read_csv(self.arg.input, sep='\t')
        infodf['Outdir'] = infodf['Outdir'].fillna(self.arg.outdir)
        infodf['WORKFLOW_DIR'] = infodf.Module.map(self.Pipeline)
        infodf[['R1','R2']] = infodf['Fastq'].str.split('[;,]',expand=True)
        infodf = infodf.astype(str).fillna('NA')

        for _g,_d in infodf.groupby('Sampleid'):
            print( {**{'SID':_g}, **_d.T.apply(lambda x: '('+' '.join(x) + ')', axis=1).to_dict()} )

        Parallel( n_jobs=-1 )\
                ( delayed( self.mkinfo )({**{'SID':_g}, **_d.T.apply(lambda x: '('+' '.join(x) + ')', axis=1).to_dict()})
                    for _g, _d in infodf.groupby('Sampleid'))

    
class Pipeline():
    def __init__(self, arg, log,  *array, **dicts):
        self.arg = arg
        self.log = log
        self.array  = array
        self.dicts  = dicts

    def Pipe(self):
        if self.arg.commands in ['mkflow', 'Auto']:
            MakeFlow(self.arg, self.log).mkflow()

def Commands():
    info ='''
>^o^<
***********************************************************
* Author : Zhou Wei                                       *
* Date   : %s                       *
* E-mail : welljoea@gmail.com                             *
* You are using The scripted by Zhou Wei.                 *
* If you find some bugs, please email to me.              *
* Please let me know and acknowledge in your publication. *
* Sincerely                                               *
* Best wishes!                                            *
***********************************************************
'''%(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))

    args = Args()
    Log = Logger( '%s/%s_log.log'%(args.outdir, args.commands) )
    os.makedirs( os.path.dirname(args.outdir) , exist_ok=True)

    Log.NI(info.strip())
    Log.NI("The argument you have set as follows:".center(59, '*'))
    for i,k in enumerate(vars(args),start=1):
        Log.NI('**%s|%-13s: %s'%(str(i).zfill(2), k, str(getattr(args, k))) )
    Log.NI(59 * '*')

    try:
        Pipeline(args, Log).Pipe()
        Log.CI('Success!!!')
    except Exception:
        Log.CW('Failed!!!')
        traceback.print_exc()
    finally:
        Log.CI('You can check your progress in log file.')
Commands()