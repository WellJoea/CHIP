#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
***********************************************************
* @File    : BamStat.py
* @Author  : Zhou Wei                                     *
* @Date    : 2020/11/24 17:23:32                          *
* @E-mail  : welljoea@gmail.com                           *
* @Version : --                                           *
* You are using the program scripted by Zhou Wei.         *
* If you find some bugs, please                           *
* Please let me know and acknowledge in your publication. *
* Thank you!                                              *
* Best wishes!                                            *
***********************************************************
'''

import os

class STSTAT():
    def __init__(self):
        self.ST = 'a'
    def CMD(self, INbam, Outpre):
        pass

class GETSTAT():
    def __init__(self):
        self.ST = 'a'
    def CMD(self, INbam, Outpre):
        pass


import argparse
def Args():
    parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                prefix_chars='-+',
                conflict_handler='resolve',
                description="",
                epilog="")

    parser.add_argument('-V','--version',action ='version',
                version='EcDNA version 0.1')

    subparsers = parser.add_subparsers(dest="commands",
                    help='models help.')
    P_Common = subparsers.add_parser('Common',conflict_handler='resolve', #add_help=False,
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                    help='The common parameters used for other models.')
    P_Common.add_argument("-f", "--infile",type=str,
                    help='''the input file or input number split by ",".''')
    P_Common.add_argument("-i", "--indir",type=str,
                    help='''the input directory.''')
    P_Common.add_argument("-o", "--outdir",type=str,default=os.getcwd(),
                    help="output file dir, default=current dir.")
    P_Common.add_argument("-n", "--njob",type=int,default=5,
                    help="The maximum number of concurrently running jobs.")
    P_Common.add_argument("-bd", "--bamdir", type=str, default='02.MiniMap',
                    help="input bam directory for fetch ")
    P_Common.add_argument("-fd", "--fetchdir", type=str,  default='03.SoftMap',
                    help="out directory for fetch")
    P_Common.add_argument("-sd", "--searchdir", type=str, default='03.SoftMap',
                    help="out directory of search")
    P_Common.add_argument("-md", "--mergedir", type=str,  default='03.SoftMap',
                    help="out directory of merge")
    P_Common.add_argument("-rd", "--regiondir", type=str, default='04.EcRegion',
                    help="out directory for region")
    P_Common.add_argument("-cd", "--checkdir", type=str, default='05.CheakBP',
                    help="out directory for check breakpoint of  plasmid")
    P_Common.add_argument("-ch", "--Chrom", action='store_true', default=True,
                    help="only keep the specified chromosome: 1-22,X,Y,MT.")
    P_Common.add_argument("-bt", "--bedtools", type=str, default='/share/home/share/software/bedtools2/bin/',
                    help="bedtools path")
    P_Common.add_argument("-st", "--samtools", type=str, default='/share/home/share/software/samtools-1.10/bin/',
                help="samtools path")
    P_Common.add_argument("-gb", "--gtfbed", type=str, default='/share/home/share/Repository/GenomeDB/Reference/Homo_Sapiens/ENSEMBL/Homo_sapiens.GRCh38.100.gtf.gene.bed',
                    help="the gene bed file used for annotation of regions")
    P_Common.add_argument("-ap", "--annopeak", type=str, default='/share/home/share/software/Homer/bin/annotatePeaks.pl',
                    help="the annotatePeaks.pl script")
    P_Common.add_argument("-gt", "--gtf", type=str, default='/share/home/share/Repository/GenomeDB/Reference/Homo_Sapiens/ENSEMBL/Homo_sapiens.GRCh38.100.gtf',
                    help="the genome gtf file.")
    P_Common.add_argument("-cb", "--checkbed", type=str, default='/share/home/zhou_wei/Workspace/11Project/02Plasmid/01analysescript/uniqueovr/BEDUniq.region.txt',
                    help="the bed file of plasmid unique region.")

    P_Autopipe = subparsers.add_parser('Auto', conflict_handler='resolve', prefix_chars='-+',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                    parents=[P_Common, P_fetch, P_search, P_merge, P_region],
                    help='the auto-processing for all.')
    P_Autopipe.add_argument("+P", "++pipeline",nargs='+',
                    help="the auto-processing: Fetch, Search, Region.")
    P_Autopipe.add_argument('+M','++MODEL' , nargs='+', type=str, default=['Standard'],
                    help='''Chose more the one models from Standard, Fselect,Fitting and Predict used for DIY pipline.''')
    args  = parser.parse_args()
    return args

import logging
class DispatchingFormatter:
    def __init__(self, formatters, default_formatter):
        self._formatters = formatters
        self._default_formatter = default_formatter

    def format(self, record):
        formatter = self._formatters.get(record.name, self._default_formatter)
        return formatter.format(record)

class Logger:
    level_dict = {
        'NOTSET'  : logging.NOTSET,
        'DEBUG'   : logging.DEBUG,
        'INFO'    : logging.INFO,
        'WARNING' : logging.WARNING,
        'ERROR'   : logging.ERROR,
        'CRITICAL': logging.CRITICAL,
    }

    ChangeFrom = DispatchingFormatter(
            { 'c' : logging.Formatter( '[%(asctime)s] [%(levelname)-4s]: %(message)s', '%Y-%m-%d %H:%M:%S'),
              'p' : logging.Formatter( '[%(levelname)-4s]: %(message)s'),
              'n' : logging.Formatter( '%(message)s' ),
            }, 
            logging.Formatter('%(message)s')
     )

    def __init__(self, outpath, filemode='w',  clevel = 'INFO', Flevel = 'INFO'):

        logging.basicConfig(
            level    = Logger.level_dict[clevel] ,
            format   = '[%(asctime)s] [%(levelname)-4s]: %(message)s',
            datefmt  = '%Y-%m-%d %H:%M:%S',
            filename = None,
        )

        File = logging.FileHandler(outpath,  mode= filemode)
        File.setLevel(Logger.level_dict[Flevel])
        File.setFormatter(Logger.ChangeFrom)
        logging.getLogger().addHandler(File)

        self.R = logging
        self.C = logging.getLogger('c')
        self.P = logging.getLogger('p')
        self.N = logging.getLogger('n')
        self.CI = logging.getLogger('c').info
        self.NI = logging.getLogger('n').info
        self.CW = logging.getLogger('c').warning
        self.NW = logging.getLogger('n').warning

import os
import time
import traceback
def Commands():
    info ='''
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
    os.makedirs(args.outdir, exist_ok=True)
    Log = Logger( '%s/%s_log.log'%(args.outdir, args.commands) )

    Log.NI(info.strip())
    Log.NI("The argument you have set as follows:".center(59, '*'))

    _ml = max([ len(i) for i in vars(args).keys()] )
    for i,(k,v) in enumerate(vars(args).items(), start=1):
        Log.NI( '**{:0>{}}|{:<{}}: {}'.format(i,2, k,_ml, v) )
    Log.NI(59 * '*')

    try:
        ecDNA(args, Log).Pipeline()
        Log.CI('Success!!!')
    except Exception:
        Log.CW('Failed!!!')
        traceback.print_exc()
    finally:
        Log.CI('You can check your progress in log file.')
Commands()
