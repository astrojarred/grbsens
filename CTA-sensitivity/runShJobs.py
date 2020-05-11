#!/usr/bin/env python
# Script name: runShJobs.py
# This program creates the sh files to run jobs (it refers to the script sensi.py)
# Author: B. Patricelli (barbara.patricelli@pi.infn.it)

# Import
import os,sys, random, math, glob
from optparse import OptionParser


if __name__ == '__main__':
    
    usg = "\033[1;31m%prog [ options] mainConfig.txt\033[1;m \n"
        
    desc = "\033[34mThis is the script to submit jobs\033[0m \n"
    
    
    # Input parameters (other parameters can be added, accordingly to sensi.py)
    parser=OptionParser(description=desc,usage=usg)
    parser.add_option("-r","--numruns",type='int',default=1,help="Total Number of runs")
    parser.add_option("-j","--numjobs",type='int',default=1,help="Total Number of jobs")
    
    parser.add_option("-l","--loweren",type='float',default=0.02,help="lower energy limit (TeV)")
    parser.add_option("-u","--upperen",type='float',default=10.0,help="upper energy limit (TeV)")
    parser.add_option("-b","--bins",type='int',default=1,help="number of energy bins")
    parser.add_option("-i","--irf",type='str',default='North_0.5h',help="IRF")
    parser.add_option("-t","--inTime",type='int',default=1,help='initial obs time')
    parser.add_option("-d","--deltat",type='int',default=1,help='delta t')    
    parser.add_option("-s","--sigma",type='float',default=5.,help='significance')
    parser.add_option("-o", "--offset",type='float',default=0.,help='offset from the pointing coordinates')
                       
    (options,args) = parser.parse_args()
    if (len(args)!=0):
        parser.error("incorrect number of arguments. Use -h for help")
                                    
    NumberRuns     = options.numruns
    NumberJobs     = options.numjobs

    emin = options.loweren
    emax = options.upperen
    bins = options.bins
    irf = options.irf
    t=options.inTime
    dt = options.deltat
    sigma=options.sigma
    offset=options.offset

    #List parameters and args                                                                                                                                                  
    print '\n**Input parameters:'
    for key, val in parser.values.__dict__.iteritems():
        print key,':', val
        if (val==None):
            print '\nERROR!',key,'is a required option! Exit.\nType ',ScriptName+'.py -h for help\n'
            sys.exit(1)



    STARTDIR=os.getcwd()

    # Creating a directory for all the sh files
    softwaredir=STARTDIR+'/software/'
    cmd0='mkdir '+softwaredir
    os.system(cmd0)

#####################################################
    for i in range(1,NumberJobs+1):
        t=t+dt

        print 'Processing Job',i,'/',(NumberJobs) 
        runfileName="runjob-%.d.sh"%(i)
        rfileName=softwaredir+runfileName
        print 'Preparing run file',runfileName
        runfile = open(rfileName,"w")

        runfile.write("#!/bin/env bash \n")      
        runfile.write("STARTDIR=%s\n"%(STARTDIR))
        runfile.write("cd $STARTDIR \n")

        cmd1 = 'python $STARTDIR/sensi.py -d '+str(t)+' -o '+str(offset)+' -i '+str(irf)+' -l '+str(emin)+' -u '+str(emax)+' -s '+str(sigma) +' \n'  
        runfile.write(cmd1)
        runfile.close()

        cmd5='sh '+str(rfileName)+' &'
        os.system(cmd5)

