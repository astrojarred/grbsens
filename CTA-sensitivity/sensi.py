#!/usr/bin/env python
#!/home/users/user_name/miniconda2/envs/my_env/bin/python
# Script name: sensi.py
# Description: this script estimates the CTA sensitivity with cssens
# Author: B. Patricelli (barbara.patricelli@pi.infn.it)
#######################################################
# Import
import os, sys
from optparse import OptionParser

import gammalib
import ctools


#Get Script name
ScriptName = os.path.split(sys.argv[0])[1].split('.')[0]

#######################################################
# Main 
#######################################################

if __name__ == '__main__':
    
  usg = "\033[1;31m%prog [ options] inputFile\033[1;m \n"

  desc = "\033[34mThis is the script for the sensitivity estimate with cssens \033[0m \n"

  # Input parameters
  parser=OptionParser(description=desc,usage=usg)
  parser.add_option("-m","--InputModel",default="grb.xml",help="Source model (xml file)")
  parser.add_option("-d","--obstime",type='float',default=1.,help="observing time (s)")
  parser.add_option("-o","--offset",type='float',default=0.0,help="offset from the pointing coordinates")
  parser.add_option("-l","--loweren",type='float',default=0.03,help="lower energy limit (TeV)")
  parser.add_option("-u","--upperen",type='float',default=10.0,help="upper energy limit (TeV)")
  parser.add_option("-b","--bins",type='int',default=1,help="number of energy bins")
  parser.add_option("-c","--caldb",type='str',default='prod2',help="caldb")
  parser.add_option("-i","--irf",type='str',default='North_0.5h',help="IRF")
  parser.add_option("-s","--sigma",type='float',default=5.0,help="significance")
  parser.add_option("-r","--radius",type='float',default=2.25,help="radius of the region of interest (deg)")
  parser.add_option("-z","--binsz",type='float',default=0.2,help="Pixel size")


  (options,args) = parser.parse_args()
  if (len(args)!=0):
      parser.error("incorrect number of arguments. Use -h for help")


  InputModel = options.InputModel
  ObsTime = options.obstime
  offset = options.offset
  emin = options.loweren
  emax = options.upperen
  bins = options.bins
  irf = options.irf
  caldb = options.caldb
  sigma = options.sigma
  rad = options.radius
  binsz = options.binsz

  #####################################################

  print '************************************'
  print '** '+ScriptName
  print '************************************'

  print '** Options:'

  #List parameters and args
  print '\n**Input parameters:'
  for key, val in parser.values.__dict__.iteritems():
    print key,':', val
    if (val==None):
        print '\nERROR!',key,'is a required option! Exit.\nType ',ScriptName+'.py -h for help\n'
        sys.exit(1)


# xml file with the input for the GRB (needed for the spectral shape)
  models = gammalib.GModels(InputModel)
  TrueModel='grb.xml'
  models.save(TrueModel)

  outfile="sensi-%dsigma_obstime-%.fs_irf-%s.txt"%(sigma,ObsTime,str(irf))

  cmd='cssens inmodel='+TrueModel +' srcname=GRB caldb='+str(caldb)+' irf='+str(irf) +'  rad='+str(rad) +' emin='+str(emin) +' emax='+str(emax) +' type=Integral sigma='+str(sigma) +' bins='+str(bins) +' duration='+str(ObsTime) +' outfile='+outfile+' binsz='+str(binsz) +' offset='+str(offset)



  print cmd
  os.system(cmd)
