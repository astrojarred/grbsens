#!/usr/bin/env python
#!/home/users/user_name/miniconda2/envs/my_env/bin/python
# Script name: ReadSensi.py
# Description: this script reads the output files produced with sensi.py and creates a unique file with the CTA sensitivities as a function of the observing times
# Author: B. Patricelli (barbara.patricelli@pi.infn.it)
#######################################################

# Import
import os.path


# Final output file
outfile = open('sensitivity-5sigma_irf-North_z20_0.5.txt',"w")
outfile.write('#\n')
outfile.write('# Obs time\tcrab_flux\tphoton_flux\tenergy_flux\tsensitivity\n')
outfile.write('#s\tCrab\tph/cm2/s\terg/cm2/s\terg/cm2/s\n')

# Loop over the observing times
for i in range(1,100001):
  Inputfile='sensi_North_z20/sensi-5sigma_obstime-%ds_irf-North_z20_0.5h.txt'%(i)

  if os.path.isfile(Inputfile)==True:
    print Inputfile
    infile=open(Inputfile,"r")

    crab=[]
    pflux=[]
    eflux=[]
    sensi=[]


    with open(Inputfile,"r") as f:
        first_line=f.readline()

        last_line=f.readlines()

        for line in last_line:
            crab=(float(line.split(',')[3]))
            pflux=(float(line.split(',')[4]))
            eflux=(float(line.split(',')[5]))
            sensi=(float(line.split(',')[6]))                


            outfile.write("%.d\t%f\t%e\t%e\t%e \n"%(i,crab,pflux,eflux,sensi))
