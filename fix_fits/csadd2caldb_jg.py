#! /usr/bin/env python
# ==========================================================================
# Add IRFs to CALDB
#
# Copyright (C) 2021 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import os
import sys
import glob
import tarfile
import gammalib
import ctools
from shutil import copyfile


# ================= #
# csadd2caldb class #
# ================= #
class csadd2caldb(ctools.cscript):
    """
    Add IRFs to CALDB
    """

    # Constructor
    def __init__(self, *argv):
        """
        Constructor
        """
        # Make sure that parfile exists
        self._parfile()

        # Initialise application by calling the base class constructor
        self._init_cscript(self.__class__.__name__, ctools.__version__, argv)

        # Return
        return


    # Private methods
    def _parfile(self):
        """
        Check if parfile exists. If parfile does not exist then create a
        default parfile. This kluge avoids shipping the cscript with a parfile.
        """
        # Set parfile name
        parfile = self.__class__.__name__+'.par'

        # Load parfile or create
        try:
            pars = gammalib.GApplicationPars(parfile)
        except:
            # Signal that parfile was not found
            print('Parfile '+parfile+' not found. Create default parfile.')

            # Create default parfile
            pars = gammalib.GApplicationPars()
            pars.append(gammalib.GApplicationPar('indir','f','a','fits','','','Input IRF folder'))
            pars.append(gammalib.GApplicationPar('outdir','f','a','$CALDB','','','Output caldb folder'))
            pars.append_standard()
            pars.save(parfile)

        # Return
        return

    def _get_parameters(self):
        """
        Get parameters from parfile
        """
        # Query parameters
        self['indir'].filename()
        self['outdir'].filename()

        #  Write input parameters into logger
        self._log_parameters(gammalib.TERSE)

        # Return
        return
    
    def _collect_fits_gz_files(self, indir):
        """
        Collect fits.gz files

        Parameters
        ----------
        indir : string
            Input directory

        Returns
        -------
        files : list of str
            List of files
        """
        # Initialise list of tarfiles
        files = []

        # Search for tarfiles in current directory
        first_files = glob.glob('%s/*.fits.gz' % indir.url())
        files.extend(first_files)

        # Search for tarfiles in 'fits' subdirectory
        second_files = glob.glob('%s/fits/*.fits.gz' % indir.url())
        files.extend(second_files)
        
        print(files)

        # Return tarfiles
        return files
        

    def _collect_tarfiles(self, indir):
        """
        Collect tarfiles

        Parameters
        ----------
        indir : string
            Input directory

        Returns
        -------
        tarfiles : list of str
            List of tarfiles
        """
        # Initialise list of tarfiles
        tarfiles = []

        # Search for tarfiles in current directory
        files = glob.glob('%s/*.fits.tar.gz' % indir.url())
        tarfiles.extend(files)

        # Search for tarfiles in 'fits' subdirectory
        files = glob.glob('%s/fits/*.FITS.tar.gz' % indir.url())
        tarfiles.extend(files)

        # Return tarfiles
        return tarfiles

    def _extract_irfs(self, file, outdir, filetype="fitsgz"):
        """
        Collect IRFs from tarfile

        Parameters
        ----------
        file : str or list
            Tarfile (or fitsgz file) with IRFs
        outdir : str
            Output directory
        filetype: str
            tarfile or fitsgz

        Returns
        -------
        irfs : list
            List of IRF dictionaries
        """
        # Initialise list of extracted IRFs
        irfs = []

        # Open tarfile
        if filetype == "tarfile":
            tar = tarfile.open(file)
            members = tar.getmembers()
        elif filetype == "fitsgz":
            if not isinstance(file, list):
                members = [file]

        # Loop over all members in tarfile
        for member in members:
            
            if filetype == "tarfile":
                name = member.name
            elif filetype == "fitsgz":
                if not isinstance(file, list):
                    name = member

            # Get subarray
            subarrays  = ''
            nsubarrays = 0
            if 'LST' in name:
                subarrays  += '_LST'
                nsubarrays += 1
            if 'MST' in name:
                subarrays  += '_MST'
                nsubarrays += 1
            if 'SST' in name:
                subarrays  += '_SST'
                nsubarrays += 1
            if nsubarrays > 1:
                subarrays = ''

            # Extract Instrument
            if 'South-' in name:
                site = 'South'
            else:
                site = 'North'
            if '180000s' in name:
                duration = '50h'
            elif '18000s' in name:
                duration = '5h'
            elif '1800s' in name:
                duration = '0.5h'
            if 'SouthAz' in name:
                orientation = '_S'
                azimuth     = 180.0
            elif 'NorthAz' in name:
                orientation = '_N'
                azimuth     = 0.0
            else:
                orientation = ''
                azimuth     = 90.0
            
            elements   = name.split('-')
            zenith     = [e for e in elements if 'deg' in e][0].strip('deg')
            instrument = '%s_z%s%s_%s%s' % (site, zenith, orientation, duration, subarrays)

            # Build output directory
            path = '%s/%s' % (outdir, instrument)

            # Write header
            self._log_header3(gammalib.TERSE, 'Extract IRFs from "%s"' % name)
            self._log_value(gammalib.NORMAL, 'Instrument', instrument)

            # Extract file
            output_path = '%s/%s' % (self['outdir'].filename().url(), path)
            
            if filetype == "tarfile":
                tar.extract(member, path=output_path)
            else:
                if not os.path.exists(output_path):
                    print(output_path , "does not exist yet... creating")
                    os.makedirs(output_path, exist_ok=True)
                print(f"Copying {member} to {output_path}")
                copyfile(member, output_path + "/" + member.split("/")[-1])

            # Create IRF record
            record = {'path': path, 'instrument': instrument, 'file': name,
                      'site': site, 'duration': duration, 'zenith': zenith, 'azimuth': azimuth,
                      'orientation': orientation}

            # Add record
            irfs.append(record)

        if filetype == "tarfile":
            # Close tarfile
            tar.close()

        # Return
        return irfs

    def _add_irf_to_cif(self, prod, hdu, irf):
        """
        Add IRF to calibration index file HDU

        Parameters
        ----------
        prod : str
            Production
        hdu : ~gammalib.GFitsBinTable
            CIF HDU table
        irf : dict
            IRF dictionary
        """
        # Set list of IRF component names
        #names = ['EA', 'PSF', 'EDISP', 'BKG']
        components = [{'name': 'EA',    'CNAM': 'EFF_AREA', 'DESC': 'CTA effective area'},
                      {'name': 'PSF',   'CNAM': 'RPSF',     'DESC': 'CTA point spread function'},
                      {'name': 'EDISP', 'CNAM': 'EDISP',    'DESC': 'CTA energy dispersion'},
                      {'name': 'BKG',   'CNAM': 'BKG',      'DESC': 'CTA background'}]

        # Initialise CIF row index
        row = hdu.nrows()

        # Append rows for all components to CIF extension
        hdu.append_rows(len(components))

        # Add information for all components
        for component in components:

            # Set calibration information
            cal_name     = 'NAME(%s)'         % irf['instrument']
            cal_version  = 'VERSION(%s)'      % prod
            cal_cut      = 'CLASS(BEST)'
            cal_analysis = 'ANALYSIS(CTA)'
            cal_zenith   = 'ZENITH(%.3f)deg'  % float(irf['zenith'])
            cal_azimuth  = 'AZIMUTH(%.3f)deg' % float(irf['azimuth'])
            cal_bounds   = [cal_name, cal_version, cal_cut, cal_analysis, \
                            cal_zenith, cal_azimuth]

            # Set generic information
            hdu['TELESCOP'][row] = 'CTA'
            hdu['INSTRUME'][row] = gammalib.toupper(prod)
            hdu['DETNAM'][row]   = 'NONE'
            hdu['FILTER'][row]   = 'NONE'
            hdu['CAL_DEV'][row]  = 'ONLINE'
            hdu['CAL_CLAS'][row] = 'BCF'
            hdu['CAL_DTYP'][row] = 'DATA'
            hdu['CAL_VSD'][row]  = '2014-01-30'
            hdu['CAL_VST'][row]  = '00:00:00'
            hdu['REF_TIME'][row] = 51544.0
            hdu['CAL_QUAL'][row] = 0            # 0=good, 1=bad, 2=dubious, ...
            hdu['CAL_DATE'][row] = '14/01/30'

            # Set component specific information
            hdu['CAL_DIR'][row]   = irf['path']
            hdu['CAL_FILE'][row]  = irf['file']
            hdu['CAL_CNAM'][row]  = component['CNAM']
            hdu['CAL_DESC'][row]  = component['DESC']
            hdu['CAL_XNO'][row]   = 1
            for i in range(9):
                if i >= len(cal_bounds):
                    hdu['CAL_CBD'][row,i] = 'NONE'
                else:
                    hdu['CAL_CBD'][row,i] = cal_bounds[i]

            # Increment row index
            row += 1

        # Return
        return

    def _create_cif_table(self):
        """
        Create Calibration Database Index File binary table

        Returns
        -------
        table : ~gammalib.GFitsBinTable
            Calibration Database Index File binary table
        """
        # Create binary table
        table = gammalib.GFitsBinTable()

        # Append columns. Reference: CAL/GEN/92-008
        table.append(gammalib.GFitsTableStringCol('TELESCOP', 0, 10))
        table.append(gammalib.GFitsTableStringCol('INSTRUME', 0, 10))
        table.append(gammalib.GFitsTableStringCol('DETNAM', 0, 20))
        table.append(gammalib.GFitsTableStringCol('FILTER', 0, 10))
        table.append(gammalib.GFitsTableStringCol('CAL_DEV', 0, 20))
        table.append(gammalib.GFitsTableStringCol('CAL_DIR', 0, 70))
        table.append(gammalib.GFitsTableStringCol('CAL_FILE', 0, 70))   # Extend beyond standard
        table.append(gammalib.GFitsTableStringCol('CAL_CLAS', 0, 3))
        table.append(gammalib.GFitsTableStringCol('CAL_DTYP', 0, 4))
        table.append(gammalib.GFitsTableStringCol('CAL_CNAM', 0, 20))
        table.append(gammalib.GFitsTableStringCol('CAL_CBD', 0, 70, 9))
        table.append(gammalib.GFitsTableShortCol('CAL_XNO', 0))
        table.append(gammalib.GFitsTableStringCol('CAL_VSD', 0, 10))
        table.append(gammalib.GFitsTableStringCol('CAL_VST', 0, 8))
        table.append(gammalib.GFitsTableDoubleCol('REF_TIME', 0))
        table.append(gammalib.GFitsTableShortCol('CAL_QUAL', 0))
        table.append(gammalib.GFitsTableStringCol('CAL_DATE', 0, 8))
        table.append(gammalib.GFitsTableStringCol('CAL_DESC', 0, 70))

        # Set keywords. Reference: CAL/GEN/92-008
        table.extname('CIF')
        table.card('CIFVERSN', '1992a', 'Version of CIF format')

        # Return table
        return table

    def _add_irfs_to_cif(self, prod, irfs):
        """
        Add IRFs to calibration index file

        Parameters
        ----------
        prod : str
            Production
        irfs : list
            List of IRFs
        """
        # Build calibration database index name
        fname = '%s/data/cta/%s/caldb.indx' % (self['outdir'].filename().url(), prod)
        
        # Open calibration database index file
        cif = gammalib.GFits(fname, True)

        # If file has no CIF extension than create it now
        if not cif.contains('CIF'):
            cif.append(self._create_cif_table())

        # Get calibration database index HDU
        cif_hdu = cif.table('CIF')

        # Loop over all IRFs and add information to CIF
        for irf in irfs:
            self._add_irf_to_cif(prod, cif_hdu, irf)

        # Save and close CIF
        cif.save(True)
        cif.close()

        # Return
        return

    def _add_irfs(self, filetype="fitsgz"):
        """
        Add IRFs to calibration database
        
        filetype = "fitsgz" or "tarfile"
        """
        # Collect tarfiles
        
        if filetype == "tarfile":
            files = self._collect_tarfiles(self['indir'].filename())
        elif filetype == "fitsgz":
            files = self._collect_fits_gz_files(self['indir'].filename())

        # Add IRFs in tarfiles
        for file in files:

            # Extract Production
            fname      = os.path.basename(file)
            elements   = fname.split('-')
            prod       = [e for e in elements if 'prod' in e][0]
            version    = [e for e in elements if 'v' in e][0]
            prod       = '%s-%s' % (prod, version)

            # Build output directory
            outdir = 'data/cta/%s/bcf' % (prod)

            # Write header and log parameters
            self._log_header2(gammalib.TERSE, 'Extract IRFs')
            self._log_value(gammalib.NORMAL, filetype.capitalize(), file)
            self._log_value(gammalib.NORMAL, 'Production', prod)
            self._log_value(gammalib.NORMAL, 'Target directory', outdir)

            # Extract IRFs from tarfile
            irfs = self._extract_irfs(file, outdir, filetype="fitsgz")

            # Add IRFs to calibration index database
            self._add_irfs_to_cif(prod, irfs)

        # Return
        return


    # Public methods
    def run(self):
        """
        Run the script
        """
        # Switch screen logging on in debug mode
        if self._logDebug():
            self._log.cout(True)

        # Get parameters
        self._get_parameters()

        # Write header
        self._log_header1(gammalib.TERSE, 'Add IRFs to CALDB')

        # Add IRFs to CALDB
        self._add_irfs()

        # Return
        return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Create instance of application
    app = csadd2caldb(sys.argv)

    # Execute application
    app.execute()
