import numpy as np
import pandas as pd
from astropy.io import ascii
from astropy.table import Table
from scipy.interpolate import interp2d


class EBL:
    def __init__(self, model):

        self.model_file = model

        self._import_model()

        self.interpolation = interp2d(
            self._model.index.to_numpy() * 1e6,
            np.asfarray(self._model.columns.to_numpy()),
            z=self.model.to_numpy().T,
            kind="linear",
        )

    def _import_model(self):

        df = pd.read_csv(
            self.model_file, comment="#", delim_whitespace=True, index_col=0
        )

        # set index name of DataFrame
        df.index.name = "energy"

        self._model = df

    @property
    def model(self):
        return self._model

    def get(self, energy, z):
        """Get EBL value at given energy and redshift"""

        return self.interpolation(energy, z)[0]

    def abs(self, energy, z):
        """Get absorption at given energy and redshift"""

        return np.exp(-1 * self.get(energy, z))


class Spectrum:
    """Base class for a spectrum"""

    def __init__(self, ebl_model=None):

        if ebl_model is not None:
            self.ebl = EBL(ebl_model)
        else:
            self.ebl = None


class PowerLaw(Spectrum):
    """Power Law Spectrum"""

    def __init__(self, prefactor, pivot, index, z=None, ebl_model=None):

        Spectrum.__init__(self, ebl_model)

        self.z = z
        self.prefactor = prefactor
        self.pivot = pivot
        self.index = index

    def S(self, energy):
        """Calculate the spectrum at given energy and optional redshift"""

        if self.z and self.ebl:
            absorption = self.ebl.abs(energy, self.z)
        else:
            absorption = 1

        return self.prefactor * ((energy / self.pivot) ** self.index) * absorption

    def export_file_function(self, output_file, n_points=100):

        energies = (
            np.logspace(
                np.log10(self.ebl._model.index.min()),
                np.log10(self.ebl._model.index.max()),
                n_points,
            )
            * 1e6
        )

        spectrum = self.S(energies)  # Model in Tev, -> MeV

        t = Table()
        t["energy [MeV]"] = energies
        t["sepctrum [ph/cm2/s/MeV]"] = spectrum

        ascii.write(t, output_file, overwrite=True)

        xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<source_library title="source library">
  <source name="GRB" type="PointSource">
    <spectrum type="FileFunction" file="{output_file}">
        <parameter name="Normalization" scale="1.0" value="1.0" min="0.0" max="1000.0" free="1"/>
    </spectrum>
    <spatialModel type="PointSource">
      <parameter name="RA" value="83.6331" scale="1" min="-360" max="360" free="0" />
      <parameter name="DEC" value="22.0145" scale="1" min="-90" max="90" free="0" />
    </spatialModel>
  </source>
  <source name="CTABackgroundModel" type="CTAIrfBackground" instrument="CTA">
    <spectrum type="PowerLaw">
      <parameter name="Prefactor" value="1" error="0" scale="1" min="0.001" max="1000" free="1" />
      <parameter name="Index" value="0" error="0" scale="1" min="-5" max="5" free="1" />
      <parameter name="PivotEnergy" value="1" scale="1000000" min="0.01" max="1000000" free="0" />
    </spectrum>
  </source>
</source_library>
"""

        # write xml to file
        with open(output_file.replace(".dat", ".xml"), "w") as f:
            f.write(xml)
