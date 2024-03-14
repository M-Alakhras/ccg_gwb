# timing_model_simulation.py
"""
Simulate timing models.
"""

from ccg_gwb.simulation.ATNF_utilities import query_ATNF, ATNF_ephemeris, parse_ephem
from ccg_gwb.simulation.timing_model_parameters import write_par_file
import os

class TimingModel_Simulator(object):

    def __init__(self, ATNF=False, ATNF_Condition='', outdir=None):

        if outdir is None:
            outdir = os.getcwd()
        
        self._ATNF = ATNF
        self._ATNF_Condition = ATNF_Condition
        self._outdir = os.path.abspath(outdir)

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

    @property
    def ATNF(self):
        return self._ATNF

    @ATNF.setter
    def ATNF(self, value):
        self._ATNF = value

    @property
    def ATNF_Condition(self):
        return self._ATNF_Condition

    @ATNF_Condition.setter
    def ATNF_Condition(self, value):
        self._ATNF_Condition = value

    @property
    def outdir(self):
        return self._outdir

    @outdir.setter
    def outdir(self, value):
        self._outdirt = os.path.abspath(value)


    def start(self):
        if self.ATNF:
            query = query_ATNF(condition=self.ATNF_Condition)
            ephems = ATNF_ephemeris(query)
            all_params = [parse_ephem(ephem, quiet=True) for ephem in ephems]
            valid_params = [params for params in all_params if len(params) > 0]
            for params in valid_params:
                PSR = [param.value for param in params if param.name == 'PSR'][0]
                write_par_file(params, outfile=self.outdir+"/"+PSR+".par")
        else:
            # TODO: simulate parameters.
            pass
        