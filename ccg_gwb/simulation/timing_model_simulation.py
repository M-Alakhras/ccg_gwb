# timing_model_simulation.py
"""
Simulate timing models.
"""

import os
import subprocess

from tqdm.auto import tqdm

from ccg_gwb import CCG_ENV
from ccg_gwb.simulation.ATNF_utilities import ATNF_ephemeris, parse_ephem, query_ATNF
from ccg_gwb.simulation.timing_model_parameters import write_par_file


class TimingModel_Simulator(object):
    def __init__(self, ATNF=False, ATNF_Condition="", outdir=None, quiet=False):

        if outdir is None:
            outdir = os.getcwd()

        self._ATNF = ATNF
        self._ATNF_Condition = ATNF_Condition
        self._outdir = os.path.abspath(outdir)
        self.quiet = quiet

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
        self._outdir = os.path.abspath(value)

    def start(self):
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        if self.ATNF:
            query = query_ATNF(condition=self.ATNF_Condition)
            ephems = ATNF_ephemeris(query)
            all_params = [parse_ephem(ephem, quiet=self.quiet) for ephem in ephems]
            valid_params = [params for params in all_params if len(params[0]) > 0]
            for params in tqdm(valid_params):
                pint_params = params[0]
                extra_params = params[1]
                PSR = [param.value for param in pint_params if param.name == "PSR"][0]
                UNITS = [param.value for param in pint_params if param.name == "UNITS"][0]
                file_name = self.outdir + "/" + PSR + ".par"
                extra_file_name = self.outdir + "/" + PSR + ".extra"
                if not os.path.exists(file_name):
                    write_par_file(pint_params, outfile=file_name)
                if not os.path.exists(extra_file_name):
                    write_par_file(extra_params, outfile=extra_file_name)
                if UNITS == "TCB":
                    if not os.path.exists(file_name.replace(".par", "_tdb.par")):
                        subprocess.run(["tcb2tdb", file_name, file_name.replace(".par", "_tdb.par")], env=CCG_ENV)
        else:
            # TODO: simulate parameters.
            pass
