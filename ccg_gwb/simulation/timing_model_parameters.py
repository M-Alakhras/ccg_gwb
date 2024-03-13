# timing_model_parameters
"""
Timing model parameters calsses.
"""

class Parameter(object):
    
    def __init__(self, name, value=None, error=None, freeze=True):
        self._name = name
        self._value = value
        self._error = error
        self._freeze = freeze
        self._valid = True

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self):
        print("Warning:: Parameter name can't be changed.")

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        self._value = value

    @property
    def error(self):
        return self._error

    @error.setter
    def error(self, value):
        self._error = value

    @property
    def freeze(self):
        return self._freeze

    @freeze.setter
    def freeze(self, value):
        self._freeze = value

    @property
    def valid(self):
        return self._valid

    @valid.setter
    def valid(self, value):
        print("Warning:: Validity of parameter is read-only.")

    def auto_detect(self):
        if self.name in ['PSR', 'EPHEM', 'CLK', 'UNITS']:
            return MiscellaneousParameter(self.name, value=self.value)
        if self.name in ['POSEPOCH', 'PX', 'RAJ', 'DECJ', 'PMRA', 'PMDEC', 'ELONG', 
                         'ELAT', 'PMELONG', 'PMELAT', 'GL', 'GB', 'PML', 'PMB']:
            return AstrometryParameter(self.name, value=self.value, error=self.error, freeze=self.freeze)
        DM_deriv = ['DM'+str(i) for i in range(1,10)]
        if self.name in ['DM', 'DMEPOCH']+DM_deriv:
            return DispersionParameter(self.name, value=self.value, error=self.error, freeze=self.freeze)
        F_deriv = ['F'+str(i) for i in range(10)]
        if self.name in ['PEPOCH']+F_deriv:
            return SpindownParameter(self.name, value=self.value, error=self.error, freeze=self.freeze)
        if self.name in ['BINARY']:
            return BinaryModel(self.value)
        FB_deriv = ['FB'+str(i) for i in range(10)]
        if self.name in ['PB', 'PBDOT', 'A1', 'A1DOT', 'ECC', 'EDOT', 
                         'T0', 'OM', 'OMDOT', 'M2', 'SINI', 'A0', 'B0', 
                         'GAMMA', 'DR', 'DTHETA', 'H3', 'H4', 'STIGMA', 'KIN', 
                         'KOM', 'K96', 'SHAPMAX', 'TASC', 'EPS1', 'EPS2', 
                         'EPS1DOT', 'EPS2DOT', 'NHARMS', 'LNEDOT']+FB_deriv:
            return BinaryParameter(self.name, binary_model='Generic', value=self.value, error=self.error, freeze=self.freeze)
        return self



class MiscellaneousParameter(Parameter):

    def __init__(self, name, value=None):
        super().__init__(name, value=value)
        self.check_validity()

    def check_validity(self):
        if self.name not in ['PSR', 'EPHEM', 'CLK', 'UNITS']:
            self._valid = False



class AstrometryParameter(Parameter):

    def __init__(self, name, value=None, error=None, freeze=True):
        super().__init__(name, value=value, error=error, freeze=freeze)
        self.check_validity()

    @property
    def astrometry(self):
        return self._astrometry

    @astrometry.setter
    def astrometry(self, value):
        print("Warning:: Astrometry type can't be changed.")
    
    def check_validity(self):
        self._astrometry = None
        if self.name in ['POSEPOCH', 'PX']:
            self._astrometry = 'General'
        if self.name in ['RAJ', 'DECJ', 'PMRA', 'PMDEC']:
            self._astrometry = 'Equatorial'
        if self.name in ['ELONG', 'ELAT', 'PMELONG', 'PMELAT']:
            self._astrometry = 'Ecliptic'
        if self.name in ['GL', 'GB', 'PML', 'PMB']:
            self._astrometry = 'Galactic'
        if self.astrometry is None:
            self._valid = False



class DispersionParameter(Parameter):

    def __init__(self, name, value=None, error=None, freeze=True):
        super().__init__(name, value=value, error=error, freeze=freeze)
        self.check_validity()

    def check_validity(self):
        DM_deriv = ['DM'+str(i) for i in range(1,10)]
        if self.name not in ['DM', 'DMEPOCH']+DM_deriv:
            self._valid = False



class SpindownParameter(Parameter):

    def __init__(self, name, value=None, error=None, freeze=True):
        super().__init__(name, value=value, error=error, freeze=freeze)
        self.check_validity()

    def check_validity(self):
        F_deriv = ['F'+str(i) for i in range(10)]
        if self.name not in ['PEPOCH']+F_deriv:
            self._valid = False

    

class BinaryModel(Parameter):

    def __init__(self, name):
        super().__init__('BINARY', value=name, error=None, freeze=True)
        self.check_validity()

    def check_validity(self):
        if self.value not in ['Generic', 'BT', 'DD', 'DDH', 'DDK', 'DDS', 
                              'ELL1', 'ELL1H', 'ELL1K']:
            self._valid = False



class BinaryParameter(Parameter):

    def __init__(self, name, binary_model=None, value=None, error=None, freeze=True):
        super().__init__(name, value=value, error=error, freeze=freeze)
        self._binary_model = BinaryModel(binary_model)
        self.check_validity()

    @property
    def binary_model(self):
        return self._binary_model

    @binary_model.setter
    def binary_model(self, value):
        self._binary_model = value
        self.check_validity()

    def check_validity(self):
        FB_deriv = ['FB'+str(i) for i in range(10)]
        if self.binary_model == 'BT':
            if self.name not in ['PB', 'PBDOT', 'A1', 'A1DOT', 'ECC', 'EDOT', 
                                 'T0', 'OM', 'OMDOT', 'GAMMA']+FB_deriv:
                self._valid = False
        if self.binary_model == 'DD':
            if self.name not in ['PB', 'PBDOT', 'A1', 'A1DOT', 'ECC', 'EDOT', 
                                 'T0', 'OM', 'OMDOT', 'M2', 'SINI', 'A0', 'B0', 
                                 'GAMMA', 'DR', 'DTHETA']+FB_deriv:
                self._valid = False
        if self.binary_model == 'DDH':
            if self.name not in ['PB', 'PBDOT', 'A1', 'A1DOT', 'ECC', 'EDOT', 
                                 'T0', 'OM', 'OMDOT', 'A0', 'B0', 'GAMMA', 'DR', 
                                 'DTHETA', 'H3', 'STIGMA']+FB_deriv:
                self._valid = False
        if self.binary_model == 'DDK':
            if self.name not in ['PB', 'PBDOT', 'A1', 'A1DOT', 'ECC', 'EDOT', 
                                 'T0', 'OM', 'OMDOT', 'M2', 'GAMMA', 'A0', 'B0', 
                                 'DR', 'DTHETA', 'KIN', 'KOM', 'K96']+FB_deriv:
                self._valid = False
        if self.binary_model == 'DDS':
            if self.name not in ['PB', 'PBDOT', 'A1', 'A1DOT', 'ECC', 'EDOT', 
                                 'T0', 'OM', 'OMDOT', 'M2', 'GAMMA', 'A0', 'B0', 
                                 'DR', 'DTHETA', 'SHAPMAX']+FB_deriv:
                self._valid = False
        if self.binary_model == 'ELL1':
            if self.name not in ['PB', 'PBDOT', 'A1', 'A1DOT', 'M2', 'SINI', 
                                 'TASC', 'EPS1', 'EPS2', 'EPS1DOT', 'EPS2DOT']+FB_deriv:
                self._valid = False
        if self.binary_model == 'ELL1H':
            if self.name not in ['PB', 'PBDOT', 'A1', 'A1DOT', 'TASC', 'EPS1', 
                                 'EPS2', 'EPS1DOT', 'EPS2DOT', 'H3', 'H4', 
                                 'STIGMA', 'NHARMS']+FB_deriv:
                self._valid = False
        if self.binary_model == 'ELL1K':
            if self.name not in ['PB', 'PBDOT', 'A1', 'A1DOT', 'M2', 'SINI', 
                                 'TASC', 'EPS1', 'EPS2', 'OMDOT', 'LNEDOT']+FB_deriv:
                self._valid = False
            



def validate_parameters(params):
    default_values = {'EPHEM': 'DE421', 'CLK': 'TT(BIPM2021)', 'UNITS': 'TDB', 
                      'POSEPOCH': '55000'}

    # Check miscellaneous
    miscellaneous_params = [param for param in params if type(param) is MiscellaneousParameter]
    miscellaneous_params_names = [param.name for param in params if type(param) is MiscellaneousParameter]
    if 'PSR' not in miscellaneous_params_names:
        return [], params
    for param_name in ['EPHEM', 'CLK', 'UNITS']:
        if param_name not in miscellaneous_params_names:
            param = Parameter(param_name, value=default_values[param_name])
            miscellaneous_params += [param.auto_detect()]
    
    # Check astrometry
    astrometry_params = [param for param in params if type(param) is AstrometryParameter]
    astrometry_params_names = [param.name for param in params if type(param) is AstrometryParameter]
    if 'POSEPOCH' not in astrometry_params_names:
        param = Parameter('POSEPOCH', value=default_values['POSEPOCH'])
        astrometry_params += [param.auto_detect()]
    if 'RAJ' in astrometry_params_names and 'DECJ' in astrometry_params_names:
        astrometry_params = [param for param in astrometry_params if param._astrometry in ['General', 'Equatorial']]
    elif 'ELONG' in astrometry_params_names and 'ELAT' in astrometry_params_names:
        astrometry_params = [param for param in astrometry_params if param._astrometry in ['General', 'Ecliptic']]
    elif 'GL' in astrometry_params_names and 'GB' in astrometry_params_names:
        astrometry_params = [param for param in astrometry_params if param._astrometry in ['General', 'Galactic']]
    else:
        return [], params

    # Check dispersion
    dispersion_params = [param for param in params if type(param) is DispersionParameter]
    dispersion_params_names = [param.name for param in params if type(param) is DispersionParameter]

    # Check spindown
    spindown_params = [param for param in params if type(param) is SpindownParameter]
    spindown_params_names = [param.name for param in params if type(param) is SpindownParameter]
    if 'PEPOCH' not in spindown_params_names:
        posepoch = [param.value for param in astrometry_params if param.name == 'POSEPOCH']
        param = Parameter('PEPOCH', value=posepoch[0])
        spindown_params += [param.auto_detect()]

    # Check binary model
    binary_model = [param for param in params if type(param) is BinaryModel]
    binary_params = [param for param in params if type(param) is BinaryParameter]
    if len(binary_model) > 0:
        for param in binary_params:
            param.binary_model = binary_model[0].value
        if not binary_model[0].valid:
            return [], params

    pint_params = miscellaneous_params + astrometry_params + dispersion_params + \
                  spindown_params + binary_model + binary_params
    extra_params = [param for param in params if param not in pint_params]
    return pint_params, extra_params