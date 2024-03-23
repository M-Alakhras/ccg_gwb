# timing_model_parameters
"""
Timing model parameters calsses.
"""
import os


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
        if self.name in ["PSR", "EPHEM", "CLK", "UNITS"]:
            return MiscellaneousParameter(self.name, value=self.value)
        if self.name in [
            "POSEPOCH",
            "PX",
            "RAJ",
            "DECJ",
            "PMRA",
            "PMDEC",
            "ELONG",
            "ELAT",
            "PMELONG",
            "PMELAT",
            "GL",
            "GB",
            "PML",
            "PMB",
        ]:
            return AstrometryParameter(self.name, value=self.value, error=self.error, freeze=self.freeze)
        DM_deriv = ["DM" + str(i) for i in range(1, 10)]
        if self.name in ["DM", "DMEPOCH"] + DM_deriv:
            return DispersionParameter(self.name, value=self.value, error=self.error, freeze=self.freeze)
        F_deriv = ["F" + str(i) for i in range(10)]
        if self.name in ["PEPOCH"] + F_deriv:
            return SpindownParameter(self.name, value=self.value, error=self.error, freeze=self.freeze)
        if self.name in ["BINARY"]:
            return BinaryModel(self.value)
        FB_deriv = ["FB" + str(i) for i in range(10)]
        if (
            self.name
            in [
                "PBDOT",
                "A1",
                "A1DOT",
                "ECC",
                "EDOT",
                "T0",
                "OM",
                "OMDOT",
                "M2",
                "SINI",
                "A0",
                "B0",
                "GAMMA",
                "DR",
                "DTHETA",
                "H3",
                "H4",
                "STIGMA",
                "KIN",
                "KOM",
                "K96",
                "SHAPMAX",
                "TASC",
                "EPS1",
                "EPS2",
                "EPS1DOT",
                "EPS2DOT",
                "NHARMS",
                "LNEDOT",
            ]
            + FB_deriv
        ):
            return BinaryParameter(
                self.name, binary_model="Generic", value=self.value, error=self.error, freeze=self.freeze
            )
        return self


class MiscellaneousParameter(Parameter):

    def __init__(self, name, value=None):
        super().__init__(name, value=value)
        self.check_validity()

    def check_validity(self):
        if self.name not in ["PSR", "EPHEM", "CLK", "UNITS"]:
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
        if self.name in ["POSEPOCH", "PX"]:
            self._astrometry = "General"
        if self.name in ["RAJ", "DECJ", "PMRA", "PMDEC"]:
            self._astrometry = "Equatorial"
        if self.name in ["ELONG", "ELAT", "PMELONG", "PMELAT"]:
            self._astrometry = "Ecliptic"
        if self.name in ["GL", "GB", "PML", "PMB"]:
            self._astrometry = "Galactic"
        if self.astrometry is None:
            self._valid = False


class DispersionParameter(Parameter):

    def __init__(self, name, value=None, error=None, freeze=True):
        super().__init__(name, value=value, error=error, freeze=freeze)
        self.check_validity()

    def check_validity(self):
        DM_deriv = ["DM" + str(i) for i in range(1, 10)]
        if self.name not in ["DM", "DMEPOCH"] + DM_deriv:
            self._valid = False


class SpindownParameter(Parameter):

    def __init__(self, name, value=None, error=None, freeze=True):
        super().__init__(name, value=value, error=error, freeze=freeze)
        self.check_validity()

    def check_validity(self):
        F_deriv = ["F" + str(i) for i in range(10)]
        if self.name not in ["PEPOCH"] + F_deriv:
            self._valid = False


class BinaryModel(Parameter):

    def __init__(self, name):
        super().__init__("BINARY", value=name, error=None, freeze=True)
        self.check_validity()

    def check_validity(self):
        if self.value not in ["Generic", "BT", "DD", "DDK", "DDS", "ELL1", "ELL1H", "ELL1K"]:
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
        FB_deriv = ["FB" + str(i) for i in range(10)]
        if self.binary_model == "BT":
            if self.name not in ["PBDOT", "A1", "A1DOT", "ECC", "EDOT", "T0", "OM", "OMDOT", "GAMMA"] + FB_deriv:
                self._valid = False
        if self.binary_model == "DD":
            if (
                self.name
                not in [
                    "PBDOT",
                    "A1",
                    "A1DOT",
                    "ECC",
                    "EDOT",
                    "T0",
                    "OM",
                    "OMDOT",
                    "M2",
                    "SINI",
                    "A0",
                    "B0",
                    "GAMMA",
                    "DR",
                    "DTHETA",
                ]
                + FB_deriv
            ):
                self._valid = False
        # if self.binary_model == 'DDH':
        #     if self.name not in ['PBDOT', 'A1', 'A1DOT', 'ECC', 'EDOT',
        #                          'T0', 'OM', 'OMDOT', 'A0', 'B0', 'GAMMA', 'DR',
        #                          'DTHETA', 'H3', 'STIGMA']+FB_deriv:
        #         self._valid = False
        if self.binary_model == "DDK":
            if (
                self.name
                not in [
                    "PBDOT",
                    "A1",
                    "A1DOT",
                    "ECC",
                    "EDOT",
                    "T0",
                    "OM",
                    "OMDOT",
                    "M2",
                    "GAMMA",
                    "A0",
                    "B0",
                    "DR",
                    "DTHETA",
                    "KIN",
                    "KOM",
                    "K96",
                ]
                + FB_deriv
            ):
                self._valid = False
        if self.binary_model == "DDS":
            if (
                self.name
                not in [
                    "PBDOT",
                    "A1",
                    "A1DOT",
                    "ECC",
                    "EDOT",
                    "T0",
                    "OM",
                    "OMDOT",
                    "M2",
                    "GAMMA",
                    "A0",
                    "B0",
                    "DR",
                    "DTHETA",
                    "SHAPMAX",
                ]
                + FB_deriv
            ):
                self._valid = False
        if self.binary_model == "ELL1":
            if (
                self.name
                not in ["PBDOT", "A1", "A1DOT", "M2", "SINI", "TASC", "EPS1", "EPS2", "EPS1DOT", "EPS2DOT"] + FB_deriv
            ):
                self._valid = False
        if self.binary_model == "ELL1H":
            if (
                self.name
                not in [
                    "PBDOT",
                    "A1",
                    "A1DOT",
                    "TASC",
                    "EPS1",
                    "EPS2",
                    "EPS1DOT",
                    "EPS2DOT",
                    "H3",
                    "H4",
                    "STIGMA",
                    "NHARMS",
                ]
                + FB_deriv
            ):
                self._valid = False
        if self.binary_model == "ELL1K":
            if (
                self.name
                not in ["PBDOT", "A1", "A1DOT", "M2", "SINI", "TASC", "EPS1", "EPS2", "OMDOT", "LNEDOT"] + FB_deriv
            ):
                self._valid = False


def validate_parameters(params, quiet=False):
    default_values = {"EPHEM": "DE421", "CLK": "TT(BIPM2021)", "UNITS": "TDB", "PEPOCH": "55000"}

    # Check miscellaneous
    if not quiet:
        print("Checking timing model miscellaneous parameters...")
    miscellaneous_params = [param for param in params if type(param) is MiscellaneousParameter]
    miscellaneous_params_names = [param.name for param in params if type(param) is MiscellaneousParameter]
    if "PSR" not in miscellaneous_params_names:
        if not quiet:
            print("Warning:: Pulsar name must be given.")
        return [], params
    for param_name in ["EPHEM", "CLK", "UNITS"]:
        if param_name not in miscellaneous_params_names:
            param = Parameter(param_name, value=default_values[param_name])
            miscellaneous_params += [param.auto_detect()]
    CLK_param = [param for param in miscellaneous_params if param.name == "CLK"][0]
    CLK_version = CLK_param.value.split("(")[-1][:-1]
    if "BIPM" in CLK_version and CLK_version == "BIPM":
        CLK_param.value = default_values["CLK"]

    # Check spindown
    if not quiet:
        print("Checking timing model spindown parameters...")
    spindown_params = [param for param in params if type(param) is SpindownParameter]
    spindown_params_names = [param.name for param in params if type(param) is SpindownParameter]
    if "F0" not in spindown_params_names:
        if not quiet:
            print("Warning:: Spindown parameters are missing.")
        return [], params
    if "PEPOCH" not in spindown_params_names:
        param = Parameter("PEPOCH", value=default_values["PEPOCH"])
        spindown_params += [param.auto_detect()]

    # Check astrometry
    if not quiet:
        print("Checking timing model astrometry parameters...")
    astrometry_params = [param for param in params if type(param) is AstrometryParameter]
    astrometry_params_names = [param.name for param in params if type(param) is AstrometryParameter]
    if "RAJ" in astrometry_params_names and "DECJ" in astrometry_params_names:
        astrometry_params = [param for param in astrometry_params if param._astrometry in ["General", "Equatorial"]]
    elif "ELONG" in astrometry_params_names and "ELAT" in astrometry_params_names:
        astrometry_params = [param for param in astrometry_params if param._astrometry in ["General", "Ecliptic"]]
    elif "GL" in astrometry_params_names and "GB" in astrometry_params_names:
        astrometry_params = [param for param in astrometry_params if param._astrometry in ["General", "Galactic"]]
    else:
        if not quiet:
            print("Warning:: Astrometry parameters are missing.")
        return [], params
    PM = False
    for param_name in astrometry_params_names:
        if "PM" in param_name:
            PM = True
    if "POSEPOCH" not in astrometry_params_names and PM:
        pepoch = [param.value for param in spindown_params if param.name == "PEPOCH"]
        param = Parameter("POSEPOCH", value=pepoch[0])
        astrometry_params += [param.auto_detect()]

    # Check dispersion
    if not quiet:
        print("Checking timing model dispersion parameters...")
    dispersion_params = [param for param in params if type(param) is DispersionParameter]
    dispersion_params_names = [param.name for param in params if type(param) is DispersionParameter]
    DM1 = False
    for param in dispersion_params:
        if param.name == "DM1":
            DM1 = float(param.value) != 0.0
    if "DMEPOCH" not in dispersion_params_names and DM1:
        pepoch = [param.value for param in spindown_params if param.name == "PEPOCH"]
        param = Parameter("DMEPOCH", value=pepoch[0])
        spindown_params += [param.auto_detect()]

    # Check binary model
    if not quiet:
        print("Checking timing model binary parameters...")
    binary_model = [param for param in params if type(param) is BinaryModel]
    binary_params = [param for param in params if type(param) is BinaryParameter]
    if len(binary_model) > 0:
        if not binary_model[0].valid:
            if not quiet:
                print("Warning:: Unknown binary model.")
            return [], params
        for param in binary_params:
            param.binary_model = binary_model[0].value
        binary_params = [param for param in params if type(param) is BinaryParameter and param.valid]
        binary_params_name = [param.name for param in params if type(param) is BinaryParameter and param.valid]

        missing = False
        if binary_model[0].value == "BT":
            for param_name in ["T0", "A1"]:
                if param_name not in binary_params_name:
                    missing = True
        if binary_model[0].value in ["DD", "DDK", "DDS"]:
            for param_name in ["T0", "A1"]:
                if param_name not in binary_params_name:
                    missing = True
        if binary_model[0].value in ["ELL1", "ELL1H", "ELL1K"]:
            for param_name in ["TASC"]:
                if param_name not in binary_params_name:
                    missing = True
        if missing:
            if not quiet:
                print("Warning:: Binary model parameters are missing.")
            return [], params

    pint_params = (
        miscellaneous_params + astrometry_params + dispersion_params + spindown_params + binary_model + binary_params
    )

    PSR_param = [param for param in params if param.name == "PSR"]
    extra_params = PSR_param + [param for param in params if param not in pint_params]
    return pint_params, extra_params


def write_par_file(params, outfile=None):

    PSR = [param.value for param in params if param.name == "PSR"][0]

    if outfile is None:
        outfile = os.getcwd() + "/" + PSR + ".par"
    lines = []
    for param in params:
        name = param.name
        value = param.value
        error = ""
        if param.error is not None:
            error = param.error
        if param.freeze:
            fit = "0"
        else:
            fit = "1"
        if type(param) in [MiscellaneousParameter, BinaryModel]:
            fit = ""
        lines.append("{:<10} {:>22} {:1} {:<22}".format(name, value, fit, error))
    lines = "\n".join(lines)
    with open(outfile, "w") as file:
        file.write(lines)
