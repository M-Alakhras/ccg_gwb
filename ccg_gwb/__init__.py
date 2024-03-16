"""Top-level package for ccg_gwb."""

import os
import sys

from appdirs import user_cache_dir

CCG_ENV = os.environ.copy()
CCG_PATH = sys.exec_prefix + "/bin"
CCG_ENV["PATH"] = CCG_PATH + ":" + CCG_ENV["PATH"]
CCG_CACHEDIR = user_cache_dir(appname="ccg_gwb", appauthor=False)

__author__ = """Mohammad Alakhras"""
__email__ = "mohammadalakhras1989@gmail.com"
__version__ = "0.1.0"
