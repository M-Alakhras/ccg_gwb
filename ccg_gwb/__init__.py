"""  *                 *                        *          *
   *             *           *           *
*           ##########################################    *      *
       *    #========================================#
            #||   ======     ======     ======     ||#
   *        #||  //         //         //          ||#            *
            #|| //         //         //           ||#    *
      *     #|| ||         ||         ||           ||#
            #|| ||         ||         ||           ||#       *
*           #|| ||         ||         ||     ====  ||#   *
      *     #|| \\\         \\\         \\\      ||   ||#
   *        #||  \\\         \\\         \\\     //   ||#       *
            #||   ======     ======     ======/    ||#
            #========================================#    *
      *     ########################################## *
            *            *             *          *            *
    *          *             *   *          *           *
CCG_GWB: Gravitational Wave Background analysis tools.
A package developed by Complex & Cosmology Group (CCG) at Shahid Beheshti University (SBU).
it uses data-driven techniques to search the pulsar timing arrays (PTAs) data for gravitational
wave background (GWB) signal.

Usage:
------
from ccg_gwb.ccg_gwb import CCG_Session

# create a new session:
my_session = CCG_Session()
# load a previous session:
my_session = CCG_Session.load(session_name)
# save a session:
my_session.save()
# delete a session:
CCG_Session.delete()
"""

import glob
import os
import sys

from appdirs import user_cache_dir
from loguru import logger as log

CCG_ENV = os.environ.copy()
CCG_PATH = sys.exec_prefix + "/bin"
CCG_ENV["PATH"] = CCG_PATH + ":" + CCG_ENV["PATH"]
CCG_CACHEDIR = user_cache_dir(appname="ccg_gwb", appauthor=False)

format = "<level>{level: <8}</level> ({name: <30}): <level>{message}</level>"
log.remove(0)
log.add(sys.stderr, level="DEBUG", format=format)

if not os.path.exists(CCG_CACHEDIR):
    os.mkdir(CCG_CACHEDIR)

__author__ = """Mohammad Alakhras"""
__email__ = "mohammadalakhras1989@gmail.com"
__version__ = "0.1.1"

available_sessions = glob.glob(CCG_CACHEDIR + "/*.session")

print(__doc__)
if len(available_sessions) > 0:
    print("Previeous sessions: ")
    for i, session in enumerate(available_sessions):
        session_name = os.path.basename(session).split(".")[0]
        print(f"   {i+1}- {session_name}")
else:
    print("No previeous sessions.")
