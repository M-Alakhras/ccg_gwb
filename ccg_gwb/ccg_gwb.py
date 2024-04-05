# ccg_gwb.py
"""Driver module for other ccg_gwb functionalities."""

import datetime
import glob
import os
import pickle

from loguru import logger as log

from ccg_gwb import CCG_CACHEDIR


class CCG_Session:
    """Driver class of ccg_gwb package."""

    @staticmethod
    def print_all_sessions():
        available_sessions = glob.glob(CCG_CACHEDIR + "/*.session")
        if len(available_sessions) > 0:
            print("Previeous sessions: ")
            for i, session in enumerate(available_sessions):
                session_name = os.path.basename(session).split(".")[0]
                print(f"   {i+1}- {session_name}")
        else:
            print("No previeous sessions.")

    @classmethod
    def delete_all_sessions(cls):
        """WARNING:: better not to use this method!"""
        available_sessions = glob.glob(CCG_CACHEDIR + "/*.session")
        if len(available_sessions) > 0:
            for i, session in enumerate(available_sessions):
                session_name = os.path.basename(session).split(".")[0]
                cls.delete(session_name)
        else:
            log.warning("No previeous sessions.")

    def __init__(self):
        self.time = datetime.datetime.now()
        log.info("A new session has been started.")

    @staticmethod
    def load(session_name):
        try:
            with open(CCG_CACHEDIR + "/" + session_name + ".session", "rb") as f:
                session = pickle.load(f)
            log.info(f"session {session_name} loaded successfully!")
            return session
        except FileNotFoundError:
            log.error(f"session {session_name} not found!")
            return None

    def save(self):
        session_name = self.time.isoformat(sep=" ").split(".")[0]
        session_file = CCG_CACHEDIR + "/" + session_name + ".session"
        with open(session_file, "wb") as f:
            pickle.dump(self, f)
        self.last_save = datetime.datetime.now()
        time_stamp = self.last_save.isoformat(sep=" ").split(".")[0]
        log.info(f"session {session_name} saved at {time_stamp}")

    @staticmethod
    def delete(session_name):
        try:
            session_file = CCG_CACHEDIR + "/" + session_name + ".session"
            os.remove(session_file)
            log.info(f"session {session_name} has been deleted!")
        except FileNotFoundError:
            log.error(f"session {session_name} not found!")

    def load_data(self, TOAs, models, TOAs_dir=None, models_dir=None):
        pass
        # if TOAs is None:
        #     if TOAs_dir is None:
        #         pass
        # self._toas = TOAs
        # self._models = models
        # TODO:
        # - create pulsars

    def explore_data(self):
        pass

    def simulate_data(self, parameter_file=None):
        pass

    def preprocessing_data(self):
        pass

    def get_analyzer(self):
        pass
