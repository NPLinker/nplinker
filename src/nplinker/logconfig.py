# Copyright 2021 The NPLinker Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import sys
from typing_extensions import Self


class LogConfig():

    active_loggers: dict[str, logging.Logger] = {}
    logfmt = '%(asctime)s [%(levelname)s] %(filename)s:%(lineno)d, %(message)s'
    default_loglevel = logging.INFO
    # default destination for new Loggers
    default_logdest = logging.StreamHandler(sys.stdout)
    # additional destinations to be added to new Loggers
    additional_logdests: list[logging.Handler] = []

    @staticmethod
    def getLogger(obj: str, level=default_loglevel, dest=default_logdest) -> logging.Logger:
        """Return a logging.Logger associated with the object <obj>.

        The Logger's level and dest values will be set to the corresponding
        parameters passed to this method.
        """
        if obj in LogConfig.active_loggers:
            return LogConfig.active_loggers[obj]

        logger = logging.Logger(obj)
        logger.setLevel(level)
        dest.setFormatter(
            logging.Formatter(LogConfig.logfmt, datefmt='%H:%M:%S'))
        logger.addHandler(dest)
        LogConfig.active_loggers[obj] = logger
        return logger

    @staticmethod
    def setLogLevel(level):
        """Apply a new log level value to all loggers created by getLogger"""
        LogConfig.default_loglevel = level
        for logger in LogConfig.active_loggers.values():
            logger.setLevel(level)

    @staticmethod
    def setLogLevelStr(level):
        """Apply a new log level value to all loggers created by getLogger

        Identical to setLogLevel but parameter is a string instead of a
        constant from the logging module (e.g. "INFO", "DEBUG")
        """
        if not hasattr(logging, level):
            raise Exception(f'Unknown/invalid loglevel "{level}"')

        LogConfig.setLogLevel(getattr(logging, level))

    @staticmethod
    def setLogDestination(dest):
        LogConfig.default_logdest = dest
        LogConfig.additional_logdests = []
        dest.setFormatter(logging.Formatter(LogConfig.logfmt))
        for logger in LogConfig.active_loggers.values():
            logger.handlers = []
            logger.addHandler(dest)

    @staticmethod
    def addLogDestination(dest):
        LogConfig.additional_logdests.append(dest)
        dest.setFormatter(logging.Formatter(LogConfig.logfmt))
        for logger in LogConfig.active_loggers.values():
            logger.addHandler(dest)
