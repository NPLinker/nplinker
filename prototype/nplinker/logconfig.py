import sys
import logging

class LogConfig(object):

    active_loggers = {}
    logfmt = '%(asctime)s [%(levelname)s] %(filename)s:%(lineno)d, %(message)s'
    default_loglevel = logging.INFO
    default_logdest = logging.StreamHandler(sys.stdout)

    @staticmethod
    def getLogger(obj, level=default_loglevel, dest=default_logdest):
        if obj in LogConfig.active_loggers:
            return LogConfig.active_loggers[obj]

        logger = logging.Logger(obj)
        logger.setLevel(level)
        dest.setFormatter(logging.Formatter(LogConfig.logfmt, datefmt='%H:%M:%S'))
        logger.addHandler(dest)
        LogConfig.active_loggers[obj] = logger
        return logger

    @staticmethod
    def setLogLevel(level):
        LogConfig.default_loglevel = level
        for logger in LogConfig.active_loggers.values():
            logger.setLevel(level)

    @staticmethod
    def setLogLevelStr(level):
        if not hasattr(logging, level):
            raise Exception('Unknown/invalid loglevel "{}"'.format(level))

        LogConfig.setLogLevel(getattr(logging, level))

    @staticmethod
    def setLogDestination(dest):
        LogConfig.default_logdest = dest
        if isinstance(dest, str):
            dest = logging.FileHandler(dest)
        dest.setFormatter(logging.Formatter(LogConfig.logfmt))
        for logger in LogConfig.active_loggers.values():
            logger.handlers = []
            logger.addHandler(dest)
