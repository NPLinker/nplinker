#!/usr/bin/env python
import sys
from setuptools import setup


if sys.version_info[:2] < (3, 9):
    raise RuntimeError("Python version >= 3.9 required.")

setup(scripts=["bin/install-nplinker-deps"])
