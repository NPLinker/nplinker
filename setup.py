#!/usr/bin/env python
import sys
from setuptools import setup


if sys.version_info[:2] < (3, 11):
    raise RuntimeError("Python version >= 3.11 required.")

setup(scripts=["bin/install-nplinker-deps"])
