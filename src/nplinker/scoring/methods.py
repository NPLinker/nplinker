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

from nplinker.logconfig import LogConfig


logger = LogConfig.getLogger(__file__)


# CG: TODO refactor this class to abstract base class
class ScoringMethod():
    """Base class of scoring methods
    """

    NAME = 'ScoringMethod'

    def __init__(self, npl):
        self.npl = npl
        self.name = self.__class__.NAME

    @staticmethod
    def setup(npl):
        """Perform any one-off initialisation required (will only be called once)"""
        pass

    def get_links(self, objects, link_collection):
        """Given a set of objects, return link information"""
        return link_collection

    def format_data(self, data):
        """Given whatever output data the method produces, return a readable string version"""
        return ''

    def sort(self, objects, reverse=True):
        """Given a list of objects, return them sorted by link score"""
        return objects
