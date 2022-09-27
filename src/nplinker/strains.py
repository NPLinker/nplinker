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


from .logconfig import LogConfig


logger = LogConfig.getLogger(__file__)


class Strain():

    def __init__(self, primary_strain_id):
        self.id = primary_strain_id
        self.aliases = set()

    def has_alias(self, alt_id):
        return alt_id in self.aliases

    def add_alias(self, alt_id):
        if len(alt_id) == 0:
            logger.warning(
                f'Refusing to add zero-length alias to strain {self}')
            return

        self.aliases.add(alt_id)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return f'Strain({self.id}) [{len(self.aliases)} aliases]'
