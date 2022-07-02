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

import math


# code to normalise peaks for spectral matching ("rosetta stone" stuff)
def sqrt_normalise(peaks):
    temp = []
    total = 0.0
    for mz, intensity in peaks:
        temp.append((mz, math.sqrt(intensity)))
        total += intensity
    norm_facc = math.sqrt(total)
    normalised_peaks = []
    for mz, intensity in temp:
        normalised_peaks.append((mz, intensity / norm_facc))
    return normalised_peaks
