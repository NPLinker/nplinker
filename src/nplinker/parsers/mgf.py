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

# coding=utf8

import os


PROTON_MASS = 1.00727645199076


class MS1():

    def __init__(self,
                 id,
                 mz,
                 rt,
                 intensity,
                 file_name,
                 scan_number=None,
                 single_charge_precursor_mass=None):
        self.id = id
        self.mz = mz
        self.rt = rt
        self.intensity = intensity
        self.file_name = file_name
        self.scan_number = scan_number
        if single_charge_precursor_mass:
            self.single_charge_precursor_mass = single_charge_precursor_mass
        else:
            self.single_charge_precursor_mass = self.mz
        self.name = f"{self.mz}_{self.rt}"

    def __str__(self):
        return 'MS1(name={}, mz={}, id={}, rt={})'.format(
            self.name, self.mz, self.rt, self.id)


# Abstract loader class
# *load_spectra* functions are too long, refactor and split when having time
class Loader():

    def __init__(self,
                 min_ms1_intensity=0.0,
                 peaklist=None,
                 isolation_window=0.5,
                 mz_tol=5,
                 rt_tol=5.0,
                 duplicate_filter_mz_tol=0.5,
                 duplicate_filter_rt_tol=16,
                 duplicate_filter=False,
                 repeated_precursor_match=None,
                 min_ms1_rt=0.0,
                 max_ms1_rt=1e6,
                 min_ms2_intensity=0.0,
                 has_scan_id=False,
                 rt_units='seconds',
                 mz_col_name='mz',
                 rt_col_name='rt',
                 csv_id_col=None,
                 id_field=None,
                 name_field=None):

        self.min_ms1_intensity = min_ms1_intensity
        self.peaklist = peaklist
        self.isolation_window = isolation_window
        self.mz_tol = mz_tol
        self.rt_tol = rt_tol
        self.duplicate_filter = duplicate_filter
        self.duplicate_filter_mz_tol = duplicate_filter_mz_tol
        self.duplicate_filter_rt_tol = duplicate_filter_rt_tol
        self.min_ms1_rt = min_ms1_rt
        self.max_ms1_rt = max_ms1_rt
        self.min_ms2_intensity = min_ms2_intensity
        if repeated_precursor_match:
            self.repeated_precursor_match = repeated_precursor_match
        else:
            self.repeated_precursor_match = 2 * self.isolation_window

        self.mz_col_name = mz_col_name
        self.rt_col_name = rt_col_name
        self.csv_id_col = csv_id_col
        self.rt_units = rt_units
        self.csv_id_col = csv_id_col
        self.id_field = id_field

        self.name_field = name_field  # only works for msp - fix for metlin people

        if not self.mz_col_name:
            self.mz_col_name = 'mz'

    def __str__(self):
        return self.__class__.__name__

    def load_spectra(self, input_set):
        raise NotImplementedError("load spectra method must be implemented")

    def _ion_masses(self, precursormass, int_charge):
        """
        Compute the parent masses. Single charge version is used for
        loss computation.
        """
        mul = abs(int_charge)
        parent_mass = precursormass * mul
        parent_mass -= int_charge * PROTON_MASS
        single_charge_precursor_mass = precursormass * mul
        if int_charge > 0:
            single_charge_precursor_mass -= (int_charge - 1) * PROTON_MASS
        elif int_charge < 0:
            single_charge_precursor_mass += (mul - 1) * PROTON_MASS
        else:
            # charge = zero - leave them all the same
            parent_mass = precursormass
            single_charge_precursor_mass = precursormass
        return parent_mass, single_charge_precursor_mass

    def _interpret_charge(self, charge):
        """
        Method to interpret the ever variable charge field in the different
        formats. Should never fail now.
        """
        if not charge:  # if it is none
            return 1
        try:
            if not type(charge) == str:
                charge = str(charge)

            # add the meat here
            # try removing any + signs
            charge = charge.replace("+", "")

            # remove trailing minus signs
            if charge.endswith('-'):
                charge = charge[:-1]
                # move the minus to the front if it
                # isn't already there
                if not charge.startswith('-'):
                    charge = '-' + charge
            # turn into an int
            int_charge = int(charge)
            return int_charge
        except:
            int_charge = 1
        return int_charge

    def _load_peak_list(self):
        """
        Modify peaklist function. Try to detect "featureid", store it in ms1_peaks
        used for mgf ms1 analysis.
        ms1_peaks: [featid, mz, rt, intensity], featid will be None if "FeatureId"
            does not exist.
        """
        self.ms1_peaks = []
        self.user_cols_names = []
        with open(self.peaklist) as f:
            heads = f.readline()

            # add this in case peaklist file is separated by ';'
            self.separator = ','
            if ';' in heads:
                self.separator = ';'

            tokens = heads.strip().split(self.separator)
            index = -1
            featid_index = None
            for i in range(len(tokens)):
                if tokens[i].lower() == self.mz_col_name.lower():
                    index = i
                elif self.csv_id_col and tokens[i].lower(
                ) == self.csv_id_col.lower():
                    featid_index = i
                # if tokens[i].lower() == "scans":
                #     featid_index = i
                if tokens[i].lower() in ['mass',
                                         'mz']:  # backwards compatibility
                    index = i
                #     break
                self.user_cols_names.append(tokens[i])

            # if any sample names missing, use "Sample_*" to replace
            empty_sample_name_id = 0
            for i in range(index + 2, len(tokens)):
                if not tokens[i]:
                    tokens[i] = "Sample_" + str(empty_sample_name_id)
                    empty_sample_name_id += 1

            self.sample_names = tokens[index + 2:]

            for line in f:
                tokens_tuple = line.strip().split(self.separator, index + 2)
                featid = None
                if featid_index is not None:
                    featid = tokens_tuple[featid_index]
                mz = tokens_tuple[index]
                rt = float(tokens_tuple[index + 1])
                if self.rt_units == 'minutes':
                    rt *= 60.0
                samples = tokens_tuple[index + 2]
                # store (featid, mz,rt,intensity)

                # record user defined index columns before "mass" column in peaklist file
                try:
                    self.ms1_peaks.append((featid, float(mz), float(rt),
                                           samples, tokens_tuple[:index]))
                except:
                    print("Failed on line: ")
                    print(line)

        # sort them by mass
        self.ms1_peaks = sorted(self.ms1_peaks, key=lambda x: x[1])
        print("Loaded {} ms1 peaks from {}".format(len(self.ms1_peaks),
                                                   self.peaklist))

    def process_peaklist(self, ms1, ms2, metadata):
        """
        Read in peaklist .csv file.
        ("..., mass, RT, samplename_1, samplename_2, ..."), delimiter: '.
        Find the most suitable ms1 hit, then update ms1, ms2 metadata.
        """
        self._load_peak_list()
        ms1 = sorted(ms1, key=lambda x: x.mz)
        new_ms1_list = []
        new_ms2_list = []
        new_metadata = {}
        # ms1_mz = [x.mz for z in ms1]
        n_peaks_checked = 0

        # generate a dict (featid_ms1_dict)to store featid: ms1 pair
        # O(N) complexisity
        # build a dict (doc_ms1)for doc_name: ms1 pair first
        doc_ms1, featid_ms1_dict = {}, {}
        for el in ms1:
            doc_name = el.name
            doc_ms1[doc_name] = el

        for k, v in metadata.items():
            if self.id_field and (self.id_field.lower() in v):
                featid = v[self.id_field.lower()]
                featid_ms1_dict[featid] = doc_ms1[k]

        # build ms1_ms2 dict, to make searching O(1) in the following loop
        # key: ms1 object
        # value: list of ms2
        ms1_ms2_dict = {}
        for el in ms2:
            ms1_ms2_dict.setdefault(el[3], [])
            ms1_ms2_dict[el[3]].append(el)

        if self.id_field and self.csv_id_col:  # if the IDs are provided, we match by that
            print("IDs provided ({},{}), using them to match")
            match_by_id = True
        else:
            print("IDs not provided, matching on m/z, rt")
            match_by_id = False

        print("Matching peaks...")
        for n_peaks_checked, peak in enumerate(self.ms1_peaks):
            if n_peaks_checked % 500 == 0:
                print(n_peaks_checked)
            featid = peak[0]
            peak_mz = peak[1]
            peak_rt = peak[2]
            peak_intensity = None if self.separator in peak[3] else float(
                peak[3])
            user_cols = peak[4]

            # first check FeatureId matching
            # if featureId not exist, then do "mz/rt matching"
            old_ms1 = None

            if match_by_id:
                if featid is not None and featid in featid_ms1_dict:
                    old_ms1 = featid_ms1_dict[featid]
            else:
                min_mz = peak_mz - self.mz_tol * peak_mz / 1e6
                max_mz = peak_mz + self.mz_tol * peak_mz / 1e6
                min_rt = peak_rt - self.rt_tol
                max_rt = peak_rt + self.rt_tol

                ms1_hits = filter(
                    lambda x: x.mz >= min_mz and x.mz <= max_mz and x.rt >=
                    min_rt and x.rt <= max_rt, ms1)

                if len(ms1_hits) == 1:
                    # Found one hit, easy
                    old_ms1 = ms1_hits[0]
                elif len(ms1_hits) > 1:
                    # Find the one with the most intense MS2 peak
                    best_ms1 = None
                    best_intensity = 0.0
                    for frag_peak in ms2:
                        if frag_peak[3] in ms1_hits:
                            if frag_peak[2] > best_intensity:
                                best_intensity = frag_peak[2]
                                best_ms1 = frag_peak[3]
                    old_ms1 = best_ms1

            # Bug fix:
            # add these two lines to avoid the case that min_ms2_intensity has been set too high,
            # then most fragments will be removed, and we cannot find a hit for ms1, which will lead to bug:
            # AttributeError: 'NoneType' object has no attribute 'id'
            if not old_ms1:
                continue

            # make a new ms1 object
            new_ms1 = MS1(old_ms1.id, peak_mz, peak_rt, peak_intensity,
                          old_ms1.file_name, old_ms1.scan_number)
            new_ms1.name = old_ms1.name
            new_ms1_list.append(new_ms1)
            new_metadata[new_ms1.name] = metadata[old_ms1.name]

            # record user index columns before "mass" column in peaklist file into metadata
            new_metadata[new_ms1.name]['user_cols'] = zip(
                self.user_cols_names, user_cols)

            if self.separator in peak[3]:
                # print "process sample", str(peak[0]), str(peak[1])
                tokens = []
                for token in peak[3].split(self.separator):
                    try:
                        token = float(token)
                    except:
                        token = None
                    if token <= 0:
                        token = None
                    tokens.append(token)
                # tokens = [float(token) for token in peak[2].split(self.separator)]
                new_metadata[new_ms1.name]['intensities'] = dict(
                    zip(self.sample_names, tokens))

            # Delete the old one so it can't be picked again - removed this, maybe it's not a good idea?
            # pos = ms1.index(old_ms1)
            # del ms1[pos]

            # Change the reference in the ms2 objects to the new ms1 object

            # Use a dictionary outside the loop to replace the following method, O(N^2) => O(N)
            # ms2_objects = filter(lambda x: x[3] == old_ms1,ms2)
            ms2_objects = []
            if old_ms1 in ms1_ms2_dict:
                ms2_objects = ms1_ms2_dict[old_ms1]

            for frag_peak in ms2_objects:
                new_frag_peak = (frag_peak[0], peak_rt, frag_peak[2], new_ms1,
                                 frag_peak[4], frag_peak[5])
                new_ms2_list.append(new_frag_peak)

        # replace the ms1,ms2 and metadata with the new versions
        ms1 = new_ms1_list
        ms2 = new_ms2_list
        metadata = new_metadata
        print(f"Peaklist filtering results in {len(ms1)} documents")
        return ms1, ms2, metadata

    def filter_ms1_intensity(self, ms1, ms2, min_ms1_intensity=1e6):
        # Use filter function to simplify code
        print("Filtering MS1 on intensity")
        # Sometimes ms1 intensity could be None
        ms1 = filter(
            lambda x: False
            if x.intensity and x.intensity < min_ms1_intensity else True, ms1)
        print(f"{len(ms1)} MS1 remaining")
        ms2 = filter(lambda x: x[3] in set(ms1), ms2)
        print(f"{len(ms2)} MS2 remaining")
        return ms1, ms2

    def filter_ms2_intensity(self, ms2, min_ms2_intensity=1e6):
        print("Filtering MS2 on intensity")
        ms2 = filter(lambda x: x[2] >= min_ms2_intensity, ms2)
        print(f"{len(ms2)} MS2 remaining")
        return ms2

    def filter_ms1(self, ms1, ms2, mz_tol=0.5, rt_tol=16):
        print("Filtering MS1 to remove duplicates")
        # Filters the loaded ms1s to reduce the number of times that the same molecule has been fragmented

        # Sort the remaining ones by intensity
        ms1_by_intensity = sorted(ms1, key=lambda x: x.intensity, reverse=True)

        final_ms1_list = []
        final_ms2_list = []
        while True:
            if len(ms1_by_intensity) == 0:
                break
            # Take the highest intensity one, find things within the window and remove them
            current_ms1 = ms1_by_intensity[0]
            final_ms1_list.append(current_ms1)
            del ms1_by_intensity[0]

            current_mz = current_ms1.mz
            mz_err = mz_tol * 1.0 * current_mz / (1.0 * 1e6)
            min_mz = current_mz - mz_err
            max_mz = current_mz + mz_err

            min_rt = current_ms1.rt - rt_tol
            max_rt = current_ms1.rt + rt_tol

            # find things inside this region
            hits = filter(
                lambda x: x.mz > min_mz and x.mz < max_mz and x.rt > min_rt and
                x.rt < max_rt, ms1_by_intensity)
            for hit in hits:
                pos = ms1_by_intensity.index(hit)
                del ms1_by_intensity[pos]

        print(f"{len(final_ms1_list)} MS1 remaining")
        for m in ms2:
            if m[3] in final_ms1_list:
                final_ms2_list.append(m)

        print(f"{len(final_ms2_list)} MS2 remaining")
        return final_ms1_list, final_ms2_list

    def process_metadata(self, ms1, metadata):
        filtered_metadata = {}
        for m in ms1:
            if m.name in metadata:
                filtered_metadata[m.name] = metadata[m.name]
        metadata = filtered_metadata

        return metadata


# A class to load spectra that sit in MGF files
class LoadMGF(Loader):

    def load_spectra(self, input_set):
        ms1 = []
        ms2 = []
        metadata = {}
        ms2_id = 0
        ms1_id = 0
        for input_file in input_set:
            # Use built-in method to get file_name
            file_name = os.path.basename(input_file)
            with open(input_file) as f:
                temp_metadata = {}
                in_doc = False
                parentmass = None
                parentintensity = None
                parentrt = None
                for line in f:
                    rline = line.rstrip()
                    if not rline or rline == "BEGIN IONS":
                        continue
                    if rline == "END IONS":
                        # finished doc, time to save
                        in_doc = False
                        temp_metadata = {}
                        parentmass = None
                        parentintensity = None
                        parentrt = None
                        new_ms1 = None
                    else:
                        if "=" in rline:
                            key, val = rline.split("=", 1)
                            key = key.lower()

                            if len(val) == 0:
                                continue

                            featid = None
                            if key in ["featureid", "feature_id"]:
                                featid = val
                                temp_metadata['featid'] = val

                            elif key == "rtinseconds":
                                # val = float(val) if isinstance(val, float) else None
                                try:
                                    val = float(val)
                                except:
                                    val = None
                                temp_metadata['parentrt'] = val
                                parentrt = val

                            elif key == "pepmass":
                                # only mass exists
                                if " " not in val:
                                    temp_metadata['precursormass'] = float(val)
                                    temp_metadata['parentintensity'] = None
                                    parentmass = float(val)
                                    parentintensity = None

                                # mass and intensity exist
                                else:
                                    parentmass, parentintensity = val.split(
                                        ' ', 1)
                                    parentmass = float(parentmass)
                                    parentintensity = float(parentintensity)
                                    temp_metadata['precursormass'] = parentmass
                                    temp_metadata[
                                        'parentintensity'] = parentintensity

                            else:
                                temp_metadata[key] = val
                        else:
                            if 'mslevel' in temp_metadata and temp_metadata[
                                    'mslevel'] == '1':
                                continue

                            if not in_doc:
                                in_doc = True
                                # Following corrects parentmass according to charge
                                # if charge is known. This should lead to better computation of neutral losses
                                single_charge_precursor_mass = temp_metadata[
                                    'precursormass']
                                precursor_mass = temp_metadata['precursormass']
                                parent_mass = temp_metadata['precursormass']

                                str_charge = temp_metadata.get('charge', '1')
                                int_charge = self._interpret_charge(str_charge)

                                parent_mass, single_charge_precursor_mass = self._ion_masses(
                                    precursor_mass, int_charge)

                                temp_metadata['parentmass'] = parent_mass
                                temp_metadata[
                                    'singlechargeprecursormass'] = single_charge_precursor_mass
                                temp_metadata['charge'] = int_charge

                                new_ms1 = MS1(ms1_id,
                                              precursor_mass,
                                              parentrt,
                                              parentintensity,
                                              file_name,
                                              single_charge_precursor_mass=
                                              single_charge_precursor_mass)
                                ms1_id += 1

                                doc_name = temp_metadata.get(
                                    self.name_field.lower(), None)
                                if not doc_name:
                                    if 'name' in temp_metadata:
                                        doc_name = temp_metadata['name']
                                    else:
                                        doc_name = f'document_{ms1_id}'
                                metadata[doc_name] = temp_metadata.copy()
                                # TODO this overrides the original format of MS1.name attribute?
                                new_ms1.name = doc_name
                                ms1.append(new_ms1)

                            tokens = rline.split()
                            if len(tokens) == 2:
                                mz = float(tokens[0])
                                intensity = float(tokens[1])
                                if intensity != 0.0:
                                    ms2.append((mz, 0.0, intensity, new_ms1,
                                                file_name, float(ms2_id)))
                                    ms2_id += 1

        # add ms1, ms2 intensity filtering for msp input
        if self.min_ms1_intensity > 0.0:
            ms1, ms2 = self.filter_ms1_intensity(
                ms1, ms2, min_ms1_intensity=self.min_ms1_intensity)
        if self.min_ms2_intensity > 0.0:
            ms2 = self.filter_ms2_intensity(
                ms2, min_ms2_intensity=self.min_ms2_intensity)

        if self.peaklist:
            ms1, ms2, metadata = self.process_peaklist(ms1, ms2, metadata)

        # Chop out filtered docs from metadata
        metadata = self.process_metadata(ms1, metadata)

        return ms1, ms2, metadata
