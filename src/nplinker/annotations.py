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

import csv
import os
from .logconfig import LogConfig


logger = LogConfig.getLogger(__file__)

GNPS_KEY = 'gnps'
GNPS_URL_FORMAT = 'https://metabolomics-usi.ucsd.edu/{}/?usi=mzspec:GNPSLIBRARY:{}'
GNPS_INDEX_COLUMN = '#Scan#'
GNPS_DATA_COLUMNS = ['Compound_Name', 'Organism', 'MQScore', 'SpectrumID']


def headers_match_gnps(headers):
    for k in GNPS_DATA_COLUMNS:
        if k not in headers:
            return False
    return True


# url_type: spectrum, png, svg, json
def gnps_url(id, url_type='spectrum'):
    return GNPS_URL_FORMAT.format(url_type, id)


def create_gnps_annotation(spec, gnps_anno):
    # also insert useful URLs
    for t in ['png', 'json', 'svg', 'spectrum']:
        gnps_anno[f'{t}_url'] = GNPS_URL_FORMAT.format(
            t, gnps_anno['SpectrumID'])

    if GNPS_KEY in spec.annotations:
        # TODO is this actually an error or can it happen normally?
        raise Exception('Multiple GNPS annotations for Spectrum {}!'.format(
            spec.spectrum_id))

    spec.set_annotations(GNPS_KEY, [gnps_anno])
    # shortcut, useful for rosetta code
    spec.gnps_id = gnps_anno['SpectrumID']


def load_annotations(root, config, spectra, spec_dict):
    if not os.path.exists(root):
        logger.debug(f'Annotation directory not found ({root})')
        return spectra

    ac = {}
    if os.path.exists(config):
        # parse annotations config file if it exists
        with open(config) as f:
            rdr = csv.reader(f, delimiter='\t')
            for row in rdr:
                # expecting 3 columns: filename, index col name, data col name(s)
                if len(row) != 3:
                    logger.warning(
                        'Malformed line in annotations configuration: {}'.
                        format(row))
                    continue
                # record the columns with filename as key
                data_cols = row[2].split(',')
                if len(data_cols) == 0:
                    logger.warning(
                        'No data columns selected in annotation file {}, skipping it!'
                        .format(row[0]))
                    continue
                ac[row[0]] = (row[1], data_cols)

    logger.debug(f'Parsed {len(ac)} annotations configuration entries')

    annotation_files = []
    for f in os.listdir(root):
        if f == os.path.split(config)[1] or not f.endswith('.tsv'):
            continue
        annotation_files.append(os.path.join(root, f))

    logger.debug('Found {} annotations .tsv files in {}'.format(
        len(annotation_files), root))
    for af in annotation_files:
        with open(af) as f:
            rdr = csv.reader(f, delimiter='\t')
            headers = next(rdr)
            filename = os.path.split(af)[1]

            if headers_match_gnps(headers):
                # assume this is our main annotations file
                logger.debug(f'Parsing GNPS annotations from {af}')

                scans_index = headers.index(GNPS_INDEX_COLUMN)

                # each line should be a different spec ID here
                for line in rdr:
                    # read the scan ID column and get the corresponding Spectrum object
                    scan_id = int(line[scans_index])
                    if scan_id not in spec_dict:
                        logger.warning(
                            'Unknown spectrum ID found in GNPS annotation file (ID={})'
                            .format(scan_id))
                        continue

                    spec = spec_dict[scan_id]
                    data_cols = set(GNPS_DATA_COLUMNS)
                    # merge in any extra columns the user has provided
                    if filename in ac:
                        data_cols.update(ac[filename][1])

                    data = {}
                    for dc in data_cols:
                        if dc not in headers:
                            logger.warning(
                                f'Column lookup failed for "{dc}"')
                            continue
                        data[dc] = line[headers.index(dc)]

                    create_gnps_annotation(spec, data)
            else:
                logger.debug(f'Parsing general annotations from {af}')
                # this is a general annotations file, so rely purely on the user-provided columns
                if filename not in ac:
                    logger.warning(
                        'Failed to parse annotations from "{}", no config info supplied in {}'
                        .format(filename, config))
                    continue

                index_col, data_cols = ac[filename]
                if index_col not in headers:
                    raise Exception(
                        'Annotation index column "{}" not found in file "{}"!'.
                        format(index_col, filename))

                spec_id_index = headers.index(index_col)

                # note that might have multiple lines for the same spec ID!
                spec_annotations = {}
                for line in rdr:
                    scan_id = int(line[spec_id_index])
                    if scan_id not in spec_dict:
                        logger.warning(
                            'Unknown spectrum ID found in annotation file "{}", ID is "{}"'
                            .format(filename, scan_id))
                        continue

                    spec = spec_dict[scan_id]
                    if spec not in spec_annotations:
                        spec_annotations[spec] = []

                    data = {}
                    for dc in data_cols:
                        if dc not in headers:
                            logger.warning(
                                f'Column lookup failed for "{dc}"')
                            continue
                        data[dc] = line[headers.index(dc)]
                    spec_annotations[spec].append(data)

                logger.debug('Adding annotation data to {} spectra'.format(
                    len(spec_annotations)))
                for spec, anno in spec_annotations.items():
                    spec.set_annotations(filename, anno)

    return spectra
