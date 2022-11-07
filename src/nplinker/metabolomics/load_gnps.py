import csv
import os
import re
from typing import Literal, Optional
from nplinker.logconfig import LogConfig
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain


logger = LogConfig.getLogger(__file__)

# compile a regex for matching .mzXML and .mzML strings
RE_MZML_MZXML = re.compile('.mzXML|.mzML')

#
# methods for parsing metabolomics data files
#

GNPS_FORMAT_UNKNOWN = 'unknown'
GNPS_FORMAT_OLD_ALLFILES = 'allfiles'
GNPS_FORMAT_OLD_UNIQUEFILES = 'uniquefiles'
GNPS_FORMAT_NEW_FBMN = 'fbmn'


def _get_headers(filename: str, delimiters=['\t', ',']) -> list[str]:
    """Function to read headers from tab or comma separated table.

    Args:
        filename(str): Path to the file to read the header from.
        delimiters (list, optional): List of delimiters to consider. Defaults to ['\t', ','].

    Returns:
        list[str]: Columns names in header.
    """
    headers = None
    with open(filename) as f:
        headers = f.readline()
        for dl in delimiters:
            if len(headers.split(dl)) < 2:
                continue
            headers = headers.split(dl)
            break
    return headers


def identify_gnps_format(filename: str, has_quant_table: bool) -> Literal['unknown', 'allfiles', 'uniquefiles', 'fbmn']:
    """Peek GNPS file format for given file.

    TODO: #89 This should be rewritten to actually return the format always based on only the file and not include the quant table in it.

     Args:
        filename(str): Path to the file to peek the format for.
        has_quant_table (bool): If a quant table is present, do return GNPS_FORMAT_NEW_FBMN.

    Returns:
        Literal['unknown', 'allfiles', 'uniquefiles', 'fbmn']: Constant, one of 
        [GNPS_FORMAT_UNKNOWN, GNPS_FORMAT_OLD_ALLFILES, GNPS_FORMAT_OLD_UNIQUEFILES, GNPS_FORMAT_NEW_FBMN].

    Examples:
        >>> 
        """

    headers = _get_headers(filename)

    if headers is None:
        return GNPS_FORMAT_UNKNOWN

    # first, check for AllFiles
    if 'AllFiles' in headers:
        # this should be an old-style dataset like Crusemann, with a single .tsv file
        # containing all the necessary info. The AllFiles column should contain pairs
        # of mzXML filenames and scan numbers in this format:
        #   filename1:scannumber1###filename2:scannumber2###...
        return GNPS_FORMAT_OLD_ALLFILES
    elif 'UniqueFileSources' in headers:
        # this is a slightly newer-style dataset, e.g. MSV000084771 on the platform
        # it still just has a single .tsv file, but AllFiles is apparently replaced
        # with a UniqueFileSources column. There is also a UniqueFileSourcesCount
        # column which just indicates the number of entries in the UniqueFileSources
        # column. If there are multiple entries the delimiter is a | character
        return GNPS_FORMAT_OLD_UNIQUEFILES
    elif has_quant_table:
        # if there is no AllFiles/UniqueFileSources, but we DO have a quantification
        # table file, that should indicate a new-style dataset like Carnegie
        # TODO check for the required header columns here too
        return GNPS_FORMAT_NEW_FBMN
    else:
        # if we don't match any of the above cases then it's not a recognised format
        return GNPS_FORMAT_UNKNOWN


def _messy_strain_naming_lookup(mzxml: str, strains: StrainCollection) -> Optional[Strain]:
    if mzxml in strains:
        # life is simple!
        return strains.lookup(mzxml)

    # 1. knock off the .mzXML and try again (using splitext should also handle
    # other extensions like .mzML here)
    mzxml_no_ext = os.path.splitext(mzxml)[0]
    if mzxml_no_ext in strains:
        return strains.lookup(mzxml_no_ext)

    # 2. if that doesn't work, try using everything up to the first -/_
    underscore_index = mzxml_no_ext.find('_')
    hyphen_index = mzxml_no_ext.find('-')
    mzxml_trunc_underscore = mzxml_no_ext if underscore_index == -1 else mzxml_no_ext[:underscore_index]
    mzxml_trunc_hyphen = mzxml_no_ext if hyphen_index == -1 else mzxml_no_ext[:hyphen_index]
    
    if underscore_index != -1 and mzxml_trunc_underscore in strains:
        return strains.lookup(mzxml_trunc_underscore)
    if hyphen_index != -1 and mzxml_trunc_hyphen in strains:
        return strains.lookup(mzxml_trunc_hyphen)

    # 3. in the case of original Crusemann dataset, many of the filenames seem to
    # match up to real strains with the initial "CN" missing ???
    for mzxml_trunc in [mzxml_trunc_hyphen, mzxml_trunc_underscore]:
        if 'CN' + mzxml_trunc in strains:
            return strains.lookup('CN' + mzxml_trunc)

    # give up
    return None


def _md_convert(val):
    """Try to convert raw metadata values from text to integer, then float if that fails"""
    try:
        return int(val)
    except (ValueError, TypeError):
        try:
            return float(val)
        except (ValueError, TypeError):
            if val.lower() == 'n/a':
                return None
    return val


def _parse_mzxml_header(hdr, strains, md_table, ext_metadata_parsing):
    # given a column header from the quantification_table file produced by GNPS,
    # this method checks if it's one that contains peak information for a given
    # strain label by searching for the 'Peak area' string.
    #    e.g. 'KRD168_ISP3.mzML Peak area'
    # it then extracts the strain label by taking the text up to the '.mzML' part
    # and attempts to match this to the set of strains provided by the user in
    # the strain_mappings file. finally it also tries to extract the growth
    # medium label as given in the metadata_table file, again using the strain
    # label which should match between the two files.

    # assume everything up to '.mz' is the identifier/label of this strain
    strain_name = hdr[:hdr.index('.mz')]
    growth_medium = None

    # check if the strain label exists in the set of strain mappings the user
    # has provided
    if strain_name not in strains:
        # if this check fails, it could mean a missing strain ID mapping. this isn't
        # fatal but will produce a warning and an entry added to the file
        # unknown_strains_met.csv in the dataset folder.
        #
        # however there could be some identifiers which are not strains and
        # these should just be ignored.
        #
        # if we have a metadata table file, and the parsed strain name does NOT
        # appear in it, this indicates it's not a strain
        if md_table is not None and strain_name not in md_table:
            # we can ignore this completely
            return (None, None, False)

        # if the strain is in the table, then it means the user probably hasn't given us
        # a valid mapping for it, so a warning will be displayed and the name added
        # to the unknown_strains.csv file later on. however if the config file option
        # extended_metadata_table_parsing has been enabled (the ext_metadata_parsing)
        # parameter in this method), it means we should take the content of the
        # ATTRIBUTE_Strain column from the table and use that as the strain, adding
        # an extra alias for it too. Additionally the ATTRIBUTE_Medium column content
        # should be recorded as the growth medium for the strain
        if md_table is not None and strain_name in md_table and ext_metadata_parsing:
            growth_medium = md_table[strain_name]['growthmedium']
            strain_col_content = md_table[strain_name]['strain']

            if strain_col_content in strains:
                strain = strains.lookup(strain_col_content)
                strain.add_alias(strain_name)
                # this merges the set of aliases back into the internal
                # lookup table in the StrainCollection
                strains.add(strain)
            else:
                # if this still fails it's an unknown strain
                logger.warning(
                    'Unknown strain identifier: {} (parsed from {})'.format(
                        strain_name, hdr))
                return (strain_name, growth_medium, True)
        else:
            # emit a warning message to indicate that the user needs to add this to
            # their strain_mappings file
            logger.warning(
                'Unknown strain identifier: {} (parsed from {})'.format(
                    strain_name, hdr))

    return (strain_name, growth_medium, strain_name not in strains)


def _load_clusterinfo_old(gnps_format, strains, filename, spec_dict):
    # each line of this file represents a metabolite.
    # columns representing strain IDs are *ignored* here in favour of parsing
    # .mz(X)ML filenames from either the AllFiles or UniqueFileSources column.
    # both of these list the .mz(X)ML files the molecule was found in (plus the scan
    # number in the file in the AllFiles case)
    unknown_strains = {}
    with open(filename) as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        clu_index_index = headers.index('cluster index')
        if gnps_format == GNPS_FORMAT_OLD_ALLFILES:
            mzxml_index = headers.index('AllFiles')
        elif gnps_format == GNPS_FORMAT_OLD_UNIQUEFILES:
            mzxml_index = headers.index('UniqueFileSources')
        else:
            raise Exception(f'Unexpected GNPS format {gnps_format}')

        for line in reader:
            # get the values of the important columns
            clu_index = int(line[clu_index_index])
            if gnps_format == GNPS_FORMAT_OLD_UNIQUEFILES:
                mzxmls = line[mzxml_index].split('|')
            else:
                mzxmls = line[mzxml_index].split('###')

            metadata = {'cluster_index': clu_index, 'files': {}}
            seen_files = set()

            for data in mzxmls:
                # TODO ok to ignore scan number if available?
                mzxml = data if gnps_format == GNPS_FORMAT_OLD_UNIQUEFILES else data.split(
                    ':')[0]

                # TODO is this OK/sensible?
                if mzxml in seen_files:
                    continue
                seen_files.add(mzxml)

                # add to the list of files for this molecule
                metadata['files'][mzxml] = mzxml

                # should have a strain alias for the mxXML file to lookup here (in theory at least)
                strain = _messy_strain_naming_lookup(mzxml, strains)
                if strain is None:
                    # TODO: how to handle this? can't just give up, simply warn?
                    if mzxml not in unknown_strains:
                        logger.warning(
                            'Unknown strain: {} for cluster index {}'.format(
                                mzxml, clu_index))
                        unknown_strains[mzxml] = 1
                    else:
                        unknown_strains[mzxml] += 1
                # else:
                #     print('{} ===> {}'.format(mzxml, strain))

                if strain is not None:
                    # TODO this need fixed somehow (missing growth medium info)
                    spec_dict[clu_index].add_strain(strain, None, 1)

                # update metadata on Spectrum object
                spec_dict[clu_index].metadata.update(metadata)

    if len(unknown_strains) > 0:
        logger.warning(
            '{} unknown strains were detected a total of {} times'.format(
                len(unknown_strains), sum(unknown_strains.values())))

    return unknown_strains


def _parse_metadata_table(filename):
    """Parses the metadata table file from GNPS"""
    if filename is None:
        return None

    table = {}
    with open(filename) as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)

        # check the format is as expected
        if len(headers) < 3 or headers[0] != 'filename':
            # expecting at least 3 columns with the first being 'filename'
            logger.warning(
                'Unrecognised metadata table format in file "{}"'.format(
                    filename))
            return None

        # find the column numbers we're interested in (can't rely on both of these being present)
        # ATTRIBUTE_SampleType, _Strain, _Medium, _Organism
        sampletype_col, medium_col, strain_col = -1, -1, -1
        col_names = {
            'sampletype': 'ATTRIBUTE_SampleType',
            'growthmedium': 'ATTRIBUTE_Medium',
            'strain': 'ATTRIBUTE_Strain'
        }

        # also: ATTRIBUTE_Strain might be useful, there's also ATTRIBUTE_Organism
        for i, hdr in enumerate(headers):
            if hdr == col_names['sampletype']:
                sampletype_col = i
            if hdr == col_names['growthmedium']:
                medium_col = i
            if hdr == col_names['strain']:
                strain_col = i

        for col, name in [(sampletype_col, col_names['sampletype']),
                          (medium_col, col_names['growthmedium']),
                          (strain_col, col_names['strain'])]:
            if col == -1:
                logger.warning('No {} column in file "{}"'.format(
                    name, filename))

        for line in reader:
            # only non-BLANK entries
            if sampletype_col != -1 and line[sampletype_col] != 'BLANK':
                # want to strip out the .mz(X)ML extension here to use as the key,
                # then record the available columns (assume SampleType is always
                # available)
                data = {
                    'sampletype': line[sampletype_col],
                    'growthmedium': None,
                    'strain': None
                }
                if medium_col != -1:
                    data['growthmedium'] = line[medium_col]
                if strain_col != -1:
                    data['strain'] = line[strain_col]
                table[RE_MZML_MZXML.sub('', line[0])] = data

    return table


def _load_clusterinfo_fbmn(strains, nodes_file, extra_nodes_file,
                           md_table_file, spec_dict, ext_metadata_parsing):
    spec_info = {}

    # parse metadata table if available
    md_table = _parse_metadata_table(md_table_file)

    # combine the information in the nodes_file (clusterinfo_summary folder) and
    # the extra_nodes_file (quantification_table folder), indexed by the "cluster index"
    # and "row ID" fields respectively to link the rows
    with open(nodes_file) as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        ci_index = headers.index('cluster index')

        # take each line and generate a dict storing each header:column value
        # and insert that dict into the spec_info dict
        for line in reader:
            tmp = {}
            for i, v in enumerate(line):
                tmp[headers[i]] = v
            spec_info[int(line[ci_index])] = tmp

    with open(extra_nodes_file) as f:
        reader = csv.reader(f, delimiter=',')
        headers = next(reader)
        ri_index = headers.index('row ID')

        # do almost the same thing again but a) match the "cluster index" from the
        # nodes_file to the "row ID" from this file, and update the per-row dict
        # with the extra columns from this file
        for line in reader:
            ri = int(line[ri_index])
            tmp = {}
            for i, v in enumerate(line):
                tmp[headers[i]] = v
            spec_info[ri].update(tmp)

    logger.info('Merged nodes data (new-style), total lines = {}'.format(
        len(spec_info)))

    unknown_strains = {}

    # TODO: at the moment this iterates over each spectrum in the outer loop, and
    # over the full set of columns in the inner loop. probably makes more sense to
    # swap that around, or at least to avoid calling parse_mzxml_header repeatedly
    # for the same column headers over and over!

    # for each spectrum
    for spec_id, spec_data in spec_info.items():
        # look up the Spectrum object with this ID (spec_dict is sourced from the MGF file)
        spectrum = spec_dict[spec_id]

        # TODO this should probably be handled by parse_metadata_table
        if 'ATTRIBUTE_SampleType' in spec_data:
            st = spec_data['ATTRIBUTE_SampleType'].split(',')
            spectrum.metadata['ATTRIBUTE_SampleType'] = st

        # TODO better way of filtering/converting all this stuff down to what's relevant?
        # could search for each strain ID in column title but would be slower?
        #
        # iterate over the column values for the spectrum
        for k, v in spec_data.items():
            # if the value is a "0", ignore immediately as we only care about nonzero peak intensities
            if v == "0":
                continue

            # ignore any non-"Peak area" columns, want "<strain ID>.mz(X)ML Peak area"
            if k.find(' Peak area') == -1:
                # record the value as an entry in the metadata dict in case it's useful
                # for some other purpose
                spectrum.metadata[k] = v
                continue

            # TODO this will probably need updating for platform data (should already
            # have a set of strain names...)
            #
            # given a "<strain ID>.mz(X)ML Peak area" column value, attempt to lookup the
            # matching Strain object from the set of mappings we have (optionally doing
            # some further stuff involving the metadata table)
            (strain_name, growth_medium,
             is_unknown) = _parse_mzxml_header(k, strains, md_table,
                                               ext_metadata_parsing)

            # if the strain name couldn't be recognised at all, ignore it
            if strain_name is None:
                continue

            # if the name seems valid (appears in the metadata table), record it as unknown
            if is_unknown:
                unknown_strains[k[:k.index('.mz')]] = spec_id
                continue

            # create a new strain object if the intensity value is a float > 0
            v = _md_convert(v)
            if strain_name in strains and isinstance(v, float) and v > 0:
                # find the strain object, and add the growth medium + intensity to it. the
                # growth medium will only be set if extended_metadata_table_parsing is
                # enabled in the config file and the metadata table file contains that info
                strain = strains.lookup(strain_name)
                spectrum.add_strain(strain, growth_medium, v)

            # record this as an entry in the metadata dict as well
            spectrum.metadata[k] = v

    return spec_info, unknown_strains


def load_gnps(strains, nodes_file, quant_table_file, metadata_table_file,
              ext_metadata_parsing, spec_dict):
    gnps_format = identify_gnps_format(nodes_file, quant_table_file
                                        is not None)
    logger.debug('Nodes_file: {}, quant_table_exists?: {}'.format(
        nodes_file, quant_table_file is None))
    if gnps_format == GNPS_FORMAT_UNKNOWN:
        raise Exception('Unknown/unsupported GNPS data format')

    # now things depend on the dataset format
    # if we don't have a quantification table, must be older-style dataset (or unknown format)
    if gnps_format != GNPS_FORMAT_NEW_FBMN and quant_table_file is None:
        logger.info('Found older-style GNPS dataset, no quantification table')
        unknown_strains = _load_clusterinfo_old(gnps_format, strains,
                                                nodes_file, spec_dict)
    else:
        logger.info('quantification table exists, new-style GNPS dataset')
        _, unknown_strains = _load_clusterinfo_fbmn(strains, nodes_file,
                                                    quant_table_file,
                                                    metadata_table_file,
                                                    spec_dict,
                                                    ext_metadata_parsing)

    return unknown_strains
