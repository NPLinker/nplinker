import csv
from typing import TextIO
from nplinker.logconfig import LogConfig
from nplinker.metabolomics.IFileMappingLoader import FileMappingLoaderBase
from nplinker.metabolomics.load_gnps import GNPS_FORMAT_NEW_FBMN
from nplinker.metabolomics.load_gnps import GNPS_FORMAT_OLD_ALLFILES
from nplinker.metabolomics.load_gnps import identify_gnps_format
from nplinker.utils import find_delimiter

logger = LogConfig.getLogger(__file__)

FILE_IDENTIFIER_FBMN = " Peak area"


class GNPSFileMappingLoader(FileMappingLoaderBase):

    def __init__(self, filename: str):
        self._filename: str = filename
        self._mapping: dict[int, list[str]] = {}
        self._gnps_format = identify_gnps_format(filename, False)

        if self._gnps_format is GNPS_FORMAT_OLD_ALLFILES:
            self._load_mapping_allfiles()
        elif self._gnps_format is GNPS_FORMAT_NEW_FBMN:
            self._load_mapping_fbmn()
        else:
            raise NotImplementedError(
                "%{gnps_format} reading not implemented.")

    def mapping(self) -> dict[int, list[str]]:
        """Return mapping from spectrum id to files in which this spectrum occurs.

        Returns:
            dict[int, list[str]]: Mapping from spectrum id to names of all files in which this spectrum occurs.
        """
        return self._mapping

    def _load_mapping_allfiles(self):
        """ Load mapping for GNPS 'AllFiles' style files. """
        with open(self._filename, mode='rt', encoding='utf-8') as file:
            reader = self._get_dict_reader(file)

            for row in reader:
                spectrum_id = int(row["cluster index"])

                occurrences = row["AllFiles"].split("###")  # split by '###'
                occurrences.pop()  # remove last empty entry
                # separate the scan position from the files
                samples = [x.split(':')[0] for x in occurrences]

                self._mapping[spectrum_id] = samples

    def _get_dict_reader(self, file: TextIO) -> csv.DictReader:
        """Get a dict reader with matching delimiter for the passed file.

        Args:
            file(TextIOWrapper): File for which to get the reader.

        Returns:
            csv.DictReader: Reader for dict style table access.
        """
        delimiter = find_delimiter(self._filename)
        reader = csv.reader(file, delimiter=delimiter)
        header: list[str] = next(reader)
        reader = csv.DictReader(file, header, delimiter=delimiter)
        return reader

    def _load_mapping_fbmn(self):
        """ Load mapping for GNPS 'fbmn' style files. """
        with open(self._filename, mode='rt', encoding='utf-8') as file:
            reader = self._get_dict_reader(file)

            for row in reader:
                spectrum_id = int(row["row ID"])

                if self._mapping.get(spectrum_id) is not None:
                    logger.warning("Found duplicated row ID: %{spectrum_id}")

                samples = []

                for col in row:
                    if FILE_IDENTIFIER_FBMN in col:
                        if float(row[col]) > 0:
                            samples.append(col.strip(FILE_IDENTIFIER_FBMN))

                self._mapping[spectrum_id] = samples