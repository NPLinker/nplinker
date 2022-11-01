import json


class MibigMetadata():

    def __init__(self, file):
        """To represent the MIBiG BGC metadata (json file)

        Args:
            file(str): Path to the json file of MIBiG BGC metadata

        Examples:
            >>> metadata = MibigMetadata("/data/BGC0000001.json")
        """
        self.file = file
        with open(self.file, "rb") as f:
            self.metadata = json.load(f)

    @property
    def mibig_accession(self):
        """Get value of metadata item 'mibig_accession'"""
        return MibigMetadata.get_mibig_accession(self.metadata)

    @property
    def biosyn_class(self):
        """Get value of metadata item 'biosyn_class'"""
        return MibigMetadata.get_biosyn_class(self.metadata)

    @staticmethod
    def get_mibig_accession(metadata):
        """Get value of metadata item 'mibig_accession'"""
        if 'general_params' in metadata:
            accession = metadata['general_params']['mibig_accession']
        else:  # version≥2.0
            accession = metadata['cluster']['mibig_accession']
        return accession

    @staticmethod
    def get_biosyn_class(metadata):
        """Get value of metadata item 'mibig_accession'"""
        if 'general_params' in metadata:
            biosyn_class = metadata['general_params']['biosyn_class'][0]
        else:  # version≥2.0
            biosyn_class = metadata['cluster']['biosyn_class'][0]
        return biosyn_class
