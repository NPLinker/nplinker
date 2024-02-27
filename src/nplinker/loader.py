import os
from pathlib import Path
from deprecated import deprecated
from nplinker.class_info.chem_classes import ChemClassPredictions
from nplinker.class_info.class_matches import ClassMatches
from nplinker.class_info.runcanopus import run_canopus
from nplinker.config import config
from nplinker.genomics import add_bgc_to_gcf
from nplinker.genomics import add_strain_to_bgc
from nplinker.genomics import generate_mappings_genome_id_bgc_id
from nplinker.genomics import get_mibig_from_gcf
from nplinker.genomics.antismash import AntismashBGCLoader
from nplinker.genomics.bigscape import BigscapeGCFLoader
from nplinker.genomics.mibig import MibigLoader
from nplinker.globals import GENOME_BGC_MAPPINGS_FILENAME
from nplinker.globals import GENOME_STATUS_FILENAME
from nplinker.globals import GNPS_FILE_MAPPINGS_FILENAME
from nplinker.globals import STRAIN_MAPPINGS_FILENAME
from nplinker.logconfig import LogConfig
from nplinker.metabolomics import add_annotation_to_spectrum
from nplinker.metabolomics import add_spectrum_to_mf
from nplinker.metabolomics import add_strains_to_spectrum
from nplinker.metabolomics.gnps import GNPSAnnotationLoader
from nplinker.metabolomics.gnps import GNPSMolecularFamilyLoader
from nplinker.metabolomics.gnps import GNPSSpectrumLoader
from nplinker.pairedomics.strain_mappings_generator import podp_generate_strain_mappings
from nplinker.strain_collection import StrainCollection
from nplinker.strain_loader import load_user_strains


try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

logger = LogConfig.getLogger(__name__)

NPLINKER_APP_DATA_DIR = files("nplinker").joinpath("data")


class DatasetLoader:
    RUN_CANOPUS_DEFAULT = False
    EXTRA_CANOPUS_PARAMS_DEFAULT = "--maxmz 600 formula zodiac structure canopus"

    OR_STRAINS = "strain_mappings_file"
    # misc files
    OR_INCLUDE_STRAINS = "include_strains_file"
    # class predictions
    OR_CANOPUS = "canopus_dir"
    OR_MOLNETENHANCER = "molnetenhancer_dir"

    def __init__(self):
        # load the config data
        self._config_docker = config.get("docker", {})
        self._config_webapp = config.get("webapp", {})
        self._config_overrides = config.dataset.get("overrides", {})
        # set private attributes
        self._root = Path(config.dataset["root"])
        self._platform_id = config.dataset["platform_id"]

        # set public attributes
        self.bgcs, self.gcfs, self.spectra, self.molfams = [], [], [], []
        self.mibig_bgcs = []
        self.mibig_strains_in_use = StrainCollection()
        self.product_types = []
        self.strains = StrainCollection()
        self.webapp_scoring_cutoff = self._config_webapp.get(
            "tables_metcalf_threshold", self.TABLES_CUTOFF_DEFAULT
        )
        self.class_matches = None
        self.chem_classes = None

    def validate(self):
        """Download data and build paths for local data."""
        self._init_paths()

    def generate_strain_mappings(self):
        generate_mappings_genome_id_bgc_id(self._root / "antismash")

        podp_project_json_file = self._root.parent.parent / (self._platform_id + ".json")
        genome_status_json_file = (
            self._root.parent.parent / "downloads" / self._platform_id / GENOME_STATUS_FILENAME
        )
        genome_bgc_mappings_file = self._root / "antismash" / GENOME_BGC_MAPPINGS_FILENAME
        gnps_file_mapping_tsv_file = self._root / GNPS_FILE_MAPPINGS_FILENAME

        podp_generate_strain_mappings(
            podp_project_json_file,
            genome_status_json_file,
            genome_bgc_mappings_file,
            gnps_file_mapping_tsv_file,
            self.strain_mappings_file,
        )

    def load(self):
        if not self._load_strain_mappings():
            return False

        if not self._load_metabolomics():
            return False

        if not self._load_genomics():
            return False

        # set self.strains with all strains from input plus mibig strains in use
        self.strains = self.strains + self.mibig_strains_in_use

        if len(self.strains) == 0:
            raise Exception(f"Failed to find *ANY* strains, missing {STRAIN_MAPPINGS_FILENAME}?")

        return True

    def _init_paths(self):
        # 1. strain mapping are used for everything else so
        self.strain_mappings_file = self._config_overrides.get(self.OR_STRAINS) or os.path.join(
            self._root, STRAIN_MAPPINGS_FILENAME
        )

        # 14. MISC: <root>/include_strains.csv / include_strains_file=<override>
        self.include_strains_file = self._config_overrides.get(
            self.OR_INCLUDE_STRAINS
        ) or os.path.join(self._root, "include_strains.csv")

        # 15. CLASS: <root>/canopus / canopus_dir=<override>
        self.canopus_dir = self._config_overrides.get(self.OR_CANOPUS) or os.path.join(
            self._root, "canopus"
        )

        # 15. CLASS: <root>/canopus / canopus_dir=<override>
        self.molnetenhancer_dir = self._config_overrides.get(
            self.OR_MOLNETENHANCER
        ) or os.path.join(self._root, "molnetenhancer")

    def _load_strain_mappings(self):
        # 1. load strain mappings
        sc = StrainCollection.read_json(self.strain_mappings_file)
        for strain in sc:
            self.strains.add(strain)
        logger.info("Loaded {} non-MiBIG Strain objects".format(len(self.strains)))

        # 2. filter user specificied strains (remove all that are not specified by user).
        # It's not allowed to specify empty list of strains, otherwise validation will fail.
        if os.path.exists(self.include_strains_file):
            logger.info(f"Loading user specified strains from file {self.include_strains_file}.")
            user_strains = load_user_strains(self.include_strains_file)
            logger.info(f"Loaded {len(user_strains)} user specified strains.")
            self.strains.filter(user_strains)

        logger.info("Loaded {} Strain objects in total".format(len(self.strains)))
        return True

    def _load_metabolomics(self):
        """Loads metabolomics data to Spectrum and MolecularFamily objects.

        The attribute of `self.spectra` is set to the loaded Spectrum objects that have Strain
        objects added (i.e. `Spectrum.strains` updated). If a Spectrum object does not have Strain
        objects, it is not added to `self.spectra`.

        The attribute of `self.molfams` is set to the loaded MolecularFamily objects that have
        Strain objects added (i.e. `MolecularFamily._strains` updated). This means only Spectra
        objects with updated strains (i.e. `self.spectra`) can be added to MolecularFamily objects.
        """
        logger.debug("\nLoading metabolomics data starts...")

        # Step 1: load all Spectrum objects
        raw_spectra = GNPSSpectrumLoader(self.gnps_mgf_file).spectra
        # Step 2: load all GNPS annotations
        raw_annotations = GNPSAnnotationLoader(self.gnps_annotations_file).annotations
        # Step 3: load all MolecularFamily objects
        raw_molfams = GNPSMolecularFamilyLoader(self.gnps_edges_file).get_mfs(keep_singleton=False)

        # Step 4: add GNPS annotations to Spectrum.gnps_annotations
        add_annotation_to_spectrum(raw_annotations, raw_spectra)
        # Step 5: add strains to Spectrum.strains
        spectra_with_strains, _ = add_strains_to_spectrum(self.strains, raw_spectra)

        # Step 6: add Spectrum objects to MolecularFamily
        mf_with_spec, _, _ = add_spectrum_to_mf(spectra_with_strains, raw_molfams)

        # Step 7: set attributes of self.spectra and self.molfams with valid objects
        self.spectra = spectra_with_strains
        self.molfams = mf_with_spec

        logger.debug("Loading metabolomics data completed\n")
        return True

    def _load_genomics(self):
        """Loads genomics data to BGC and GCF objects.

        The attribute of `self.bgcs` is set to the loaded BGC objects that have the Strain object
        added (i.e. `BGC.strain` updated). If a BGC object does not have the Strain object, it is
        not added to `self.bgcs`. For MIBiG BGC objects, only those in use are added to `self.bgcs`.

        The attribute of `self.gcfs` is set to the loaded GCF objects that have the Strain objects
        added (i.e. `GCF._strains` updated). This means only BGC objects with updated Strain objects
        (i.e. `self.bgcs`) can be added to GCF objects.
        """
        logger.debug("\nLoading genomics data starts...")

        # Step 1: load antismash BGC objects & add strain info
        logger.debug("Parsing AntiSMASH directory...")
        antismash_bgcs = AntismashBGCLoader(self.antismash_dir).get_bgcs()
        antismash_bgcs_with_strain, _ = add_strain_to_bgc(self.strains, antismash_bgcs)

        # Step 2: load mibig BGC objects (having strain info)
        if self._use_mibig:
            self.mibig_bgcs = MibigLoader(self.mibig_json_dir).get_bgcs()

        # Step 3: get all BGC objects with strain info
        all_bgcs_with_strain = antismash_bgcs_with_strain + self.mibig_bgcs

        # Step 4: load all GCF objects
        # TODO: create a config for "bigscape_cluster_file" and discard "bigscape_dir" and "bigscape_cutoff"?
        bigscape_cluster_file = (
            Path(self.bigscape_dir) / "mix" / f"mix_clustering_c0.{self._bigscape_cutoff:02d}.tsv"
        )
        raw_gcfs = BigscapeGCFLoader(bigscape_cluster_file).get_gcfs()

        # Step 5: add BGC objects to GCF
        all_gcfs_with_bgc, _, _ = add_bgc_to_gcf(all_bgcs_with_strain, raw_gcfs)

        # Step 6: get mibig bgcs and strains in use from GCFs
        mibig_strains_in_use = StrainCollection()
        if self._use_mibig:
            mibig_bgcs_in_use, mibig_strains_in_use = get_mibig_from_gcf(all_gcfs_with_bgc)
        else:
            mibig_bgcs_in_use = []

        # Step 7: set attributes with valid objects
        self.bgcs = antismash_bgcs_with_strain + mibig_bgcs_in_use
        self.gcfs = all_gcfs_with_bgc
        self.mibig_strains_in_use = mibig_strains_in_use

        logger.debug("Loading genomics data completed\n")
        return True

    @deprecated(reason="To be refactored. It was used in the `self.load` method before.")
    def _load_class_info(self):
        """Load class match info (based on mibig) and chemical class predictions.

        Run CANOPUS if asked for. First sirius is run through docker, if this
        fails, it is run with a version present on the path.

        Return:
            True if everything completes
        """
        # load Class_matches with mibig info from data
        mibig_class_file = NPLINKER_APP_DATA_DIR.joinpath(
            "MIBiG2.0_compounds_with_AS_BGC_CF_NPC_classes.txt"
        )
        self.class_matches = ClassMatches(mibig_class_file)

        # run canopus if canopus_dir does not exist
        should_run_canopus = self._config_docker.get("run_canopus", self.RUN_CANOPUS_DEFAULT)
        extra_canopus_parameters = self._config_docker.get(
            "extra_canopus_parameters", self.EXTRA_CANOPUS_PARAMS_DEFAULT
        )
        if should_run_canopus:
            # don't run canopus when canopus dir exists already
            if not os.path.isdir(self.canopus_dir):
                logger.info(
                    'Running CANOPUS! extra_canopus_parameters="{}"'.format(
                        extra_canopus_parameters
                    )
                )
                try:
                    run_canopus(self.gnps_mgf_file, self.canopus_dir, extra_canopus_parameters)
                except Exception as e:
                    logger.warning(
                        'Failed to run CANOPUS on mgf file with docker, error was "{}"'.format(e)
                    )
                    logger.info("Trying to run CANOPUS again using SIRIUS from path")
                    try:
                        run_canopus(self.gnps_mgf_file, self.canopus_dir, extra_canopus_parameters)
                    except Exception as e:
                        logger.warning(
                            'Again failed to run CANOPUS on mgf file using sirius from path, error was "{}"'.format(
                                e
                            )
                        )
            else:
                logger.info("Found CANOPUS dir, CANOPUS not run again!")

        # load Chem_class_predictions (canopus, molnetenhancer are loaded)
        chem_classes = ChemClassPredictions(self.canopus_dir, self.molnetenhancer_dir, self._root)
        # if no molfam classes transfer them from spectra (due to old style MN)
        if not chem_classes.canopus.molfam_classes and chem_classes.canopus.spectra_classes:
            logger.debug("Added chemical compound classes for MFs")
            chem_classes.canopus.transfer_spec_classes_to_molfams(self.molfams)
        # include them in loader
        self.chem_classes = chem_classes
        return True
