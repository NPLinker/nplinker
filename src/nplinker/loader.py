import glob
import os
from pathlib import Path
from deprecated import deprecated
from nplinker.class_info.chem_classes import ChemClassPredictions
from nplinker.class_info.class_matches import ClassMatches
from nplinker.class_info.runcanopus import run_canopus
from nplinker.genomics import add_bgc_to_gcf
from nplinker.genomics import add_strain_to_bgc
from nplinker.genomics import generate_mappings_genome_id_bgc_id
from nplinker.genomics.antismash import AntismashBGCLoader
from nplinker.genomics.bigscape import BigscapeGCFLoader
from nplinker.genomics.mibig import MibigLoader
from nplinker.globals import GENOME_BGC_MAPPINGS_FILENAME
from nplinker.globals import GENOME_STATUS_FILENAME
from nplinker.globals import GNPS_FILE_MAPPINGS_FILENAME
from nplinker.globals import PFAM_PATH
from nplinker.globals import STRAIN_MAPPINGS_FILENAME
from nplinker.logconfig import LogConfig
from nplinker.metabolomics import add_annotation_to_spectrum
from nplinker.metabolomics import add_spectrum_to_mf
from nplinker.metabolomics import add_strains_to_spectrum
from nplinker.metabolomics.gnps import GNPSAnnotationLoader
from nplinker.metabolomics.gnps import GNPSMolecularFamilyLoader
from nplinker.metabolomics.gnps import GNPSSpectrumLoader
from nplinker.pairedomics.downloader import PODPDownloader
from nplinker.pairedomics.runbigscape import run_bigscape
from nplinker.pairedomics.strain_mappings_generator import podp_generate_strain_mappings
from nplinker.strain import Strain
from nplinker.strain_collection import StrainCollection
from nplinker.strain_loader import load_user_strains


try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

logger = LogConfig.getLogger(__name__)

NPLINKER_APP_DATA_DIR = files("nplinker").joinpath("data")


class DatasetLoader:
    ANTISMASH_DELIMITERS_DEFAULT = [".", "_", "-"]
    ANTISMASH_IGNORE_SPACES_DEFAULT = False

    TABLES_CUTOFF_DEFAULT = 2.0

    BIGSCAPE_CUTOFF_DEFAULT = 30
    EXTENDED_METADATA_TABLE_PARSING_DEFAULT = False
    USE_MIBIG_DEFAULT = True
    MIBIG_VERSION_DEFAULT = "3.1"

    RUN_BIGSCAPE_DEFAULT = True

    # https://git.wageningenur.nl/medema-group/BiG-SCAPE/-/wikis/parameters#mibig
    # BigScape mibig default version is 3.1
    EXTRA_BIGSCAPE_PARAMS_DEFAULT = "--mibig --clans-off --mix --include_singletons"

    RUN_CANOPUS_DEFAULT = False
    EXTRA_CANOPUS_PARAMS_DEFAULT = "--maxmz 600 formula zodiac structure canopus"

    # keys for overriding metabolomics data elements
    OR_NODES = "nodes_file"
    OR_EDGES = "edges_file"
    OR_EXTRA_NODES = "extra_nodes_file"
    OR_MGF = "mgf_file"
    OR_METADATA = "metadata_table_file"
    OR_QUANT = "quantification_table_file"
    OR_ANNO = "annotations_dir"
    OR_ANNO_CONFIG = "annotations_config_file"
    # and the same for genomics data
    OR_ANTISMASH = "antismash_dir"
    OR_BIGSCAPE = "bigscape_dir"
    OR_MIBIG_JSON = "mibig_json_dir"
    OR_STRAINS = "strain_mappings_file"
    # misc files
    OR_INCLUDE_STRAINS = "include_strains_file"
    # class predictions
    OR_CANOPUS = "canopus_dir"
    OR_MOLNETENHANCER = "molnetenhancer_dir"

    BIGSCAPE_PRODUCT_TYPES = [
        "NRPS",
        "Others",
        "PKSI",
        "PKS-NRP_Hybrids",
        "PKSother",
        "RiPPs",
        "Saccharides",
        "Terpene",
    ]

    def __init__(self, config_data):
        # load the config data
        self._config_dataset = config_data["dataset"]
        self._config_docker = config_data.get("docker", {})
        self._config_webapp = config_data.get("webapp", {})
        self._config_antismash = config_data.get("antismash", {})
        self._config_overrides = self._config_dataset.get("overrides", {})
        # set private attributes
        self._antismash_delimiters = self._config_antismash.get(
            "antismash_delimiters", self.ANTISMASH_DELIMITERS_DEFAULT
        )
        self._antismash_ignore_spaces = self._config_antismash.get(
            "ignore_spaces", self.ANTISMASH_IGNORE_SPACES_DEFAULT
        )
        self._bigscape_cutoff = self._config_dataset.get(
            "bigscape_cutoff", self.BIGSCAPE_CUTOFF_DEFAULT
        )
        self._extended_metadata_table_parsing = self._config_dataset.get(
            "extended_metadata_table_parsing", self.EXTENDED_METADATA_TABLE_PARSING_DEFAULT
        )
        self._use_mibig = self._config_dataset.get("use_mibig", self.USE_MIBIG_DEFAULT)
        self._mibig_version = self._config_dataset.get("mibig_version", self.MIBIG_VERSION_DEFAULT)
        self._root = Path(self._config_dataset["root"])
        self._platform_id = self._config_dataset["platform_id"]
        self._remote_loading = len(self._platform_id) > 0

        # set public attributes
        self.dataset_id = (
            os.path.split(self._root)[-1] if not self._remote_loading else self._platform_id
        )
        self.bgcs, self.gcfs, self.spectra, self.molfams = [], [], [], []
        self.mibig_bgcs = []
        self._mibig_strain_bgc_mapping = {}
        self.product_types = []
        self.strains = StrainCollection()
        self.webapp_scoring_cutoff = self._config_webapp.get(
            "tables_metcalf_threshold", self.TABLES_CUTOFF_DEFAULT
        )
        self.class_matches = None
        self.chem_classes = None

        logger.debug(
            "DatasetLoader({}, {}, {})".format(self._root, self.dataset_id, self._remote_loading)
        )

    def __repr__(self):
        return "Root={}\n   MGF={}\n   EDGES={}\n   NODES={}\n   BIGSCAPE={}\n   ANTISMASH={}\n".format(
            self._root,
            self.mgf_file,
            self.edges_file,
            self.nodes_file,
            self.bigscape_dir,
            self.antismash_dir,
        )

    def validate(self):
        """Download data and build paths for local data."""
        # if remote loading mode, need to download the data here
        # CG: for PODP workflow, strain_mappings.json is generated in the download step
        if self._remote_loading:
            self._start_downloads()

        # construct the paths and filenames required to load everything else and check
        # they all seem to exist (but don't parse anything yet)
        # TODO CG: the logics of _init_paths and _validate_paths are not clear
        # 1. after downloading (manual preparation), some files alreay exist, some not
        # 2. get the default, constructed or real path for each file/dir (need refactoring)
        #   - self._config_overrides.get()
        #   - os.path.join(self._root, 'strain_mappings.json')
        #   - find_via_glob() --> actually check if the file/dir exists
        # 3. check if (some) file/dir exists
        self._init_paths()
        self._validate_paths()

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
        if self._use_mibig:
            if not self._load_mibig():
                return False

        if not self._load_strain_mappings():
            return False

        if not self._load_metabolomics():
            return False

        if not self._load_genomics():
            return False

        if len(self.strains) == 0:
            raise Exception(f"Failed to find *ANY* strains, missing {STRAIN_MAPPINGS_FILENAME}?")

        return True

    def _start_downloads(self):
        downloader = PODPDownloader(self._platform_id)
        # TODO CG: this step generates the real path for _root. Should generate
        # it before loading process starts. Otherwise, npl.root_dir will get
        # wrong value if loading from local data or not using download.
        self._root = Path(downloader.project_results_dir)
        logger.debug("remote loading mode, configuring root=%s", self._root)
        # CG: to download both MET and GEN data
        # CG: Continue to understand how strain_mappings.json is generated
        downloader.get(
            self._config_docker.get("run_bigscape", self.RUN_BIGSCAPE_DEFAULT),
            self._config_docker.get(
                "extra_bigscape_parameters", self.EXTRA_BIGSCAPE_PARAMS_DEFAULT
            ),
            self._use_mibig,
            self._mibig_version,
        )

    def _init_paths(self):
        # 1. strain mapping are used for everything else so
        self.strain_mappings_file = self._config_overrides.get(self.OR_STRAINS) or os.path.join(
            self._root, STRAIN_MAPPINGS_FILENAME
        )

        self._init_metabolomics_paths()

        self._init_genomics_paths()

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

    def _init_metabolomics_paths(self):
        # 2. MET: <root>/clusterinfo_summary/<some UID>.tsv (or .clustersummary apparently...) / nodes_file=<override>
        self.nodes_file = self._config_overrides.get(self.OR_NODES) or find_via_glob_alts(
            [
                os.path.join(self._root, "file_mappings.tsv"),
                os.path.join(self._root, "clusterinfo*", "*.tsv"),
                os.path.join(self._root, "clusterinfo*", "*.clustersummary"),
            ],
            self.OR_NODES,
        )

        # 3. MET: <root>/networkedges_selfloop/<some UID>.selfloop (new) or .pairsinfo (old) / edges_file=<override>
        self.edges_file = self._config_overrides.get(self.OR_EDGES) or find_via_glob_alts(
            [
                os.path.join(self._root, "*.pairsinfo"),
                os.path.join(self._root, "networkedges_selfloop", "*.pairsinfo"),
                os.path.join(self._root, "networkedges_selfloop", "*.selfloop"),
            ],
            self.OR_EDGES,
        )

        # 4. MET: <root>/*.csv / extra_nodes_file=<override>
        # TODO is the glob input OK?
        # => wait for updated dataset with latest output format
        # NOTE: only optional for Crusemann or Crusemann-like dataset format!
        self.extra_nodes_file = self._config_overrides.get(self.OR_EXTRA_NODES) or find_via_glob(
            os.path.join(self._root, "quantification_table_reformatted", "*.csv"),
            self.OR_EXTRA_NODES,
            optional=True,
        )

        # 5. MET: <root>/spectra/*.mgf (or <root>/*.mgf)/ mgf_file=<override>
        self.mgf_file = self._config_overrides.get(self.OR_MGF) or find_via_glob_alts(
            [os.path.join(self._root, "*.mgf"), os.path.join(self._root, "spectra", "*.mgf")],
            self.OR_MGF,
        )

        # 6. MET: <root>/metadata_table/metadata_table-<number>.txt / metadata_table_file=<override>
        self.metadata_table_file = self._config_overrides.get(self.OR_METADATA) or find_via_glob(
            os.path.join(self._root, "metadata_table", "metadata_table*.txt"),
            self.OR_METADATA,
            optional=True,
        )

        # 7. MET: <root>/quantification_table/quantification_table-<number>.csv / quantification_table_file=<override>
        self.quantification_table_file = self._config_overrides.get(self.OR_QUANT) or find_via_glob(
            os.path.join(self._root, "quantification_table", "quantification_table*.csv"),
            self.OR_QUANT,
            optional=True,
        )

        # 8. MET: <root>/DB_result/*.tsv (new) or <root>/result_specnets_DB/*.tsv (old) / annotations_dir=<override>
        if Path.is_file(Path(self._root) / "annotations.tsv"):
            self.annotations_dir = str(self._root)
            self.annotations_config_file = os.path.join(self._root, "annotations.tsv")
        else:
            self.annotations_dir = self._config_overrides.get(
                self.OR_ANNO
            ) or find_via_glob_alts_dir(
                [
                    os.path.join(self._root, "DB_result"),
                    os.path.join(self._root, "result_specnets_DB"),
                ],
                self.OR_ANNO,
                optional=False,
            )
            if self.annotations_dir is not None:
                self.annotations_config_file = self._config_overrides.get(
                    self.OR_ANNO_CONFIG
                ) or os.path.join(self.annotations_dir, "annotations.tsv")

    def _init_genomics_paths(self):
        # 9. GEN: <root>/antismash / antismash_dir=<override>
        self.antismash_dir = self._config_overrides.get(self.OR_ANTISMASH) or os.path.join(
            self._root, "antismash"
        )
        self.antismash_cache = {}

        # 10. GEN: <root>/bigscape / bigscape_dir=<override>
        self.bigscape_dir = self._config_overrides.get(self.OR_BIGSCAPE) or os.path.join(
            self._root, "bigscape"
        )
        # what we really want here is the subdirectory containing the network/annotation files,
        # but in case this is the path to the top level bigscape output, try to figure it out automatically
        if not os.path.exists(os.path.join(self.bigscape_dir, "mix")):
            new_bigscape_dir = find_bigscape_dir(self.bigscape_dir)
            if new_bigscape_dir is not None:
                logger.info(
                    "Updating bigscape_dir to discovered location {}".format(new_bigscape_dir)
                )
                self.bigscape_dir = new_bigscape_dir

        # 11. GEN: <root>/mibig_json / mibig_json_dir=<override>
        self.mibig_json_dir = self._config_overrides.get(self.OR_MIBIG_JSON) or os.path.join(
            self._root, "mibig_json"
        )

    def _validate_paths(self):
        """Validates that the required files and directories exist before loading starts."""
        required_paths = [self.nodes_file, self.edges_file, self.mgf_file, self.antismash_dir]
        optional_paths = [self.annotations_dir]

        for f in required_paths:
            if not os.path.exists(str(f)):
                raise FileNotFoundError(f'File/directory "{f}" does not exist.')

        for f in optional_paths:
            if not os.path.exists(str(f)):
                logger.warning('Optional file/directory "%s" does not exist', f)

    def _load_mibig(self):
        mibig_bgc_loader = MibigLoader(self.mibig_json_dir)
        self.mibig_bgcs = mibig_bgc_loader.get_bgcs()
        self._mibig_strain_bgc_mapping = mibig_bgc_loader.get_strain_bgc_mapping()
        return True

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

        # 3. load MiBIG strain mappings
        if self._mibig_strain_bgc_mapping:
            for k, v in self._mibig_strain_bgc_mapping.items():
                strain = Strain(k)
                strain.add_alias(v)
                self.strains.add(strain)
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
        raw_spectra = GNPSSpectrumLoader(self.mgf_file).spectra
        # Step 2: load all GNPS annotations
        raw_annotations = GNPSAnnotationLoader(self.annotations_config_file).annotations
        # Step 3: load all MolecularFamily objects
        raw_molfams = GNPSMolecularFamilyLoader(self.edges_file).get_mfs(keep_singleton=False)

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
        not added to `self.bgcs`.

        The attribute of `self.gcfs` is set to the loaded GCF objects that have the Strain objects
        added (i.e. `GCF._strains` updated). This means only BGC objects with updated Strain objects
        (i.e. `self.bgcs`) can be added to GCF objects.
        """
        logger.debug("\nLoading genomics data starts...")

        # Step 1: load all BGC objects
        logger.debug("Parsing AntiSMASH directory...")
        antismash_bgcs = AntismashBGCLoader(self.antismash_dir).get_bgcs()
        raw_bgcs = antismash_bgcs + self.mibig_bgcs

        # Step 2: load all GCF objects
        # TODO: create a config for "bigscape_cluster_file" and discard "bigscape_dir" and "bigscape_cutoff"?
        bigscape_cluster_file = (
            Path(self.bigscape_dir) / "mix" / f"mix_clustering_c0.{self._bigscape_cutoff:02d}.tsv"
        )
        raw_gcfs = BigscapeGCFLoader(bigscape_cluster_file).get_gcfs()

        # Step 3: add Strain object to BGC
        bgc_with_strain, _ = add_strain_to_bgc(self.strains, raw_bgcs)

        # Step 4: add BGC objects to GCF
        gcf_with_bgc, _, _ = add_bgc_to_gcf(bgc_with_strain, raw_gcfs)

        # Step 5: set attributes of self.bgcs and self.gcfs with valid objects
        self.bgcs = bgc_with_strain
        self.gcfs = gcf_with_bgc

        logger.debug("Loading genomics data completed\n")
        return True

    # TODO CG: run bigscape before loading and after downloading
    def _run_bigscape(self):
        # Check for spaces in the folder names under <dataset>/antismash and
        # rename them by replacing spaces with underscores
        if not self._antismash_ignore_spaces:
            logger.debug("Checking for spaces in antiSMASH folder names...")
            for root, dirs, _ in os.walk(self.antismash_dir):
                for d in dirs:
                    if d.find(" ") != -1:
                        os.rename(os.path.join(root, d), os.path.join(root, d.replace(" ", "_")))
                        logger.warn(
                            'Renaming antiSMASH folder "{}" to "{}" to remove spaces! (suppress with ignore_spaces = true in config file)'.format(
                                d, d.replace(" ", "_")
                            )
                        )

        if not os.path.exists(self.bigscape_dir):
            should_run_bigscape = self._config_docker.get("run_bigscape", self.RUN_BIGSCAPE_DEFAULT)
            extra_bigscape_parameters = self._config_docker.get(
                "extra_bigscape_parameters", self.EXTRA_BIGSCAPE_PARAMS_DEFAULT
            )
            if should_run_bigscape:
                logger.info(
                    'Running BiG-SCAPE! extra_bigscape_parameters="{}"'.format(
                        extra_bigscape_parameters
                    )
                )
                try:
                    run_bigscape(
                        "bigscape.py",
                        os.path.join(self._root, "antismash"),
                        os.path.join(self._root, "bigscape"),
                        PFAM_PATH,
                        extra_params=extra_bigscape_parameters,
                    )
                except Exception as e:
                    logger.warning(
                        'Failed to run BiG-SCAPE on antismash data, error was "{}"'.format(e)
                    )

                self.bigscape_dir = find_bigscape_dir(self.bigscape_dir)

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
                    run_canopus(self.mgf_file, self.canopus_dir, extra_canopus_parameters)
                except Exception as e:
                    logger.warning(
                        'Failed to run CANOPUS on mgf file with docker, error was "{}"'.format(e)
                    )
                    logger.info("Trying to run CANOPUS again using SIRIUS from path")
                    try:
                        run_canopus(self.mgf_file, self.canopus_dir, extra_canopus_parameters)
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


def find_via_glob(path, file_type, optional=False):
    try:
        filename = glob.glob(path)[0]
        return filename
    except (OSError, IndexError):
        if not optional:
            # "from None" suppresses the traceback for the original exception, which isn't really needed
            raise Exception(
                'ERROR: unable to find {} in path "{}"'.format(file_type, path)
            ) from None

        logger.warn('WARNING: unable to find {} in path "{}"'.format(file_type, path))
        return None


def find_via_glob_alts_dir(paths, file_type, optional=False):
    path = None
    for p in paths:
        if os.path.exists(p):
            path = p
            break

    if path is None and not optional:
        raise Exception(
            "ERROR: unable to find {} in {} paths: ({})".format(file_type, len(paths), paths)
        )
    elif path is None:
        logger.warning(
            "WARNING: unable to find {} in {} paths: ({})".format(file_type, len(paths), paths)
        )

    return path


def find_via_glob_alts(paths, file_type, optional=False):
    filename = None
    for path in paths:
        try:
            filename = glob.glob(path)[0]
            break
        except (OSError, IndexError):
            continue

    if filename is None and not optional:
        raise Exception(
            "ERROR: unable to find {} in {} paths: ({})".format(file_type, len(paths), paths)
        )
    elif filename is None:
        logger.warning(
            "WARNING: unable to find {} in {} paths: ({})".format(file_type, len(paths), paths)
        )

    return filename


def find_bigscape_dir(broot):
    logger.info(f"Trying to discover correct bigscape directory under {broot}")
    for root, _, files in os.walk(broot):
        if "Network_Annotations_Full.tsv" in files:
            logger.info(f"Found network files directory: {root}")
            return root

    return None
