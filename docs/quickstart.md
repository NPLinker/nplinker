Now, let's start the jorney with NPLinker.

## 1. Decide which mode to use
NPLinker allows you to run in two modes:

=== "`local` mode"
    The `local` mode assumes that the data required by NPLinker is available on your local machine.

    The required input data includes:

    - GNPS molecular networking data from one of the following GNPS workflows
        - `METABOLOMICS-SNETS`,
        - `METABOLOMICS-SNETS-V2`
        - `FEATURE-BASED-MOLECULAR-NETWORKING`
    - AntiSMASH BGC data
    - BigScape data (optional)


=== "`podp` mode"
    The `podp` mode assumes that you use an identifier of
    [Paired Omics Data Platform](https://pairedomicsdata.bioinformatics.nl/) (PODP) as the input for
    NPLinker. Then NPLinker will download and prepare all data necessary based on the PODP id which
    refers to the metadata of the dataset.


So, which mode will you use? The answer is important for the next steps.


## 2. Create a working directory
The working directory is used to store all input and output data for NPLinker. You can name this
directory as you like, for example `nplinker_quickstart`:

```bash title="Create a working directory"
mkdir nplinker_quickstart
```

!!! warning "Important"
    Before going to the next step, make sure you get familiar with how NPLinker organizes data in the
    working directory, see [Working Directory Structure](./concepts/working_dir_structure.md) page.


## 3. Prepare input data (`local` mode only)
If you choose to use the `local` mode, meaning you have input data of NPLinker stored on your local
machine, you need to move the input data to the working directory created in the previous step.
Skip this step if you choose to use the `podp` mode.

### GNPS data
NPLinker accepts data from the output of the following GNPS workflows:

- `METABOLOMICS-SNETS`
- `METABOLOMICS-SNETS-V2`
- `FEATURE-BASED-MOLECULAR-NETWORKING`.

NPLinker provides the tools [`GNPSDownloader`][nplinker.metabolomics.gnps.GNPSDownloader] and
[`GNPSExtractor`][nplinker.metabolomics.gnps.GNPSExtractor] to download and extract the GNPS data
with ease. What you need to give is a valid GNPS task ID, referring to a task of the GNPS workflows
supported by NPLinker.

??? example "GNPS task id and workflow"
    Given an example of GNPS task at https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=c22f44b14a3d450eb836d607cb9521bb,
    the task id is the last part of this url, i.e. `c22f44b14a3d450eb836d607cb9521bb`. Open this link,
    you can find the worklow info at the row "Workflow" of the table "Job Status", for this case,
    it is `METABOLOMICS-SNETS`.

```python title="Download & Extract GNPS data"
from nplinker.metabolomics.gnps import GNPSDownloader, GNPSExtractor

# Go to the working directory
cd nplinker_quickstart

# Download GNPS data & get the path to the downloaded archive
downloader = GNPSDownloader("gnps_task_id", "downloads") # (1)!
downloaded_archive = downloader.download().get_download_file()

# Extract GNPS data to `gnps` directory
extractor = GNPSExtractor(downloaded_archive, "gnps") # (2)!
```

1. If you already have the downloaded archive of GNPS data, you can skip the download steps.
2. Replace `downloaded_archive` with the actuall path to your GNPS data archive if you skipped the download steps.

The required data for NPLinker will be extracted to the `gnps` subdirectory of the working directory.

!!! info
    Not all GNPS data are required by NPLinker, and only the necessary data will be extracted.
    During the extraction, these data will be renamed to the standard names used by NPLinker.
    See the page [GNPS Data](./concepts/gnps_data.md) for more information.

??? tip "Prepare GNPS data manually"
    If you have GNPS data but it is not the archive format as downloaded from GNPS, it's recommended
    to re-download the data from GNPS.

    If (re-)downloading is not possible, you could manually prepare data for the `gnps` directory.
    In this case, you must make sure that the data is organized as expected by NPLinker.
    See the page [GNPS Data](./concepts/gnps_data.md) for examples of how to prepare the data.

### AntiSMASH data
NPLinker requires AntiSMASH BGC data as input, which are organized in the `antismash` subdirectory of 
the working directory.

For each output of AntiSMASH run, the BGC data must be stored in a subdirectory named after the NCBI
accession number (e.g. `GCF_000514975.1`). And only the `*.region*.gbk` files are required by NPLinker.

When manually preparing AntiSMASH data for NPLinker, you must make sure that the data is organized as
expected by NPLinker. See the page [Working Directory Structure](./concepts/working_dir_structure.md)
for more information.

### BigScape data (optional)
It is optional to provide the output of BigScape to NPLinker. If the output of BigScape is not provided,
NPLinker will run BigScape automatically to generate the data using the AntiSMASH BGC data.

If you have the output of BigScape, you can put its `mix_clustering_c{cutoff}.tsv` file in the
`bigscape` subdirectory of the NPLinker working directory, where `{cutoff}` is the cutoff value used
in the BigScape run.

### Strain mappings file



## 4. Prepare confg file


=== "`local` mode"

    You can use local mode to run the application.

    ```bash
    local
    ```


=== "`podp` mode"

    you can use podp mode.


## 5. Run NPlinker
