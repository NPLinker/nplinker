NPLinker can run BigScape automatically if the `bigscape` directory does not exist in the working directory.

To run BigScape, NPLinker requires the following BigScape parameters:

- `--mix`
- `--include_singletons`
- `--cutoffs`

And the following parameters are not allowed:

- `--inputdir`
- `--outputdir`
- `--pfam_dir`

If BigScape parameter `--mibig` is set, make sure setting the  `mibig.to_use` to true in your config file `nplinker.toml` and `mibig.version` to the version of mibig used by bigscape.


See the [default configurations](./config_file.md#default-configurations) for the default 
parameters of BigScape.