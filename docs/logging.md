NPLinker uses the standard library [logging](https://docs.python.org/3/library/logging.html#module-logging) 
module for managing log messages and the python library [rich](https://rich.readthedocs.io/en/latest/logging.html) 
to colorize the log messages. Depending on how you use NPLinker, you can set up logging in different ways.

## NPLinker as an application
If you're using NPLinker as an application, you're running the whole workflow of NPLinker as 
described in the [Quickstart](./quickstart.md). In this case, you can set up logging in the nplinker 
configuration file `nplinker.toml`. 

## NPLinker as a library
If you're using NPLinker as a library, you're using only some functions and classes of NPLinker in 
your script. By default, NPLinker will not log any messages. However, you can set up logging in your
script to log messages. 

```python title="Set up logging in 'your_script.py'"
# Set up logging configuration first
from nplinker import setup_logging

setup_logging(level="DEBUG", file="nplinker.log", use_console=True) # (1)!

# Your business code here
# e.g. download and extract nplinker example data
from nplinker.utils import download_and_extract_archive

download_and_extract_archive(
    url="https://zenodo.org/records/10822604/files/nplinker_local_mode_example.zip",
    download_root=".",
)
```

1. The `setup_logging` function sets up the logging configuration. The `level` argument sets the 
   logging level. The `file` argument sets the log file. The `use_console` argument sets whether to 
   log messages to the console.


The log messages will be written to the log file `nplinker.log` and displayed in the console with a 
format like this: `[Date Time] Level Log-message Module:Line`.

```shell title="Run your script in a terminal"
# Run your script
$ python your_script.py
Downloading nplinker_local_mode_example.zip ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100.0% • 195.3/195.3 MB • 2.6 MB/s • 0:00:00 • 0:01:02 # (1)!
[2024-05-10 15:14:48] INFO     Extracting nplinker_local_mode_example.zip to .                      utils.py:401

# Check the log file
$ cat nplinker.log
[2024-05-10 15:14:48] INFO     Extracting nplinker_local_mode_example.zip to .                      utils.py:401
```

1. This is a progress bar but not a log message.