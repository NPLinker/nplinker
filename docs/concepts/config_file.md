## Configuration Template

```toml
--8<-- "config/nplinker.toml"
```

## Example Configuration

For a full example of a configuration file, see [here](../quickstart.md#3-prepare-config-file).

## Config loader

You can load the configuration file using the [load_config](../api/nplinker.md#nplinker.config.load_config) function.

```python
from nplinker.config import load_config
config = load_config('path/to/nplinker.toml')
```

When you use NPLinker as an application, you can get access to the configuration object directly:

```python
from nplinker import NPLinker
npl = NPLinker('path/to/nplinker.toml')
print(npl.config)
```