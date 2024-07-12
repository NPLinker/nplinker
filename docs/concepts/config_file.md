## Configuration Template

```toml
--8<-- "src/nplinker/data/nplinker.toml"
```


## Default Configurations
The default configurations are automatically used by NPLinker if you don't set them in your config file.

```toml
--8<-- "src/nplinker/nplinker_default.toml"
```

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