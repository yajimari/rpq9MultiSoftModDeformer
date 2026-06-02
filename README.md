# rpq9MultiSoftModDeformer

![License: MIT](https://img.shields.io/badge/license-MIT-brightgreen.svg)

**rpq9MultiSoftModDeformer** is a custom deformer node for Autodesk Maya that computes multiple softMod deformers within a single node.  


## 📝Features

- Each softMod is composited like a blendShape, so vertices where multiple softMods interfere are deformed intuitively.

- When running multiple softMods, it will be faster than the standard softMod.

- A dedicated handle node for the deformer is not required.

## ⚠️ Limitations

- Supported in Maya 2025 and later versions. Not available in earlier versions.

- Falloff weight using curves is not available.

- Rotation interpolation always takes shortest path in quaternion space.

## 🚀Quick Start

1. Download the files in this repository from the Releases page.
2. Copy `rpq9MultiSoftModDeformer` and `rpq9MultiSoftModDeformer.mod` in the package below to the directory set in `MAYA_MODULE_PATH`.(You can check the destination directory by running the following code.)
```python
import os
from pathlib import Path

appDir = Path(os.environ.get('MAYA_APP_DIR'))
moduleDir = appDir.joinpath('modules')
print('Installation directory: ', moduleDir.as_posix())
```
3. start Maya and run the code below. 
```python
import rpq9MultiSoftMod.sample

rpq9MultiSoftMod.sample.createAnimetedMultiSoftModScene()
```

## 📘Documentation
Detailed documentation is provided in a separate file.  
Please refer to `./reference.md`.

## 🧩Plugin Directory Structure

The correct file placement when using the plugin is as follows (Windows).
```
 %MAYA_MODULE_PATH%
 └──rpq9MultiSoftModDeformer
    └──scripts
    └──plug-ins
```

## ⚙️For developers

- When building the plugin yourself, you can build it in the same way as the Maya devkit examples.

- If you use the plugin with modifications, you need to change the values of the variables in <u>nodeRegisterData.h</u>, as they could conflict with the base plugin.