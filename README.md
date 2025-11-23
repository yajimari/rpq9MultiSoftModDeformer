# rpq9MultiSoftModDeformer

![License: MIT](https://img.shields.io/badge/license-MIT-brightgreen.svg)

**rpq9MultiSoftModDeformer** is a custom deformer node for Autodesk Maya that computes multiple softMod deformers within a single node.  


## 📝Features

- Each softMod is composited like a blendShape, so vertices where multiple softMods interfere are deformed intuitively.

- When running multiple softMods, it will be faster than the standard softMod.(In the scene created with test_rpq9MultiSoftMod.py, it ran faster than the standard version when there were around five to six or more softMod calculations.)

- A dedicated handle node for the deformer is not required.

## ⚠️ Limitations

- Supported in Maya 2025 and later versions. Not available in earlier versions.

- Falloff weight using curves is not available.

- Rotation interpolation always takes shortest path in quaternion space.

## 🚀Quick Start

1. Download the files in this repository from the Releases page.
2. Extract the downloaded file and drag and drop **quick_install.py** into the Maya you want to install.
3. Copy the contents of test/test.py into the python tab of the script editor and run it.
4. If you get a plugin not found error, Please restart Maya or manually load the plugin from the Plugin Manage.

## 🧩Plugin Directory Structure

The correct file placement when using the plugin is as follows (Windows).
```
 %MAYA_PLUG_IN_PATH%
 └──rpq9MultiSoftModDeformer.mll
```

## ⚙️For developers

- When building the plugin yourself, you can build it in the same way as the Maya devkit examples.

- If you use the plugin with modifications, you need to change the values of the variables in <u>nodeRegisterData.h</u>, as they could conflict with the base plugin.