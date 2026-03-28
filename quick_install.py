'''
MIT License

Copyright (c) 2025 Ryoya Yajima

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''
import os
from pathlib import Path
import shutil

import maya.cmds as cmds
import maya.api.OpenMaya as om2

CURRENT_DIR = os.path.dirname(__file__)

def onMayaDroppedPythonFile(*args, **kwargs):
    installToCurrentVersion()


def installToCurrentVersion():
    version = cmds.about(majorVersion=True)
    if int(version) < 2025:
        raise RuntimeError('Unsupported Maya version.')

    currentDir = Path(CURRENT_DIR)

    packageDir = currentDir.joinpath('package')
    if not packageDir.exists():
        raise FileNotFoundError(f'Not found dir {packageDir.as_posix()}.')

    appDir = Path(os.environ.get("MAYA_APP_DIR"))
    moduleDirPath = appDir.joinpath('modules')
    os.makedirs(moduleDirPath, exist_ok=True)

    shutil.copytree(packageDir, moduleDirPath, dirs_exist_ok=True)

    om2.MGlobal.displayInfo(f'===== Finish install to {moduleDirPath.as_posix()}. Please restart Maya. =====')