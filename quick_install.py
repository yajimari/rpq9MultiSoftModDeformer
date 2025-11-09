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
    platform = cmds.about(os=True)
    platformDir = ''
    if platform == 'win64':
        platformDir = 'windows'
    elif platform == 'mac64':
        platformDir = 'osx'
    elif platform == 'linux64':
        platformDir = 'linux'
    else:
        raise RuntimeError('Unknown platform.')

    pluginDir = currentDir.joinpath('plug-ins', platformDir, version)
    if not pluginDir.exists():
        raise FileNotFoundError(f'Not builded Maya{version} plugin.')

    appDir = Path(os.environ.get("MAYA_APP_DIR"))
    currentVerPluginPath = appDir.joinpath(cmds.about(majorVersion=True), 'plug-ins')
    os.makedirs(currentVerPluginPath, exist_ok=True)

    for file in [path for path in pluginDir.glob("*") if path.is_file()]:
        shutil.copy2(file, currentVerPluginPath)

    om2.MGlobal.displayInfo(f'===== Finish install to {currentVerPluginPath.as_posix()} =====')