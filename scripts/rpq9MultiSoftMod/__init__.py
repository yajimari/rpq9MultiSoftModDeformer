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
import importlib.util
from .constants import PLUGIN_NAME, PLUGIN_NODE_NAME
from .model import Matrix4x4, DeformerWeights, DeformerWeightList, Rpq9MultiSoftModData


if importlib.util.find_spec('maya') is not None:
    import maya.cmds as cmds
    from . import core, deformerWeight, deformerWeightUI, utils

    utils.loadPlugin()

    if not cmds.runTimeCommand('CreateRpq9MultiSoftMod', exists=True):
        cmds.runTimeCommand(
            'CreateRpq9MultiSoftMod',
            annotation='Execute CreateRpq9MultiSoftMod',
            command='import rpq9MultiSoftMod.command;cmd = rpq9MultiSoftMod.command.Command();cmd.executeWithPreferences()',
            default=True)

    if not cmds.runTimeCommand('CreateRpq9MultiSoftModOptions', exists=True):
        cmds.runTimeCommand(
            'CreateRpq9MultiSoftModOptions',
            annotation='Show CreateRpq9MultiSoftModOptions',
            command='import rpq9MultiSoftMod.command;cmd = rpq9MultiSoftMod.command.Command();cmd.createDialog(optionVarOverrideDict=None, saveOptionVars=True)',
            default=True)