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
import maya.cmds as cmds
import maya.api.OpenMaya as om2

from .constants import PLUGIN_NAME, PLUGIN_NODE_NAME
from .controller import createController

def loadPlugin() -> None:
    if not cmds.pluginInfo(PLUGIN_NAME, q=True, loaded=True):
        cmds.loadPlugin(PLUGIN_NAME, qt=True)

def getDeformerGeometries(deformer:str) -> list[str]:
    shapes = cmds.deformer(deformer, q=True, g=True)
    return [cmds.listRelatives(shape, p=True)[0] for shape in shapes]


def isMultiSoftModNode(node:str) -> bool:
    return cmds.nodeType(node) == PLUGIN_NODE_NAME


def isMultiSoftModNodeMObj(mobj:om2.MObject) -> bool:
    mfn = om2.MFnDependencyNode(mobj)
    return mfn.typeName == PLUGIN_NODE_NAME


def getMultiSoftModNodes(node:str) -> list[str]:
    return cmds.ls(cmds.listHistory(node), type=PLUGIN_NODE_NAME)


def getMPlug(attr:str) -> om2.MPlug:
    selList= om2.MSelectionList()
    selList.add(attr)
    return selList.getPlug(0)


def createRpq9MultiSoftMod(geos:list[str]|str, controllerNum:int=0, **kwargs):
    kwargs.pop('type', None)
    kwargs.pop('typ', None)
    kwargs['type'] = PLUGIN_NODE_NAME
    deformer = cmds.deformer(geos, **kwargs)[0]
    controllers = []
    for i in range(controllerNum):
        centerController, modifyController = createController()
        cmds.connectAttr(f'{centerController}.envelope', f'{deformer}.inputData[{i}].localEnvelope', f=True)
        cmds.connectAttr(f'{centerController}.wm[0]', f'{deformer}.inputData[{i}].centerMatrix', f=True)
        cmds.connectAttr(f'{modifyController}.wm[0]', f'{deformer}.inputData[{i}].modifyMatrix', f=True)
        cmds.connectAttr(f'{centerController}.falloffRadius', f'{deformer}.inputData[{i}].falloffRadius', f=True)
        controllers.append(centerController)
        controllers.append(modifyController)
    res = [deformer] + controllers
    return res