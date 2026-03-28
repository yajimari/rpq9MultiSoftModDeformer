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

from .constants import PLUGIN_NODE_NAME


def getDeformerGeometries(deformer:str) -> list[str]:
    shapes = cmds.deformer(deformer, q=True, g=True)
    return [cmds.listRelatives(shape, p=True)[0] for shape in shapes]


def isMultiSoftModNode(mobj:om2.MObject) -> bool:
    mfn = om2.MFnDependencyNode(mobj)
    return mfn.typeName == PLUGIN_NODE_NAME


def getMultiSoftModNodes(node:str) -> list[str]:
    return cmds.ls(cmds.listHistory(node), type='rpq9MultiSoftMod')