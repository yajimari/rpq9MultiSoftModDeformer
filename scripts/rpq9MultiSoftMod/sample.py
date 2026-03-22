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
import math
import random
import maya.cmds as cmds

from .constants import PLUGIN_NAME, PLUGIN_NODE_NAME
from .controller import createController


def createTestMultiSoftMod(num:int) -> list[list[str]]:
    if not cmds.pluginInfo(PLUGIN_NAME, q=True, loaded=True):
        cmds.loadPlugin(PLUGIN_NAME, qt=True)
    plane = cmds.polyPlane()
    deformer = cmds.deformer(plane, type=PLUGIN_NODE_NAME)[0]
    controllers = []
    for i in range(num):
        centerController, modifyController = createController()
        cmds.connectAttr(f'{centerController}.envelope', f'{deformer}.inputData[{i}].localEnvelope', f=True)
        cmds.connectAttr(f'{centerController}.wm[0]', f'{deformer}.inputData[{i}].centerMatrix', f=True)
        cmds.connectAttr(f'{modifyController}.wm[0]', f'{deformer}.inputData[{i}].modifyMatrix', f=True)
        cmds.connectAttr(f'{centerController}.falloffRadius', f'{deformer}.inputData[{i}].falloffRadius', f=True)
        controllers.append([centerController, modifyController])
    return controllers


def createTestMultiSoftModSeries(num:int) -> list[list[str]]:
    if not cmds.pluginInfo(PLUGIN_NAME, q=True, loaded=True):
        cmds.loadPlugin(PLUGIN_NAME, qt=True)
    plane = cmds.polyPlane()
    controllers = []
    for i in range(num):
        deformer = cmds.deformer(plane, type='rpq9MultiSoftMod')[0]
        centerController, modifyController = createController()
        cmds.connectAttr(f'{centerController}.envelope', f'{deformer}.inputData[0].localEnvelope', f=True)
        cmds.connectAttr(f'{centerController}.wm[0]', f'{deformer}.inputData[0].centerMatrix', f=True)
        cmds.connectAttr(f'{modifyController}.wm[0]', f'{deformer}.inputData[0].modifyMatrix', f=True)
        cmds.connectAttr(f'{centerController}.falloffRadius', f'{deformer}.inputData[0].falloffRadius', f=True)
        controllers.append([centerController, modifyController])
    return controllers


def createTestSoftMod(num:int) -> list[list[str]]:
    plane = cmds.polyPlane()
    controllers = []
    for i in range(num):
        deformer, handle = cmds.softMod(plane)
        centerController, modifyController = createController()
        cmds.parent(handle, centerController)
        cmds.connectAttr(f'{centerController}.translate', f'{deformer}.falloffCenter', f=True)
        cmds.connectAttr(f'{centerController}.im', f'{deformer}.bindPreMatrix', f=True)
        cmds.connectAttr(f'{modifyController}.translate', f'{handle}.translate', f=True)
        cmds.connectAttr(f'{modifyController}.rotate', f'{handle}.rotate', f=True)
        cmds.connectAttr(f'{modifyController}.scale', f'{handle}.scale', f=True)
        cmds.connectAttr(f'{modifyController}.shear', f'{handle}.shear', f=True)
        cmds.connectAttr(f'{centerController}.falloffRadius', f'{deformer}.falloffRadius', f=True)
        for attr in 'trs':
            for axis in 'xyz':
                cmds.setAttr(f'{handle}.{attr}{axis}', l=True, k=False, cb=False)
        cmds.setAttr(f'{handle}.v', 0, l=True, k=False, cb=False)
        controllers.append([centerController, modifyController])
    return controllers


def setTestRandAnim(controllers:list[list[str]], animAttr:str='ty') -> None:
    random.seed(0)

    for centerController, modifyController in controllers:
        cmds.setAttr(f'{centerController}.translate', random.random()-0.5, random.random()-0.5, random.random()-0.5, type='double3')
        cmds.setKeyframe(modifyController, attribute=animAttr, time=1, value=random.random())
        cmds.setKeyframe(modifyController, attribute=animAttr, time=100, value=random.random())


def createAnimetedMultiSoftModScene(softModNum:int=5):
    cmds.file(new=True, f=True)
    controllers = createTestMultiSoftMod(softModNum)
    setTestRandAnim(controllers)


def createAnimetedMultiSoftModSeriesScene(softModNum:int=5):
    cmds.file(new=True, f=True)
    controllers = createTestMultiSoftModSeries(softModNum)
    setTestRandAnim(controllers)


def createAnimetedSoftModScene(softModNum:int=5):
    cmds.file(new=True, f=True)
    controllers = createTestSoftMod(softModNum)
    setTestRandAnim(controllers)