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
import maya.cmds as cmds
import maya.api.OpenMaya as om2


def createSphereCurve() -> str:
    numSides = 8
    allPoints = []
    xyPoints = []
    xzPoints = []
    yzPoints = []

    for i in range(numSides):
        ang = math.radians(i * (360.0 / numSides))
        sin_value = math.sin(ang)
        cos_value = math.cos(ang)
        xzPoints.append((cos_value, 0.0, sin_value))
        xyPoints.append((cos_value, sin_value, 0.0))
        yzPoints.append((0.0, cos_value, sin_value))

    allPoints = xzPoints + xyPoints + xyPoints[:2] + yzPoints + [yzPoints[0], xyPoints[1], xyPoints[0]]
    crv = cmds.curve(p=allPoints, degree=1)
    cmds.closeCurve(crv, replaceOriginal=True)
    return crv


def createCrossCurve() -> str:
    allPoints = []
    for i in range(3):
        for v in (1.0, -1.0):
            pos = [0.0, 0.0, 0.0]
            pos[i] = v
            allPoints.append(pos)
        allPoints.append([0.0, 0.0, 0.0])
    crv = cmds.curve(p=allPoints, degree=1)
    return crv


def createController() -> tuple[str, str]:
    centerController = createSphereCurve()
    modifyController = createCrossCurve()
    cluster, handle = cmds.cluster(centerController, modifyController, rel=True)
    cmds.parent(modifyController, handle, centerController)
    cmds.addAttr(centerController, ln='envelope', min=0.0, max=1.0, dv=1.0, k=True)
    cmds.addAttr(modifyController, ln='envelope', pxy=f'{centerController}.envelope')
    cmds.addAttr(centerController, ln='falloffRadius', min=0.0, dv=1.0, k=True)
    cmds.addAttr(modifyController, ln='falloffRadius', pxy=f'{centerController}.falloffRadius')

    for axis in 'xyz':
        cmds.connectAttr(f'{centerController}.falloffRadius', f'{handle}.s{axis}', f=True)
        cmds.setAttr(f'{centerController}.s{axis}', l=True, k=False, cb=False)

    for attr in 'trs':
        for axis in 'xyz':
            cmds.setAttr(f'{handle}.{attr}{axis}', l=True, k=False, cb=False)

    cmds.setAttr(f'{centerController}.v', l=True, k=False, cb=False)
    cmds.setAttr(f'{modifyController}.v', l=True, k=False, cb=False)
    cmds.setAttr(f'{handle}.v', 0, l=True, k=False, cb=False)

    return centerController, modifyController