import math
import random
import maya.cmds as cmds

PLUGIN_NAME = 'rpq9MultiSoftModDeformer'
PLUGIN_NODE_NAME = 'rpq9MultiSoftMod'

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


def createController() -> tuple[str]:
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



if __name__ == '__main__':
    controllers = createTestMultiSoftMod(5)
    setTestRandAnim(controllers)