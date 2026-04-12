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
from .utils import isMultiSoftModNode
from .jsonIO import readJson, writeJson
from .model import Rpq9MultiSoftModData, Matrix4x4
from .deformerWeight import getVertexWeightAttributeData, setVertexWeightData


def getRpq9MultiSoftModData(node:str) -> Rpq9MultiSoftModData:
    if not isMultiSoftModNode(node):
        raise TypeError(f'"{node}" is not "{PLUGIN_NODE_NAME}" type.')

    kwargs = { 
            'name': node,
            'geometries': cmds.deformer(node, q=True, geometry=True),
            'geometryIndices': cmds.deformer(node, q=True, geometryIndices=True),
            'envelope': cmds.getAttr(f'{node}.envelope'),
            'falloffMode': cmds.getAttr(f'{node}.falloffMode'),
            'localEnvelopes': {},
            'centerMatrices': {},
            'modifyMatrices': {},
            'falloffRadii': {},
            'weightList': None,
            'localWeightList': {}
            }

    kwargs['weightList'] = getVertexWeightAttributeData(f'{node}.weightList')

    kwargs['localWeightList'] = {}
    inputDataMultiIndices = cmds.getAttr(f'{node}.inputData', mi=True) or []
    for inputDataMultiIndex in inputDataMultiIndices:
        kwargs['localEnvelopes'][str(inputDataMultiIndex)] = cmds.getAttr(f'{node}.inputData[{inputDataMultiIndex}].localEnvelope')
        kwargs['centerMatrices'][str(inputDataMultiIndex)] = Matrix4x4(*cmds.getAttr(f'{node}.inputData[{inputDataMultiIndex}].centerMatrix'))
        kwargs['modifyMatrices'][str(inputDataMultiIndex)] = Matrix4x4(*cmds.getAttr(f'{node}.inputData[{inputDataMultiIndex}].modifyMatrix'))
        kwargs['falloffRadii'][str(inputDataMultiIndex)] = cmds.getAttr(f'{node}.inputData[{inputDataMultiIndex}].falloffRadius')
        kwargs['localWeightList'][str(inputDataMultiIndex)] = getVertexWeightAttributeData(f'{node}.inputData[{inputDataMultiIndex}].localWeightList')

    return Rpq9MultiSoftModData(**kwargs)


def saveRpq9MultiSoftModData(node:str, filePath:str) -> None:
    data = getRpq9MultiSoftModData(node)
    writeJson(filePath, data.toDict())
    om2.MGlobal.displayInfo(f'Saved "{node}" rpq9MultiSoftMod data: {filePath}')


def applyRpq9MultiSoftModData(softModData:Rpq9MultiSoftModData):
    softModData.reindexingGeometry()
    deformer = cmds.deformer(softModData.geometries, type=PLUGIN_NODE_NAME, n=softModData.name)[0]
    cmds.setAttr(f'{deformer}.envelope', softModData.envelope)
    cmds.setAttr(f'{deformer}.falloffMode', softModData.falloffMode)

    for matrixIndex in softModData.centerMatrices.keys():
        cmds.setAttr(f'{deformer}.inputData[{matrixIndex}].localEnvelope', softModData.localEnvelopes[matrixIndex])
        cmds.setAttr(f'{deformer}.inputData[{matrixIndex}].centerMatrix', *softModData.centerMatrices[matrixIndex].toList(), type='matrix')
        cmds.setAttr(f'{deformer}.inputData[{matrixIndex}].modifyMatrix', *softModData.modifyMatrices[matrixIndex].toList(), type='matrix')
        cmds.setAttr(f'{deformer}.inputData[{matrixIndex}].falloffRadius', softModData.falloffRadii[matrixIndex])

    setVertexWeightData(deformer, softModData)

    return deformer


def applyRpq9MultiSoftModDataFile(filePath:str):
    data = readJson(filePath)
    softModData = Rpq9MultiSoftModData.fromJsonDict(data)
    return applyRpq9MultiSoftModData(softModData)


def loadVertexWeightData(filePath:str) -> None:
    data = readJson(filePath)
    softModData = Rpq9MultiSoftModData.fromJsonDict(data)
    deformerName = softModData.name
    if not cmds.objExists(deformerName):
        raise RuntimeError(f'"{deformerName}" is not found.')
    setVertexWeightData(deformerName, softModData)


def loadVertexWeightDataToNode(node:str, filePath:str):
    data = readJson(filePath)
    softModData = Rpq9MultiSoftModData.fromJsonDict(data)
    if not cmds.objExists(node):
        raise RuntimeError(f'"{node}" is not found.')
    setVertexWeightData(node, softModData)


def isValidRpq9MultiSoftModAttribute(node:str, verbose:bool=False) -> bool:
    res = True
    softModData = getRpq9MultiSoftModData(node)

    if not softModData.isValid(verbose):
        if verbose:
            res = False
        else:
            return False

    vertexCountData = {str(geoIndex): cmds.polyEvaluate(geo, v=True) for geo, geoIndex in zip(softModData.geometries, softModData.geometryIndices)}

    for geoIndex, weights in softModData.weightList.items():
        for vertexIndex, weight in weights.items():
            if int(vertexIndex) >= vertexCountData[geoIndex]:
                if verbose:
                    om2.MGlobal.displayInfo(f'Using an index larger than the number of vertices in weightList.(geoIndex: {geoIndex}, vertexIndex: {vertexIndex})')
                    res = False
                else:
                    return False

    for matrixIndex, geoWeightData in softModData.localWeightList.items():
        for geoIndex, weights in geoWeightData.items():
            for vertexIndex, weight in weights.items():
                if int(vertexIndex) >= vertexCountData[geoIndex]:
                    if verbose:
                        om2.MGlobal.displayInfo(f'Using an index larger than the number of vertices in weightList.(inputDataIndex: {matrixIndex}, geoIndex: {geoIndex}, vertexIndex: {vertexIndex})')
                        res = False
                    else:
                        return False

    return res