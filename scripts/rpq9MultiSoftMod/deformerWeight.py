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
import copy

import maya.cmds as cmds
import maya.api.OpenMaya as om2

from .constants import PLUGIN_NAME, PLUGIN_NODE_NAME
from .utils import isMultiSoftModNode
from .jsonIO import readJson, writeJson



STASH_ATTR_LIST_NAME = 'rpq9MultiSoftMod_stashWeightList'
STASH_ATTR_WEIGHT_NAME = 'rpq9MultiSoftMod_stashWeights'



def isVertexWeightAttriubte(plug:om2.MPlug) -> bool:
    if not plug.isArray or not plug.isCompound or plug.numChildren() != 1:
        return False
    
    elementPlug = plug.elementByLogicalIndex(0)
    childPlug = elementPlug.child(0)
    if not childPlug.isArray or childPlug.attribute().hasFn(om2.MFn.kNumericAttribute):
        return False
    return True


def copyVertexWeightAttribute(source:str, target:str) -> None:
    if not cmds.objExists(source) or not cmds.objExists(target):
        raise RuntimeError('Not exist argument attribute.')

    selList= om2.MSelectionList()
    selList.add(source)
    selList.add(target)

    sourcePlug = selList.getPlug(0)
    targetPlug = selList.getPlug(1)

    if isVertexWeightAttriubte(sourcePlug):
        raise TypeError(f'"{source}" is not vertex weight attriubte.')
    
    if isVertexWeightAttriubte(targetPlug):
        raise TypeError(f'"{target}" is not vertex weight attriubte.')

    sourceChild = sourcePlug.elementByLogicalIndex(0).child(0).name().split('.')[-1]
    targetChild = targetPlug.elementByLogicalIndex(0).child(0).name().split('.')[-1]

    cmds.removeMultiInstance(target, all=True)

    geoMultiIndices = cmds.getAttr(source, mi=True) or []

    for geoindex in geoMultiIndices:
        vertexIndices = cmds.getAttr(f'{source}[{geoindex}].{sourceChild}', mi=True) or []
        for vertexIndex in vertexIndices:
            value = cmds.getAttr(f'{source}[{geoindex}].{sourceChild}[{vertexIndex}]')

            cmds.setAttr(f'{target}[{geoindex}].{targetChild}[{vertexIndex}]', value)


def stashWeightAttr(deformer:str):
    if not cmds.attributeQuery(STASH_ATTR_LIST_NAME, node=deformer, ex=True):
        cmds.addAttr(deformer, ln=STASH_ATTR_LIST_NAME, numberOfChildren=1, at='compound', multi=True, disconnectBehaviour=0)
        cmds.addAttr(deformer, ln=STASH_ATTR_WEIGHT_NAME, at='float', parent=STASH_ATTR_LIST_NAME , multi=True, disconnectBehaviour=0)

    copyVertexWeightAttribute(f'{deformer}.weightList', f'{deformer}.{STASH_ATTR_LIST_NAME}')


def unstashWeightAttr(deformer:str):
    if not cmds.attributeQuery(STASH_ATTR_LIST_NAME, node=deformer, ex=True):
        raise RuntimeError('Not exist stash attr.')

    copyVertexWeightAttribute(f'{deformer}.{STASH_ATTR_LIST_NAME}', f'{deformer}.weightList')
    cmds.deleteAttr(f'{deformer}.{STASH_ATTR_LIST_NAME}')


def getConnectedMultiSoftModInputDataPlugs(node:str) -> list[om2.MPlug]:
    selList = om2.MSelectionList()
    selList.add(node)
    mobj = selList.getDependNode(0)
    iter = om2.MItDependencyGraph(mobj, om2.MFn.kPluginDeformerNode, om2.MItDependencyGraph.kDownstream, om2.MItDependencyGraph.kBreadthFirst)
    plugs = []
    while not iter.isDone():
        if isMultiSoftModNode(iter.currentNode()):
            plugs.append(iter.currentPlug())
        iter.next()
    return plugs


def getMultiSoftModInputDataIndex(node:str, deformer:str='') -> int:
    plugs = getConnectedMultiSoftModInputDataPlugs(node)
    if not plugs:
        return -1

    for plug in plugs:
        plugName = plug.name()

        if not 'inputData' in plugName:
            continue
        if deformer and plugName.split('.')[0] != deformer:
            continue

        return plug.parent().logicalIndex()
    return -1


def copyMultiSoftModDeformerWeights(source:str, target:str) -> None:
    if cmds.nodeType(source) != PLUGIN_NODE_NAME or cmds.nodeType(target) != PLUGIN_NODE_NAME:
        raise TypeError('Not valid deformer type.')
    source_geos = cmds.deformer(source, q=True, g=True)
    target_geos = cmds.deformer(target, q=True, g=True)
    if not len(source_geos) == len(target_geos):
        raise RuntimeError('Not same deformer geometry num.')
    
    stashWeightAttr(source)
    for source_geo, target_geo in zip(source_geos, target_geos):
        mi = cmds.getAttr(f'{source}.inputData', mi=True) or []
        for currentIndex in mi:
            copyVertexWeightAttribute(f'{source}.inputData[{currentIndex}].localWeightList', f'{source}.weightList')
            cmds.copyDeformerWeights(sourceDeformer=source, sourceShape=source_geo, destinationDeformer=target, destinationShape=target_geo, noMirror=True, surfaceAssociation='closestComponent')
            copyVertexWeightAttribute(f'{target}.weightList', f'{target}.inputData[{currentIndex}].localWeightList')

        unstashWeightAttr(source)
        cmds.copyDeformerWeights(sourceDeformer=source, sourceShape=source_geo, destinationDeformer=target, destinationShape=target_geo, noMirror=True, surfaceAssociation='closestComponent')


def getVertexWeightAttributeData(attr:str) -> dict[str, dict[str, float]] | None:
    selList= om2.MSelectionList()
    selList.add(attr)
    plug = selList.getPlug(0)

    if isVertexWeightAttriubte(plug):
        raise TypeError(f'"{attr}" is not vertex weight attriubte.')
    res = {}
    geoMultiIndices = plug.getExistingArrayAttributeIndices()
    for geoindex in geoMultiIndices:
        child = plug.elementByLogicalIndex(geoindex).child(0)
        res[geoindex] = {}
        vertexIndices = child.getExistingArrayAttributeIndices()
        for vertexIndex in vertexIndices:
            elementPlug = child.elementByLogicalIndex(vertexIndex)
            res[geoindex][vertexIndex] = elementPlug.asDouble()

    return res


def setVertexWeightAttributeData(attr:str, data:dict[str, dict[str, float]]) -> None:
    selList= om2.MSelectionList()
    selList.add(attr)
    plug = selList.getPlug(0)

    if isVertexWeightAttriubte(plug):
        raise TypeError(f'"{attr}" is not vertex weight attriubte.')

    for geoIndex, currentGeoData in data.items():
        geoPlug = plug.elementByLogicalIndex(int(geoIndex))
        childPlug = geoPlug.child(0)
        for vertexId, value in currentGeoData.items():
            vertexPlug = childPlug.elementByLogicalIndex(int(vertexId))
            cmds.setAttr(vertexPlug.name(), value)


def getVertexWeightData(node:str) -> dict | None:
    if cmds.nodeType(node) != PLUGIN_NODE_NAME:
        return None

    res = { 'name': node,
            'geometry': cmds.deformer(node, q=True, g=True),
            'geometryIndices': cmds.deformer(node, q=True, geometryIndices=True),
            'centerMatrices': {}}

    res['weightList'] = getVertexWeightAttributeData(f'{node}.weightList')

    res['localWeightList'] = {}
    inputDataMultiIndices = cmds.getAttr(f'{node}.inputData', mi=True)
    for inputDataMultiIndex in inputDataMultiIndices:
        res['localWeightList'][inputDataMultiIndex] = getVertexWeightAttributeData(f'{node}.inputData[{inputDataMultiIndex}].localWeightList')
        res['centerMatrices'][inputDataMultiIndex] = cmds.getAttr(f'{node}.inputData[{inputDataMultiIndex}].centerMatrix')

    return res


def setVertexWeightData(node:str, data:str) -> None:
    if cmds.nodeType(node) != PLUGIN_NODE_NAME:
        return

    setVertexWeightAttributeData(f'{node}.weightList', data['weightList'])

    for inputDataindex, currentIndexData in data['localWeightList'].items():
        setVertexWeightAttributeData(f'{node}.inputData[{inputDataindex}].localWeightList', currentIndexData)


def saveVertexWeightData(node:str, filePath:str) -> None:
    data = getVertexWeightData(node)
    writeJson(filePath, data)
    om2.MGlobal.displayInfo(f'Saved "{node}" vertex weights: {filePath}')


def loadVertexWeightData(filePath:str) -> None:
    data = readJson(filePath)
    deformerName = data['name']
    if not cmds.objExists(deformerName):
        raise RuntimeError(f'"{deformerName}" is not found.')
    setVertexWeightData(deformerName, data)


def loadVertexWeightDataToNode(node:str, filePath:str):
    data = readJson(filePath)
    if not cmds.objExists(node):
        raise RuntimeError(f'"{node}" is not found.')
    setVertexWeightData(node, data)


def applyDeformerData(filePsth:str):
    data = readJson(filePath)
    deformerName = data['name']
    geos = data['geometry']
    deformer = cmds.deformer(geos, type=PLUGIN_NODE_NAME, n=deformerName)

    setVertexWeightData(deformer, data)
    for index, value in data['centerMatrices'].items():
        cmds.setAttr(f'{deformer}.inputData[{index}].centerMatrix', value, type='matrix')
        cmds.setAttr(f'{deformer}.inputData[{index}].modifyMatrix', value, type='matrix')


def removeGeometryFromDefromerData(data:dict, removeGeo:str) -> int:
    geos = data['geometry']
    geoIndices = data['geometryIndices']
    geoIndex = geoIndices[geos.index(removeGeo)]

    del data['weightList'][str(geoIndex)]
    for inputDataindex, currentIndexData in data['localWeightList'].items():
        del currentIndexData[str(geoIndex)]

    return geoIndex