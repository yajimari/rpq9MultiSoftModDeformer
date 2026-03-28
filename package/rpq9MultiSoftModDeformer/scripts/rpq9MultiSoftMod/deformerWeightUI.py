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

from .constants import PLUGIN_NAME, PLUGIN_NODE_NAME
from .utils import getDeformerGeometries, getMultiSoftModNodes
from .deformerWeight import (
    copyVertexWeightAttribute,
    stashWeightAttr,
    unstashWeightAttr,
    getMultiSoftModInputDataIndex,
    saveVertexWeightData,
    loadVertexWeightDataToNode
)


UI_NAME = 'rpq9MultiSoftMod_DeformerWeightUI'
DEFORMER_MENU_NAME = 'rpq9MultiSoftMod_DeformerWeightUI_DeformerMenu'
LOCAL_INDEX_MENU_NAME = 'rpq9MultiSoftMod_DeformerWeightUI_LocalIndexMenu'


#---UI commands
def getCurrentDeformer():
    currentDeformer = cmds.optionMenu(DEFORMER_MENU_NAME, q=True, v=True)
    if not currentDeformer:
        raise RuntimeError('No set deformer.')
    return currentDeformer


def selectDeformerGeometries():
    currentDeformer = getCurrentDeformer()
    geos = getDeformerGeometries(currentDeformer)
    cmds.select(geos)


def setDeformerMenu():
    cmds.optionMenu(DEFORMER_MENU_NAME, e=True, dai=True)
    sel = cmds.ls(sl=True)
    if not sel:
        raise RuntimeError('No selected.')
    multiSoftMods = getMultiSoftModNodes(sel[0])
    if not multiSoftMods:
        raise RuntimeError(f'Not found "{PLUGIN_NODE_NAME}."')
    for multiSoftMod in multiSoftMods:
        cmds.menuItem(label=multiSoftMod, p=DEFORMER_MENU_NAME)

    setLocalIndexMenu()


def setLocalIndexMenu():
    cmds.optionMenu(LOCAL_INDEX_MENU_NAME, e=True, dai=True)
    currentDeformer = getCurrentDeformer()
    mi = cmds.getAttr(f'{currentDeformer}.inputData', mi=True) or []
    for index in mi:
        cmds.menuItem(label=index, p=LOCAL_INDEX_MENU_NAME)


def filePathFromDialog(dialogMode:int) -> str:
    DIALOG_CAPTION = ('Save As', 'Open')
    fileDialogStyle = cmds.optionVar(q='FileDialogStyle')
    filePath = cmds.fileDialog2(ff='JSON Files(*.json)', ds=fileDialogStyle, cap=DIALOG_CAPTION[dialogMode], fm=dialogMode)
    return filePath[0] if filePath else None


def exportVertexWeightDataCB():
    currentDeformer = getCurrentDeformer()
    filePath = filePathFromDialog(dialogMode=0)
    if filePath is None:
        return
    saveVertexWeightData(currentDeformer, filePath)


def importVertexWeightDataCB():
    currentDeformer = getCurrentDeformer()
    filePath = filePathFromDialog(dialogMode=1)
    if filePath is None:
        return
    loadVertexWeightDataToNode(currentDeformer, filePath)


def copyVertexWeightAttributeCB(reverse:bool=False):
    currentDeformer = getCurrentDeformer()
    currentIndex = cmds.optionMenu(LOCAL_INDEX_MENU_NAME, q=True, v=True)

    copyAttributes = [f'{currentDeformer}.weightList', f'{currentDeformer}.inputData[{currentIndex}].localWeightList']
    if reverse:
        copyAttributes.reverse()
    copyVertexWeightAttribute(*copyAttributes)


def setLocalIndexMenuValueFromSel():
    sel = cmds.ls(sl=True)
    if not sel:
        raise RuntimeError('No selected.')

    currentDeformer = getCurrentDeformer()
    index = getMultiSoftModInputDataIndex(sel[0], currentDeformer)
    if index >= 0:
        cmds.optionMenu(LOCAL_INDEX_MENU_NAME, e=True, v=index)


def clearWeightAttributeValue():
    currentDeformer = getCurrentDeformer()
    cmds.removeMultiInstance(f'{currentDeformer}.weightList', all=True)


def stashWeightAttrCB():
    currentDeformer = getCurrentDeformer()
    stashWeightAttr(currentDeformer)


def unstashWeightAttrCB():
    currentDeformer = getCurrentDeformer()
    unstashWeightAttr(currentDeformer)


def showWindow():
    if cmds.window(UI_NAME, ex=True):
        cmds.deleteUI(UI_NAME)

    window = cmds.window(UI_NAME, t=UI_NAME)
    masterForm = cmds.formLayout()
    column = cmds.columnLayout(adj=True)
    cmds.separator(h=10, st='none')
    cmds.rowLayout(numberOfColumns=3, adjustableColumn3=2)
    cmds.text(l='deformer: ')
    deformerMenu = cmds.optionMenu(DEFORMER_MENU_NAME, cc=lambda *args: setLocalIndexMenu())
    cmds.button(l='Set', w=50, c=lambda *args: setDeformerMenu())
    cmds.setParent('..')

    cmds.separator(h=5, st='none')
    cmds.button(l='Select Geometry', h=30, c=lambda *args: selectDeformerGeometries())
    cmds.separator(h=10, st='single')

    cmds.rowLayout(numberOfColumns=3, adjustableColumn3=2)

    cmds.columnLayout(adj=True)
    cmds.text(l='weight')
    cmds.button(l='clear Weight', c=lambda *args: clearWeightAttributeValue())
    cmds.setParent('..')

    cmds.columnLayout(adj=True)
    cmds.button(l='<-', c=lambda *args: copyVertexWeightAttributeCB(True))
    cmds.button(l='->', c=lambda *args: copyVertexWeightAttributeCB(False))
    cmds.setParent('..')

    cmds.columnLayout(adj=True)
    cmds.rowLayout(numberOfColumns=2, adjustableColumn2=2)
    cmds.text(l='localIndex: ')
    localIndexMenu = cmds.optionMenu(LOCAL_INDEX_MENU_NAME)
    cmds.setParent('..')
    cmds.button(l='Set Sel', c=lambda *args: setLocalIndexMenuValueFromSel())
    cmds.setParent('..')
    cmds.setParent('..')

    cmds.separator(h=10, st='none')
    cmds.setParent('..')

    form = cmds.formLayout()
    stashButton = cmds.button(l='stash Weight', c=lambda *args: stashWeightAttrCB())
    unstashButton = cmds.button(l='unstash weight', c=lambda *args: unstashWeightAttrCB())
    importButton = cmds.button(l='import', c=lambda *args: importVertexWeightDataCB())
    exportButton = cmds.button(l='export', c=lambda *args: exportVertexWeightDataCB())
    cmds.formLayout(form, e=True,
        af=[
            [stashButton, 'top', 0],
            [stashButton, 'left', 0],
            [unstashButton, 'top', 0],
            [unstashButton, 'right', 0],
            [importButton, 'bottom', 0],
            [importButton, 'left', 0],
            [exportButton, 'bottom', 0],
            [exportButton, 'right', 0]
        ],
        ap=[
            [stashButton, 'right', 2, 50],
            [stashButton, 'bottom', 2, 50],
            [unstashButton, 'left', 2, 50],
            [unstashButton, 'bottom', 2, 50],
            [importButton, 'top', 2, 50],
            [importButton, 'right', 2, 50],
            [exportButton, 'top', 2, 50],
            [exportButton, 'left', 2, 50]
        ]
    )
    cmds.setParent('..')

    cmds.formLayout(masterForm, e=True,
        af=[
            [column, 'top', 0],
            [column, 'left', 0],
            [column, 'right', 0],
            [form, 'left', 0],
            [form, 'right', 0],
            [form, 'bottom', 0]
        ],
        an=[[column, 'bottom']],
        ac=[[form, 'top', 0, column]]
    )
    cmds.setParent('..')
    cmds.showWindow(window)