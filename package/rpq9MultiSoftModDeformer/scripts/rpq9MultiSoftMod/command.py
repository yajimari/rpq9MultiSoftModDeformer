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

import maya.internal.common.cmd.deformer as cmddeformer
import maya.internal.common.utils.ui as ui_utils

from .model import FalloffMode
from .utils import createRpq9MultiSoftMod


class Command(cmddeformer.Command):
    def __init__(self):
        super(Command, self).__init__()

        self.commandName = 'Rpq9MultiSoftMod'
        self.commandHelpTag	= 'Rpq9MultiSoftMod'
        self.commandDescription = 'Compute multiple softMod in one node.'
        self.commandTitle = 'Rpq9MultiSoftMod Option'
        self.optionVarCategory = 'Deform.Rpq9MultiSoftMod'
        self.optionVarPrefix = 'Rpq9MultiSoftMod'
        self.optionVarDefaults.update( {
            'falloffMode': FalloffMode.kSmooth.value,
            'controllerNum': 5
            } )

    @classmethod
    def command(cls, falloffMode=FalloffMode.kSmooth.value, controllerNum=5, **kwargs):
        returnedNodes = list()

        kw = Command.getDeformerCommandArgs(**kwargs)

        sel = cmds.ls(os=True)
        result = createRpq9MultiSoftMod(sel, controllerNum, **kw)

        if result is not None:
            deformerNode = result[0]
            cmds.setAttr(f'{deformerNode}.falloffMode', falloffMode)
            returnedNodes.append(deformerNode)

        return cls.finalizeCommand(returnedNodes)

    def addBasicDeformerDialogWidgets(self):
        widgetDict = {}

        with ui_utils.AttributeLayoutManager():
            options=[
                ('None', FalloffMode.kNone.value),
                ('Smooth', FalloffMode.kSmooth.value),
                ('Spline', FalloffMode.kSpline.value),
                ('SmoothStep', FalloffMode.kSmoothStep.value),
                ('EaseInOut', FalloffMode.kEaseInOut.value),
                ]
            widget, lookup = ui_utils.createOptionMenu('Falloff Mode', options=options)
            self.optionMenuGrp_labelToEnum['falloffMode'] = lookup
            widgetDict['falloffMode'] = (cmds.optionMenuGrp, widget)

            widget = ui_utils.intSliderGrp('Controller Num', 0.0, 20.0, ulim=False)
            widgetDict['controllerNum'] = (cmds.intSliderGrp, widget)

        return widgetDict