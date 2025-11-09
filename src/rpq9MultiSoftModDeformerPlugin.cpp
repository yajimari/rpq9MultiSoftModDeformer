/*
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
*/
#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>
#include "rpq9MultiSoftModDeformerNode.h"
#include "embeddedMel.h"
#include "embeddedKernel.h"


static MStatus evalEmbeddedMel(){
    MString code;
    code.setUTF8(kEmbeddedMel);
    return MGlobal::executeCommand(code, /*display*/false, /*undo*/false);
}


MStatus initializePlugin( MObject obj ){
    MStatus status;
    MFnPlugin plugin( obj, "Ryoya Yajima", "1.0", "Any");
    status = plugin.registerNode(   Rpq9MultiSoftModDeformer::kPluginNodeName, Rpq9MultiSoftModDeformer::id, Rpq9MultiSoftModDeformer::creator,
                                    Rpq9MultiSoftModDeformer::initialize, MPxNode::kDeformerNode );
    if(!status)  return status;
    MString nodeClassName(Rpq9MultiSoftModDeformer::kPluginNodeName);
    MString registrantId(Rpq9MultiSoftModGPUDeformer::kRegistrantId);
    MGPUDeformerRegistry::registerGPUDeformerCreator(
        nodeClassName,
        registrantId,
        Rpq9MultiSoftModGPUDeformer::getGPUDeformerInfo());

    Rpq9MultiSoftModDeformer::kernelSource = kEmbeddedKernel;
    Rpq9MultiSoftModDeformer::kernelProgramName = "rpq9MultiSoftModDeformerKernels";

    status = evalEmbeddedMel();
    return status;
}


MStatus uninitializePlugin( MObject obj){
    MStatus status;
    MString nodeClassName(Rpq9MultiSoftModDeformer::kPluginNodeName);
    MString registrantId(Rpq9MultiSoftModGPUDeformer::kRegistrantId);
    MGPUDeformerRegistry::deregisterGPUDeformerCreator(
        nodeClassName,
        registrantId);

    MFnPlugin plugin(obj);
    status = plugin.deregisterNode(Rpq9MultiSoftModDeformer::id);

    return status;
}