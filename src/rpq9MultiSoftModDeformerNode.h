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
#pragma once
#include<array>
#include<map>
#include <maya/MPxDeformerNode.h>
#include <maya/MItGeometry.h>

#include <maya/MStatus.h>
#include <maya/MTypeId.h>
#include <maya/MPlug.h>

#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>

#include <maya/MMatrix.h>
#include <maya/MString.h>

#include <maya/MPxGPUStandardDeformer.h>
#include <maya/MGPUDeformerRegistry.h>
#include <maya/MOpenCLBuffer.h>
#include <maya/MOpenCLAutoPtr.h>
#include <maya/MOpenCLUtils.h>

#include "nodeRegisterData.h"


struct SoftModData {
    std::vector<float> localEnvelopeValues;
    std::vector<cl_float4> centerPositions;
    std::vector<cl_float4> translations;
    std::vector<cl_float4> quaternions;
    std::vector<cl_float4> scales;
    std::vector<cl_float4> shears;
    std::vector<float> falloffRadiusValues;
};


class Rpq9MultiSoftModDeformer : public MPxDeformerNode
{
public:
    Rpq9MultiSoftModDeformer();
    ~Rpq9MultiSoftModDeformer() override;

    static  void*       creator();
    static  MStatus     initialize();

    MStatus prepareDeform(MDataBlock& block,
                          unsigned int multiIndex,
                          unsigned int changeFlags) override;

    MStatus preEvaluation(  const MDGContext& context,
                            const MEvaluationNode& evaluationNode) override;

    MStatus	deform(MDataBlock& block,
                    MItGeometry& iter,
                    const MMatrix& mat,
                    unsigned int	multiIndex) override;

    MObject& accessoryAttribute() const override;
    MStatus accessoryNodeSetup(MDagModifier& cmd) override;

    static unsigned int getInputDataData(   MDataBlock& block,
                                            SoftModData& data,
                                            MStatus& status);

    static unsigned int getInputDataData(   MDataBlock& block,
                                            SoftModData& data,
                                            std::vector<float>& localWeightValues,
                                            unsigned int vertexNum,
                                            unsigned int multiIndex,
                                            MStatus& status);

public:
    static MTypeId id;

    static MObject inputData;
    static MObject localEnvelope;
    static MObject centerMatrix;
    static MObject modifyMatrix;
    static MObject falloffRadius;
    static MObject localWeightList;
    static MObject localWeights;
    static MObject falloffMode;

    static MString kernelSource;
    static MString kernelProgramName;
    static std::array<MString, 4> kernelNames;

    static const char* kPluginNodeName;

    enum FalloffMode {
        kNone = 0,
        kLinear = 1,
        kSmooth = 2,
        kEaseInOut = 3,
    };
private:
    std::map<unsigned int, SoftModData> softModDataCache;
    std::map<unsigned int, std::vector<float>> localWeightValuesCache;
    bool isLocalWeightDirty = false;
};



class Rpq9MultiSoftModGPUDeformer : public MPxGPUStandardDeformer
{
public:
    Rpq9MultiSoftModGPUDeformer();
    ~Rpq9MultiSoftModGPUDeformer() override;

    MPxGPUDeformer::DeformerStatus evaluate(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& outputPlug, const MPlugArray& inputPlugs, const MGPUDeformerData& inputData, MGPUDeformerData& outputData) override;
    void terminate() override;

    static MGPUDeformerRegistrationInfo* getGPUDeformerInfo();
    static bool validateNodeInGraph(MDataBlock& block, const MEvaluationNode&, const MPlug& plug, MStringArray* messages);
    static bool validateNodeValues(MDataBlock& block, const MEvaluationNode&, const MPlug& plug, MStringArray* messages);
    bool passThroughWithZeroEnvelope() const override;

private:
    void prepareKernels();

public:
    static const char* kRegistrantId;

private:
    MOpenCLBuffer localEnvelopeBuffer;
    MOpenCLBuffer centerPositionsBuffer;
    MOpenCLBuffer translationsBuffer;
    MOpenCLBuffer quaternionsBuffer;
    MOpenCLBuffer scalesBuffer;
    MOpenCLBuffer shearsBuffer;
    MOpenCLBuffer localWeightBuffer;
    MOpenCLBuffer falloffRadiusValuesBuffer;
    SoftModData softModDataCache;

    std::vector<MAutoCLKernel> kernelInfoArray;
};



class Rpq9MultiSoftModDeformerNodeGPUDeformerInfo : public MGPUDeformerRegistrationInfo
{
public:
    Rpq9MultiSoftModDeformerNodeGPUDeformerInfo();
    ~Rpq9MultiSoftModDeformerNodeGPUDeformerInfo() override;

    MPxGPUDeformer* createGPUDeformer() override;

    bool validateNodeInGraph(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages) override;
    bool validateNodeValues(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages) override;
};