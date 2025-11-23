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
#include <cassert>
#include <cmath>
#include <algorithm>

#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>

#include <maya/MFnMatrixData.h>
#include <maya/MPointArray.h>
#include <maya/MVector.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <maya/MOpenCLInfo.h>

#include "rpq9MultiSoftModDeformerNode.h"


namespace {
    using FalloffWeightFunc = float(*)(float, float);

    //--- interporation functions
    template <typename T>
    T remap1to0(T value, T maxValue) noexcept {
        static_assert(std::is_floating_point_v<T>, "remap1to0 requires an floating point type");

        if (maxValue == T{ 0 }) return T{ 0 };
        return std::max(T{ 0 }, T{ 1 } - value / maxValue);
    }


    template <typename T>
    T smoothQuadratic(T value) noexcept {
        static_assert(std::is_floating_point_v<T>, "smoothQuadratic requires a floating point type");
        return value * value * (T{ 3 } - T{ 2 } *value);
    }


    template <typename T>
    T easeInOutQuadratic(T value) noexcept {
        static_assert(std::is_floating_point_v<T>, "easeInOutQuadratic requires a floating point type");
        if (value < T{ 0.5 }) {
            return T{ 2 } *value * value;
        }
        else {
            const T t = T{ 2 } - T{ 2 } *value;
            return T{ 1 } - (t * t) / T{ 2 };
        }
    }
    //---

    //column-major matrix
    struct alignas(64) Mat4 {
        cl_float4 c0, c1, c2, c3;
    };


    constexpr inline cl_float4 operator+(const cl_float4& a, const cl_float4& b) noexcept {
        return { a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w };
    }

    constexpr inline cl_float4 operator-(const cl_float4& a, const cl_float4& b) noexcept {
        return { a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w };
    }

    constexpr inline cl_float4 operator*(const cl_float4& a, float s) noexcept {
        return { a.x * s, a.y * s, a.z * s, a.w * s };
    }

    constexpr inline cl_float4 operator*(float s, const cl_float4& a) noexcept {
        return a * s;
    }

    constexpr inline cl_float4 operator/(const cl_float4& a, float s) noexcept {
        return { a.x / s, a.y / s, a.z / s, a.w / s };
    }


    constexpr inline cl_float4 multiply(const Mat4& mat, cl_float4 vec) {
        return mat.c0 * vec.x + mat.c1 * vec.y + mat.c2 * vec.z + mat.c3 * vec.w;
    }

    inline float length3(cl_float4 vec) noexcept {
        return std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    }

    inline float length(cl_float4 vec) noexcept {
        return std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z + vec.w * vec.w);
    }


    template <typename T>
    constexpr inline cl_float4 to_cl_float4(T v0, T v1, T v2, T v3) noexcept {
        cl_float4 res;
        res.x = static_cast<float>(v0);
        res.y = static_cast<float>(v1);
        res.z = static_cast<float>(v2);
        res.w = static_cast<float>(v3);
        return res;
    }


    constexpr inline Mat4 identity() noexcept {
        return Mat4{
            cl_float4{1.0f, 0.0f, 0.0f, 0.0f},
            cl_float4{0.0f, 1.0f, 0.0f, 0.0f},
            cl_float4{0.0f, 0.0f, 1.0f, 0.0f},
            cl_float4{0.0f, 0.0f, 0.0f, 1.0f}};
    }


    constexpr inline Mat4 multiply(const Mat4& a, const Mat4& b) noexcept {
        return Mat4{
            a.c0 * b.c0.x + a.c1 * b.c0.y + a.c2 * b.c0.z + a.c3 * b.c0.w,
            a.c0 * b.c1.x + a.c1 * b.c1.y + a.c2 * b.c1.z + a.c3 * b.c1.w,
            a.c0 * b.c2.x + a.c1 * b.c2.y + a.c2 * b.c2.z + a.c3 * b.c2.w,
            a.c0 * b.c3.x + a.c1 * b.c3.y + a.c2 * b.c3.z + a.c3 * b.c3.w};
    }


    constexpr inline cl_float4 quatMultiply(const cl_float4& a, const cl_float4& b) noexcept {
        return {
                a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
                a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
                a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
                a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
    }


    constexpr inline void setTranslationToMatrix(Mat4& mat, const cl_float4 translate) noexcept {
        mat.c3.x = translate.x;
        mat.c3.y = translate.y;
        mat.c3.z = translate.z;
    }


    inline void setQuatToMatrix(Mat4& mat, const cl_float4 quat) noexcept {
        const float len = length(quat);
        float xn = 0.0f, yn = 0.0f, zn = 0.0f, wn = 1.0f;
        if (len > 0.0f) {
            const float inv = 1.0f / len;
            xn = quat.x * inv; yn = quat.y * inv; zn = quat.z * inv; wn = quat.w * inv;
        }

        const float xx = xn * xn, yy = yn * yn, zz = zn * zn;
        const float xy = xn * yn, xz = xn * zn, yz = yn * zn;
        const float wx = wn * xn, wy = wn * yn, wz = wn * zn;

        const float m00 = 1.0f - 2.0f * (yy + zz);
        const float m01 = 2.0f * (xy - wz);
        const float m02 = 2.0f * (xz + wy);

        const float m10 = 2.0f * (xy + wz);
        const float m11 = 1.0f - 2.0f * (xx + zz);
        const float m12 = 2.0f * (yz - wx);

        const float m20 = 2.0f * (xz - wy);
        const float m21 = 2.0f * (yz + wx);
        const float m22 = 1.0f - 2.0f * (xx + yy);

        mat.c0 = cl_float4{ m00, m10, m20, 0.0f };
        mat.c1 = cl_float4{ m01, m11, m21, 0.0f };
        mat.c2 = cl_float4{ m02, m12, m22, 0.0f };
    }

 
    constexpr inline void setScaleToMatrix(Mat4& mat, const cl_float4 scale) noexcept {
        mat.c0.x = scale.x;
        mat.c1.y = scale.y;
        mat.c2.z = scale.z;
    }


    constexpr inline void setShearToMatrix(Mat4& mat, const cl_float4 shear) noexcept {
        const float shXY = shear.x, shXZ = shear.y, shYZ = shear.z;
        mat.c1.x = shXY;
        mat.c2.x = shXZ;
        mat.c2.y = shYZ;
    }


    inline Mat4 composeMatrix(cl_float4 translate,
        cl_float4 quaternion,
        cl_float4 scale,
        cl_float4 shear) noexcept
    {
        Mat4 TR = identity(), S = identity(), Sh = identity();
        setTranslationToMatrix(TR, translate);
        setQuatToMatrix(TR, quaternion);
        setScaleToMatrix(S, scale);
        setShearToMatrix(Sh, shear);
        return multiply(multiply(TR, S), Sh);
    }


    constexpr inline float dot(cl_float4 a,  cl_float4 b) noexcept {
        return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
    }


    constexpr inline cl_float4 negate( cl_float4 vec) noexcept {
        return cl_float4({ -vec.x, -vec.y, -vec.z, -vec.w });
    }


    constexpr float EPS_F = 1e-6f;

    inline cl_float4 normalize( cl_float4 q) noexcept {
        float n2 = dot(q, q);
        if (n2 < EPS_F * EPS_F) return cl_float4({ 0.f, 0.f, 0.f, 1.f });
        float inv = 1.0f / std::sqrt(n2);
        return q * inv;
    }


    inline float safe_acosf(float value) noexcept {
        value = std::max(-1.0f, std::min(1.0f, value));
        return std::acos(value);
    }


    cl_float4 lerp(cl_float4 a, cl_float4 b, float t) noexcept {
        cl_float4 r;
        r.x = std::fma(t, b.x - a.x, a.x);
        r.y = std::fma(t, b.y - a.y, a.y);
        r.z = std::fma(t, b.z - a.z, a.z);
        r.w = std::fma(t, b.w - a.w, a.w);
        return r;
    }


    cl_float4 slerp(cl_float4 a, cl_float4 b, float t){
        a = normalize(a);
        b = normalize(b);

        float cos_omega = dot(a, b);
        if (cos_omega < 0.0f) {
            cos_omega = -cos_omega;
            b = negate(b);
        }

        float theta = safe_acosf(cos_omega);

        if (std::fabs(theta) < EPS_F) {
            return normalize(a * (1.0f - t) +  b * t);
        }

        float sTheta = std::sin(theta);
        if (std::fabs(sTheta) < EPS_F) {
            return normalize(a * (1.0f - t) + b * t);
        }

        float s0 = std::sin((1.0f - t) * theta) / sTheta;
        float s1 = std::sin(t * theta) / sTheta;

        return normalize(a * s0 + b * s1);
    }
}


// -----------------------------------------------------------------------------
// Rpq9MultiSoftModGPUDeformer
// -----------------------------------------------------------------------------
MTypeId Rpq9MultiSoftModDeformer::id(RPQ9_MULTI_SOFT_MOD_DEFORMER_NODE_ID);
MObject Rpq9MultiSoftModDeformer::inputData;
MObject Rpq9MultiSoftModDeformer::localEnvelope;
MObject Rpq9MultiSoftModDeformer::centerMatrix;
MObject Rpq9MultiSoftModDeformer::modifyMatrix;
MObject Rpq9MultiSoftModDeformer::falloffRadius;
MObject Rpq9MultiSoftModDeformer::falloffMode;
MString Rpq9MultiSoftModDeformer::kernelSource;
MString Rpq9MultiSoftModDeformer::kernelProgramName;
std::array<MString, 4> Rpq9MultiSoftModDeformer::kernelNames = {"noneSoftMod", "linearSoftMod", "smoothSoftMod", "easeInOutSoftMod"};
const char* Rpq9MultiSoftModDeformer::kPluginNodeName = RPQ9_MULTI_SOFT_MOD_DEFORMER_NODE_NAME;

Rpq9MultiSoftModDeformer::Rpq9MultiSoftModDeformer() {}
Rpq9MultiSoftModDeformer::~Rpq9MultiSoftModDeformer() {}

void* Rpq9MultiSoftModDeformer::creator(){
    return new Rpq9MultiSoftModDeformer();
}


MStatus Rpq9MultiSoftModDeformer::initialize(){
    MFnTypedAttribute tAttr;
    MFnNumericAttribute nAttr;
    MFnEnumAttribute eAttr;
    MFnCompoundAttribute cAttr;

    localEnvelope = nAttr.create("localEnvelope", "le", MFnNumericData::kFloat, 1.0f);
    nAttr.setKeyable(true);
    nAttr.setReadable(false);
    nAttr.setMin(0.0);
    nAttr.setMax(1.0);
    nAttr.setDisconnectBehavior(MFnNumericAttribute::kNothing);

    centerMatrix = tAttr.create("centerMatrix", "cm", MFnData::kMatrix, MFnMatrixData().create());
    tAttr.setKeyable(true);
    tAttr.setReadable(false);
    tAttr.setDisconnectBehavior(MFnNumericAttribute::kNothing);

    modifyMatrix = tAttr.create("modifyMatrix", "mm", MFnData::kMatrix, MFnMatrixData().create());
    tAttr.setKeyable(true);
    tAttr.setReadable(false);
    tAttr.setDisconnectBehavior(MFnNumericAttribute::kNothing);

    falloffRadius = nAttr.create("falloffRadius", "fr", MFnNumericData::kFloat, 1.0f);
    nAttr.setKeyable(true);
    nAttr.setReadable(false);
    nAttr.setDisconnectBehavior(MFnNumericAttribute::kNothing);

    inputData = cAttr.create("inputData", "id");
    cAttr.setKeyable(true);
    cAttr.setReadable(false);
    cAttr.setArray(true);
    cAttr.setDisconnectBehavior(MFnCompoundAttribute::kDelete);
    cAttr.addChild(localEnvelope);
    cAttr.addChild(centerMatrix);
    cAttr.addChild(modifyMatrix);
    cAttr.addChild(falloffRadius);

    falloffMode = eAttr.create("falloffMode", "fm", Rpq9MultiSoftModDeformer::kSmooth);
    eAttr.setKeyable(true);
    eAttr.setReadable(false);
    eAttr.addField("None", Rpq9MultiSoftModDeformer::kNone);
    eAttr.addField("Linear", Rpq9MultiSoftModDeformer::kLinear);
    eAttr.addField("Smooth", Rpq9MultiSoftModDeformer::kSmooth);
    eAttr.addField("EaseInOut", Rpq9MultiSoftModDeformer::kEaseInOut);

    addAttribute(inputData);
    addAttribute(falloffMode);

    attributeAffects(inputData, outputGeom);
    attributeAffects(falloffMode, outputGeom);

    return MStatus::kSuccess;
}


MStatus Rpq9MultiSoftModDeformer::prepareDeform(MDataBlock& block, unsigned int multiIndex, unsigned int changeFlags){
    return MS::kSuccess;
}


unsigned int Rpq9MultiSoftModDeformer::getInputDataData(MDataBlock& block,
                                                        std::vector<float>& localEnvelopeValues,
                                                        std::vector<cl_float4>& centerPositions,
                                                        std::vector<cl_float4>& translations,
                                                        std::vector<cl_float4>& quaternions,
                                                        std::vector<cl_float4>& scales,
                                                        std::vector<cl_float4>& shears,
                                                        std::vector<float>& falloffRadiusValues,
                                                        MStatus& status) 
{
    MArrayDataHandle inputDataArrayDataHandle = block.inputArrayValue(Rpq9MultiSoftModDeformer::inputData, &status);
    if (MS::kSuccess != status) return 0;

    unsigned int inputDataNum = inputDataArrayDataHandle.elementCount();
    if (inputDataNum <= 0) return 0;

    localEnvelopeValues.resize(inputDataNum);
    centerPositions.resize(inputDataNum);
    translations.resize(inputDataNum);
    quaternions.resize(inputDataNum);
    scales.resize(inputDataNum);
    shears.resize(inputDataNum);
    falloffRadiusValues.resize(inputDataNum);


    for (unsigned int i = 0; i < inputDataNum; ++i) {
        MDataHandle inputDataDataHandle = inputDataArrayDataHandle.inputValue();
        MDataHandle childHandle = inputDataDataHandle.child(Rpq9MultiSoftModDeformer::localEnvelope);
        localEnvelopeValues[i] = childHandle.asFloat();

        childHandle = inputDataDataHandle.child(Rpq9MultiSoftModDeformer::centerMatrix);
        MMatrix centerMatrixValue = MFnMatrixData(childHandle.data()).matrix();
        centerPositions[i] = to_cl_float4(centerMatrixValue(3, 0), centerMatrixValue(3, 1), centerMatrixValue(3, 2), 1.0);

        childHandle = inputDataDataHandle.child(Rpq9MultiSoftModDeformer::modifyMatrix);
        MMatrix matrix = MFnMatrixData(childHandle.data()).matrix();

        MMatrix localCenterMatrixInverseValue = centerMatrixValue.inverse();
        localCenterMatrixInverseValue(3, 0) = 0.0;
        localCenterMatrixInverseValue(3, 1) = 0.0;
        localCenterMatrixInverseValue(3, 2) = 0.0;

        matrix(3, 0) -= centerPositions[i].x;
        matrix(3, 1) -= centerPositions[i].y;
        matrix(3, 2) -= centerPositions[i].z;

        MMatrix transformMatrix = localCenterMatrixInverseValue * matrix;
        MTransformationMatrix transformationMatrix(transformMatrix);

        double temp3[3] = { 0.0f, 0.0f, 0.0f };
        MVector translate = transformationMatrix.getTranslation(MSpace::kTransform);
        translate.get(temp3);
        translations[i] = to_cl_float4(temp3[0], temp3[1], temp3[2], 1.0);

        transformationMatrix.getScale(temp3, MSpace::kTransform);
        scales[i] = to_cl_float4(temp3[0], temp3[1], temp3[2], 1.0);

        transformationMatrix.getShear(temp3, MSpace::kTransform);
        shears[i] = to_cl_float4(temp3[0], temp3[1], temp3[2], 1.0);

        double temp4[4] = { 0.0f, 0.0f, 0.0f, 1.0f };
        transformationMatrix.getRotationQuaternion(temp4[0], temp4[1], temp4[2], temp4[3]);
        quaternions[i] = to_cl_float4(temp4[0], temp4[1], temp4[2], temp4[3]);

        childHandle = inputDataDataHandle.child(Rpq9MultiSoftModDeformer::falloffRadius);
        falloffRadiusValues[i] = childHandle.asFloat();

        inputDataArrayDataHandle.next();
    }
    status = MS::kSuccess;
    return inputDataNum;
}


MStatus Rpq9MultiSoftModDeformer::deform(MDataBlock& block, MItGeometry& iter, const MMatrix& /*m*/, unsigned int multiIndex){
    MStatus status;

    MDataHandle envDataHandle = block.inputValue(envelope, &status);
    if (MS::kSuccess != status) return status;
    float deformerEnvelope = envDataHandle.asFloat();

    std::vector<float> localEnvelopeValues;
    std::vector<cl_float4> centerPositions;
    std::vector<cl_float4> translations;
    std::vector<cl_float4> quaternions;
    std::vector<cl_float4> scales;
    std::vector<cl_float4> shears;
    std::vector<float> falloffRadiusValues;

    unsigned int inputDataNum = Rpq9MultiSoftModDeformer::getInputDataData(block, localEnvelopeValues, centerPositions, translations, quaternions, scales, shears, falloffRadiusValues, status);

    MDataHandle falloffModeDataHandle = block.inputValue(falloffMode, &status);
    if (MS::kSuccess != status) return status;
    short falloffModeValue = falloffModeDataHandle.asShort();

    FalloffWeightFunc falloffWeightFunc = nullptr;
    const Rpq9MultiSoftModDeformer::FalloffMode mode = static_cast<Rpq9MultiSoftModDeformer::FalloffMode>(falloffModeValue);
    switch (mode) {
        case Rpq9MultiSoftModDeformer::kNone:
            falloffWeightFunc = [](float, float) {return 1.0f;};
            break;
        case Rpq9MultiSoftModDeformer::kLinear:
            falloffWeightFunc = [](float distanceFromCenter, float falloffRadius) {return remap1to0(distanceFromCenter, falloffRadius);};
            break;
        case Rpq9MultiSoftModDeformer::kSmooth:
            falloffWeightFunc = [](float distanceFromCenter, float falloffRadius) {return smoothQuadratic(remap1to0(distanceFromCenter, falloffRadius));};
            break;
        case Rpq9MultiSoftModDeformer::kEaseInOut:
            falloffWeightFunc = [](float distanceFromCenter, float falloffRadius) {return easeInOutQuadratic(remap1to0(distanceFromCenter, falloffRadius));};
            break;
        default:
            falloffWeightFunc = [](float, float) {return 1.0f;};
            break;
    }

    unsigned int numWeights;
    const MIndexMapper& mapper = indexMapper(multiIndex);
    const float* vertWeights = envelopeWeights(block, multiIndex, &numWeights);

    MPointArray points;
    iter.allPositions(points);
    unsigned int affectCount = mapper.affectCount();
    assert(nullptr == vertWeights || numWeights == affectCount);
    assert(points.length() == affectCount);

    constexpr cl_float4 zeroVector = { 0.0f, 0.0f, 0.0f, 1.0f };
    constexpr cl_float4 oneVector = { 1.0f, 1.0f, 1.0f, 1.0f };

    auto deformOneVertex = [&](unsigned int aid) {const unsigned int uId = static_cast<unsigned int>(aid);
        cl_float4 clPt = to_cl_float4(points[uId].x, points[uId].y, points[uId].z, 1.0);
        cl_float4 addPt{ 0.0f, 0.0f, 0.0f, 0.0f };

        for(unsigned int j = 0; j < inputDataNum; ++j){
            cl_float4 distanceVector = centerPositions[j] - clPt;
            float distanceFromCenter = length3(distanceVector);
            if (distanceFromCenter > falloffRadiusValues[j]) continue;

            const float falloffWeight = falloffWeightFunc(distanceFromCenter, falloffRadiusValues[j]);
            const cl_float4 weightedTranslate  = lerp(zeroVector, translations[j],  falloffWeight);
            const cl_float4 weightedScale      = lerp(oneVector,  scales[j],        falloffWeight);
            const cl_float4 weightedShear      = lerp(zeroVector, shears[j],        falloffWeight);
            const cl_float4 weightedQuaternion = slerp(zeroVector, quaternions[j],  falloffWeight);
            const Mat4 transMat = composeMatrix(weightedTranslate, weightedQuaternion, weightedScale, weightedShear);

            cl_float4 localPt = clPt - centerPositions[j];
            localPt.w = 1.0f;

            const cl_float4 newPtBase = multiply(transMat, localPt);
            addPt = addPt + (newPtBase - localPt) * localEnvelopeValues[j];
        }

        const float weight = (vertWeights ? vertWeights[aid] * deformerEnvelope : deformerEnvelope);
        const MVector mDelta(static_cast<double>(addPt.x), static_cast<double>(addPt.y), static_cast<double>(addPt.z));
        points[uId] += weight * mDelta;
    };

    constexpr unsigned int kMinVertsForParallel = 256u;

    if(affectCount < kMinVertsForParallel) {
        for (unsigned int aid = 0; aid < affectCount; ++aid) {
            deformOneVertex(aid);
        }
    }else{
        tbb::parallel_for(tbb::blocked_range<unsigned int>(0, affectCount), [&](const tbb::blocked_range<unsigned int>& r){
                for (unsigned int aid = r.begin(); aid != r.end(); ++aid){
                    deformOneVertex(aid);
                }
            }
        );
    }

    iter.setAllPositions(points);
    return status;
}


MObject& Rpq9MultiSoftModDeformer::accessoryAttribute() const {
    return MObject::kNullObj;
}


MStatus Rpq9MultiSoftModDeformer::accessoryNodeSetup(MDagModifier& cmd){
    return MS::kSuccess;
}


// -----------------------------------------------------------------------------
// Rpq9MultiSoftModGPUDeformer
// -----------------------------------------------------------------------------
const char* Rpq9MultiSoftModGPUDeformer::kRegistrantId = RPQ9_MULTI_SOFT_MOD_GPU_DEFORMER_REGISTRANT_ID;

MGPUDeformerRegistrationInfo* Rpq9MultiSoftModGPUDeformer::getGPUDeformerInfo()
{
    static Rpq9MultiSoftModDeformerNodeGPUDeformerInfo theOne;
    return &theOne;
}

Rpq9MultiSoftModGPUDeformer::Rpq9MultiSoftModGPUDeformer() : MPxGPUStandardDeformer() {}

Rpq9MultiSoftModGPUDeformer::~Rpq9MultiSoftModGPUDeformer(){
    terminate();
}


bool Rpq9MultiSoftModGPUDeformer::validateNodeInGraph(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages){
    return true;
}


bool Rpq9MultiSoftModGPUDeformer::validateNodeValues(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages){
    return true;
}


void Rpq9MultiSoftModGPUDeformer::prepareKernels(){
    const std::size_t kernelNum = Rpq9MultiSoftModDeformer::kernelNames.size();
    if (kernelNum == kernelInfoArray.size()) return;
    kernelInfoArray.resize(kernelNum);
    for(std::size_t i=0; i<kernelNum; ++i){
        kernelInfoArray[i] = MOpenCLInfo::getOpenCLKernelFromString(Rpq9MultiSoftModDeformer::kernelSource,
                                                                    Rpq9MultiSoftModDeformer::kernelProgramName,
                                                                    Rpq9MultiSoftModDeformer::kernelNames[i]);
    }
}


MPxGPUDeformer::DeformerStatus Rpq9MultiSoftModGPUDeformer::evaluate(
    MDataBlock& block,
    const MEvaluationNode& evaluationNode,
    const MPlug& plug,
    const MPlugArray & inputPlugs,
    const MGPUDeformerData& inputData,
    MGPUDeformerData& outputData)
{
    MStatus status;
    MAutoCLEvent clEvent;
    MAutoCLEventList clEventList;

    MPxGPUDeformer::DeformerStatus deformerStatus = prepareEvaluation(block, evaluationNode, plug, inputPlugs, inputData, outputData, clEvent);
    if (deformerStatus != MPxGPUDeformer::kDeformerSuccess) return deformerStatus;
    
    prepareAffectMapBuffer();
    prepareWeightsBuffer(evaluationNode);

    std::vector<float> localEnvelopeValues;
    std::vector<cl_float4> centerPositions;
    std::vector<cl_float4> translations;
    std::vector<cl_float4> quaternions;
    std::vector<cl_float4> scales;
    std::vector<cl_float4> shears;
    std::vector<float> falloffRadiusValues;

    unsigned int inputDataNum = Rpq9MultiSoftModDeformer::getInputDataData(block, localEnvelopeValues, centerPositions, translations, quaternions, scales, shears, falloffRadiusValues, status);

    MDataHandle falloffModeDataHandle = block.inputValue(Rpq9MultiSoftModDeformer::falloffMode, &status);
    if (MS::kSuccess != status) return MPxGPUDeformer::kDeformerFailure;
    short falloffModeValue = falloffModeDataHandle.asShort();

    MOpenCLUtils::uploadToGPU(localEnvelopeValues, localEnvelopeBuffer, MOpenCLUtils::kBlocking);
    MOpenCLUtils::uploadToGPU(centerPositions, centerPositionsBuffer, MOpenCLUtils::kBlocking);
    MOpenCLUtils::uploadToGPU(translations, translationsBuffer, MOpenCLUtils::kBlocking);
    MOpenCLUtils::uploadToGPU(quaternions, quaternionsBuffer, MOpenCLUtils::kBlocking);
    MOpenCLUtils::uploadToGPU(scales, scalesBuffer, MOpenCLUtils::kBlocking);
    MOpenCLUtils::uploadToGPU(shears, shearsBuffer, MOpenCLUtils::kBlocking);
    MOpenCLUtils::uploadToGPU(falloffRadiusValues, falloffRadiusValuesBuffer, MOpenCLUtils::kBlocking);

    cl_int err = CL_SUCCESS;
    unsigned int count = affectCount();
    prepareKernels();
    MAutoCLKernel& currentKernel = kernelInfoArray[falloffModeValue];

    err = initializeOutputPositions(clEvent);
    if ( err != CL_SUCCESS ) return MPxGPUDeformer::kDeformerFailure;

    MAutoCLEvent syncInputEvent = clEvent;
    clEvent = MAutoCLEvent();
    clEventList.add(syncInputEvent);

    unsigned int parameterId = 0;
    err = MOpenCLUtils::setKernelArgBuffer(currentKernel.get(), parameterId++, inputPositions().buffer(), err);
    err = MOpenCLUtils::setKernelArgBuffer(currentKernel.get(), parameterId++, localEnvelopeBuffer, err);
    err = MOpenCLUtils::setKernelArgBuffer(currentKernel.get(), parameterId++, centerPositionsBuffer, err);
    err = MOpenCLUtils::setKernelArgBuffer(currentKernel.get(), parameterId++, translationsBuffer, err);
    err = MOpenCLUtils::setKernelArgBuffer(currentKernel.get(), parameterId++, quaternionsBuffer, err);
    err = MOpenCLUtils::setKernelArgBuffer(currentKernel.get(), parameterId++, scalesBuffer, err);
    err = MOpenCLUtils::setKernelArgBuffer(currentKernel.get(), parameterId++, shearsBuffer, err);
    err = MOpenCLUtils::setKernelArgBuffer(currentKernel.get(), parameterId++, falloffRadiusValuesBuffer, err);
    err = MOpenCLUtils::setKernelArgBuffer(currentKernel.get(), parameterId++, weightsBuffer(), err, hasEnvelopeWeights());
    err = MOpenCLUtils::setKernelArgBuffer(currentKernel.get(), parameterId++, affectMapBuffer(), err, !isIdentityMap());
    err = MOpenCLUtils::setKernelArgBuffer(currentKernel.get(), parameterId++, outputPositions().buffer(), err);
    err = MOpenCLUtils::setKernelArg<cl_float>(currentKernel.get(), parameterId++, envelope(), err);
    err = MOpenCLUtils::setKernelArg<cl_uint>(currentKernel.get(), parameterId++, count, err);
    err = MOpenCLUtils::setKernelArg<cl_uint>(currentKernel.get(), parameterId++, inputDataNum, err);

    size_t workGroupSize;
    size_t retSize;
    err = clGetKernelWorkGroupInfo(
        currentKernel.get(),
        MOpenCLInfo::getOpenCLDeviceId(),
        CL_KERNEL_WORK_GROUP_SIZE,
        sizeof(size_t),
        &workGroupSize,
        &retSize);
    MOpenCLInfo::checkCLErrorStatus(err);

    size_t localWorkSize = 256;
    if (retSize > 0) localWorkSize = workGroupSize;
    size_t globalWorkSize = ((count + localWorkSize - 1) / localWorkSize) * localWorkSize;

    err = clEnqueueNDRangeKernel(
        MOpenCLInfo::getMayaDefaultOpenCLCommandQueue(),
        currentKernel.get(),
        1,
        NULL,
        &globalWorkSize,
        &localWorkSize,
        clEventList.size(),
        clEventList.array(),
        clEvent.getReferenceForAssignment());

    if ( err != CL_SUCCESS ) return MPxGPUDeformer::kDeformerFailure;

    return finishEvaluation(clEvent, outputData);
}


void Rpq9MultiSoftModGPUDeformer::terminate(){
    localEnvelopeBuffer.reset();
    centerPositionsBuffer.reset();
    translationsBuffer.reset();
    quaternionsBuffer.reset();
    scalesBuffer.reset();
    shearsBuffer.reset();
    falloffRadiusValuesBuffer.reset();
    for(std::size_t i=0; i<kernelInfoArray.size(); ++i){
        MOpenCLInfo::releaseOpenCLKernel(kernelInfoArray[i]);
    }
    kernelInfoArray.clear();
    MPxGPUStandardDeformer::terminate();
}


// -----------------------------------------------------------------------------
// Rpq9MultiSoftModDeformerNodeGPUDeformerInfo
// -----------------------------------------------------------------------------

Rpq9MultiSoftModDeformerNodeGPUDeformerInfo::Rpq9MultiSoftModDeformerNodeGPUDeformerInfo() {}
Rpq9MultiSoftModDeformerNodeGPUDeformerInfo::~Rpq9MultiSoftModDeformerNodeGPUDeformerInfo() {}

MPxGPUDeformer* Rpq9MultiSoftModDeformerNodeGPUDeformerInfo::createGPUDeformer(){
    return new Rpq9MultiSoftModGPUDeformer();
}


bool Rpq9MultiSoftModDeformerNodeGPUDeformerInfo::validateNodeInGraph(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages){
    return Rpq9MultiSoftModGPUDeformer::validateNodeInGraph(block, evaluationNode, plug, messages);
}


bool Rpq9MultiSoftModDeformerNodeGPUDeformerInfo::validateNodeValues(MDataBlock& block, const MEvaluationNode& evaluationNode, const MPlug& plug, MStringArray* messages){
    return Rpq9MultiSoftModGPUDeformer::validateNodeValues(block, evaluationNode, plug, messages);
}