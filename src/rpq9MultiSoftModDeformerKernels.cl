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

 //--- interporation functions
inline float remap1to0(float value, float maxValue){
    // maxValue != 0 ? 1/maxValue : 0
    float inv = select(0.0f, native_recip(maxValue),
                       isnotequal(maxValue, 0.0f));

    float y = fma(-value, inv, 1.0f); // 1.0f - value * inv
    return fmax(0.0f, y);
}


inline float smoothQuadratic(float value){
    return value * value * (3.0f - 2.0f * value);
}


inline float easeInOutQuadratic(float value){
    float a = 2.0f * value * value;
    float t = 2.0f - 2.0f * value;
    float b = 1.0f - 0.5f * t * t;

    return select(a, b, isgreaterequal(value, 0.5f));
}
//-------------


// -------------------------------------------------------------
// column-major
// -------------------------------------------------------------

typedef struct __attribute__((aligned(16))){
    float4 c0, c1, c2, c3; // columns
} Mat4;


#define EPS_F (1e-6f)

inline float safe_acosf(float x){
    return acos(clamp(x, -1.0f, 1.0f));
}


inline float4 quat_normalize_safe(float4 quat){
    float n2 = dot(quat, quat);
    if (n2 < (EPS_F * EPS_F)) return (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    return fast_normalize(quat);
}


inline float4 quatMultiply(float4 a, float4 b){
    return (float4)(
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z
    );
}


inline Mat4 identity(){
    Mat4 mat;
    mat.c0 = (float4)(1.0f, 0.0f, 0.0f, 0.0f);
    mat.c1 = (float4)(0.0f, 1.0f, 0.0f, 0.0f);
    mat.c2 = (float4)(0.0f, 0.0f, 1.0f, 0.0f);
    mat.c3 = (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    return mat;
}


inline float4 mat4MulVec4(const Mat4 mat, float4 vec){
    return mat.c0 * vec.x + mat.c1 * vec.y + mat.c2 * vec.z + mat.c3 * vec.w;
}


inline Mat4 mat4MulMat4(const Mat4 a, const Mat4 b){
    Mat4 res;
    res.c0 = a.c0 * b.c0.x + a.c1 * b.c0.y + a.c2 * b.c0.z + a.c3 * b.c0.w;
    res.c1 = a.c0 * b.c1.x + a.c1 * b.c1.y + a.c2 * b.c1.z + a.c3 * b.c1.w;
    res.c2 = a.c0 * b.c2.x + a.c1 * b.c2.y + a.c2 * b.c2.z + a.c3 * b.c2.w;
    res.c3 = a.c0 * b.c3.x + a.c1 * b.c3.y + a.c2 * b.c3.z + a.c3 * b.c3.w;
    return res;
}


inline void setTranslationToMatrix(__private Mat4* mat, float4 translate){
    mat->c3.x = translate.x;
    mat->c3.y = translate.y;
    mat->c3.z = translate.z;
}


inline void setQuatToMatrix(__private Mat4* mat, float4 quat){
    quat = quat_normalize_safe(quat);
    float xn = quat.x, yn = quat.y, zn = quat.z, wn = quat.w;

    float xx = xn * xn, yy = yn * yn, zz = zn * zn;
    float xy = xn * yn, xz = xn * zn, yz = yn * zn;
    float wx = wn * xn, wy = wn * yn, wz = wn * zn;

    float m00 = 1.f - 2.f * (yy + zz);
    float m01 = 2.f * (xy - wz);
    float m02 = 2.f * (xz + wy);

    float m10 = 2.f * (xy + wz);
    float m11 = 1.f - 2.f * (xx + zz);
    float m12 = 2.f * (yz - wx);

    float m20 = 2.f * (xz - wy);
    float m21 = 2.f * (yz + wx);
    float m22 = 1.f - 2.f * (xx + yy);

    mat->c0 = (float4)(m00, m10, m20, 0.f);
    mat->c1 = (float4)(m01, m11, m21, 0.f);
    mat->c2 = (float4)(m02, m12, m22, 0.f);
    mat->c3.w = 1.f;
}


inline void setScaleToMatrix(__private Mat4* mat, float4 scale){
    mat->c0.x = scale.x;
    mat->c1.y = scale.y;
    mat->c2.z = scale.z;
}


inline void setShearToMatrix(__private Mat4* mat, float4 shear){
    float shXY = shear.x, shXZ = shear.y, shYZ = shear.z;
    mat->c1.x = shXY;
    mat->c2.x = shXZ;
    mat->c2.y = shYZ;
}


// compose: T * R * S * Sh
inline Mat4 composeMatrix(float4 translate,
                          float4 quaternion,
                          float4 scale,
                          float4 shear)
{
    Mat4 TR = identity();
    Mat4 S = identity();
    Mat4 Sh = identity();

    setTranslationToMatrix(&TR, translate);
    setQuatToMatrix(&TR, quaternion);
    setScaleToMatrix(&S, scale);
    setShearToMatrix(&Sh, shear);

    return mat4MulMat4(mat4MulMat4(TR, S), Sh);
}


inline float4 slerp(float4 a, float4 b, float t){
    a = quat_normalize_safe(a);
    b = quat_normalize_safe(b);

    float cos_omega = dot(a, b);
    if (cos_omega < 0.0f) {
        cos_omega = -cos_omega;
        b = -b;
    }

    float theta = safe_acosf(cos_omega);
    if (fabs(theta) < EPS_F) {
        return quat_normalize_safe(mix(a, b, t));
    }

    float sTheta = native_sin(theta);
    if (fabs(sTheta) < EPS_F) {
        return quat_normalize_safe(mix(a, b, t));
    }

    float s0 = native_sin((1.0f - t) * theta) / sTheta;
    float s1 = native_sin(t * theta) / sTheta;

    return quat_normalize_safe(a * s0 + b * s1);
}



/*
    rpq9MultiSoftModDeformer kernels
*/

__kernel void noneSoftMod(
    __global const float*   initPositions,
    __global const float*   localEnvelopes,
    __global const float4*  deformCenterPositions,
    __global const float4*  translations,
    __global const float4*  quaternions,
    __global const float4*  scales,
    __global const float4*  shears,
    __global const float*   falloffRadiuses,
    __global const float*   localWeights,
    __global const float*   vertWeights,
    __global const uint*    affectMap,
    __global float*         deformedPositions,
    const float             deformerEnvelope,
    const uint              numVertices,
    const uint              numControlMatrices)
{
    const uint gid = get_global_id(0);
    if (gid >= numVertices) return;

    const uint vid = (affectMap ? affectMap[gid] : gid);

    float3 pt = vload3(vid, initPositions);
    float4 clPt = (float4)(pt.x, pt.y, pt.z, 1.0f);

    float4 addPt = (float4)(0.0f, 0.0f, 0.0f, 0.0f);

    const float4 zeroVec = (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    const float4 oneVec  = (float4)(1.0f, 1.0f, 1.0f, 1.0f);

    for (uint j = 0; j < numControlMatrices; ++j) {
        float4 center = deformCenterPositions[j];
        float distanceFromCenter = fast_distance(clPt, center);
        float radius = falloffRadiuses[j];

        if (distanceFromCenter > radius) continue;

        Mat4 transMat = composeMatrix(translations[j], quaternions[j], scales[j], shears[j]);
        float4 localPt = clPt - center;
        localPt.w = 1.0f;
        float4 newPtBase = mat4MulVec4(transMat, localPt);
        addPt += (newPtBase - localPt) * localEnvelopes[j] * localWeights[j * numVertices + vid];
    }

    float weight = (vertWeights ? vertWeights[vid] * deformerEnvelope : deformerEnvelope);
    float3 outPt = pt + weight * (float3)(addPt.x, addPt.y, addPt.z);

    vstore3(outPt, vid, deformedPositions);
}


__kernel void linearSoftMod(
    __global const float*   initPositions,
    __global const float*   localEnvelopes,
    __global const float4*  deformCenterPositions,
    __global const float4*  translations,
    __global const float4*  quaternions,
    __global const float4*  scales,
    __global const float4*  shears,
    __global const float*   falloffRadiuses,
    __global const float* localWeights,
    __global const float*   vertWeights,
    __global const uint*    affectMap,
    __global float*         deformedPositions,
    const float             deformerEnvelope,
    const uint              numVertices,
    const uint              numControlMatrices)
{
    const uint gid = get_global_id(0);
    if (gid >= numVertices) return;

    const uint vid = (affectMap ? affectMap[gid] : gid);

    float3 pt = vload3(vid, initPositions);
    float4 clPt = (float4)(pt.x, pt.y, pt.z, 1.0f);

    float4 addPt = (float4)(0.0f, 0.0f, 0.0f, 0.0f);

    const float4 zeroVec = (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    const float4 oneVec  = (float4)(1.0f, 1.0f, 1.0f, 1.0f);

    for (uint j = 0; j < numControlMatrices; ++j) {
        float4 center = deformCenterPositions[j];
        float distanceFromCenter = fast_distance(clPt, center);
        float radius = falloffRadiuses[j];

        if (distanceFromCenter > radius) continue;

        float falloffWeight = remap1to0(distanceFromCenter, radius);

        float4 weightedTranslate  = mix(zeroVec, translations[j], falloffWeight);
        float4 weightedScale      = mix(oneVec,  scales[j], falloffWeight);
        float4 weightedShear      = mix(zeroVec, shears[j], falloffWeight);
        float4 weightedQuaternion = slerp(zeroVec, quaternions[j], falloffWeight);

        Mat4 transMat = composeMatrix(weightedTranslate, weightedQuaternion, weightedScale, weightedShear);
        float4 localPt = clPt - center;
        localPt.w = 1.0f;
        float4 newPtBase = mat4MulVec4(transMat, localPt);
        addPt += (newPtBase - localPt) * localEnvelopes[j] * localWeights[j * numVertices + vid];
    }

    float weight = (vertWeights ? vertWeights[vid] * deformerEnvelope : deformerEnvelope);
    float3 outPt = pt + weight * (float3)(addPt.x, addPt.y, addPt.z);

    vstore3(outPt, vid, deformedPositions);
}


__kernel void smoothSoftMod(
    __global const float*   initPositions,
    __global const float*   localEnvelopes,
    __global const float4*  deformCenterPositions,
    __global const float4*  translations,
    __global const float4*  quaternions,
    __global const float4*  scales,
    __global const float4*  shears,
    __global const float*   falloffRadiuses,
    __global const float* localWeights,
    __global const float*   vertWeights,
    __global const uint*    affectMap,
    __global float*         deformedPositions,
    const float             deformerEnvelope,
    const uint              numVertices,
    const uint              numControlMatrices)
{
    const uint gid = get_global_id(0);
    if (gid >= numVertices) return;

    const uint vid = (affectMap ? affectMap[gid] : gid);

    float3 pt = vload3(vid, initPositions);
    float4 clPt = (float4)(pt.x, pt.y, pt.z, 1.0f);

    float4 addPt = (float4)(0.0f, 0.0f, 0.0f, 0.0f);

    const float4 zeroVec = (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    const float4 oneVec  = (float4)(1.0f, 1.0f, 1.0f, 1.0f);

    for (uint j = 0; j < numControlMatrices; ++j) {
        float4 center = deformCenterPositions[j];
        float distanceFromCenter = fast_distance(clPt, center);
        float radius = falloffRadiuses[j];

        if (distanceFromCenter > radius) continue;

        float falloffWeight = smoothQuadratic(remap1to0(distanceFromCenter, radius));

        float4 weightedTranslate  = mix(zeroVec, translations[j], falloffWeight);
        float4 weightedScale      = mix(oneVec,  scales[j], falloffWeight);
        float4 weightedShear      = mix(zeroVec, shears[j], falloffWeight);
        float4 weightedQuaternion = slerp(zeroVec, quaternions[j], falloffWeight);

        Mat4 transMat = composeMatrix(weightedTranslate, weightedQuaternion, weightedScale, weightedShear);
        float4 localPt = clPt - center;
        localPt.w = 1.0f;
        float4 newPtBase = mat4MulVec4(transMat, localPt);
        addPt += (newPtBase - localPt) * localEnvelopes[j] * localWeights[j * numVertices + vid];
    }

    float weight = (vertWeights ? vertWeights[vid] * deformerEnvelope : deformerEnvelope);
    float3 outPt = pt + weight * (float3)(addPt.x, addPt.y, addPt.z);

    vstore3(outPt, vid, deformedPositions);
}



__kernel void easeInOutSoftMod(
    __global const float*   initPositions,
    __global const float*   localEnvelopes,
    __global const float4*  deformCenterPositions,
    __global const float4*  translations,
    __global const float4*  quaternions,
    __global const float4*  scales,
    __global const float4*  shears,
    __global const float*   falloffRadiuses,
    __global const float* localWeights,
    __global const float*   vertWeights,
    __global const uint*    affectMap,
    __global float*         deformedPositions,
    const float             deformerEnvelope,
    const uint              numVertices,
    const uint              numControlMatrices)
{
    const uint gid = get_global_id(0);
    if (gid >= numVertices) return;

    const uint vid = (affectMap ? affectMap[gid] : gid);

    float3 pt = vload3(vid, initPositions);
    float4 clPt = (float4)(pt.x, pt.y, pt.z, 1.0f);

    float4 addPt = (float4)(0.0f, 0.0f, 0.0f, 0.0f);

    const float4 zeroVec = (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    const float4 oneVec  = (float4)(1.0f, 1.0f, 1.0f, 1.0f);

    for (uint j = 0; j < numControlMatrices; ++j) {
        float4 center = deformCenterPositions[j];
        float distanceFromCenter = fast_distance(clPt, center);
        float radius = falloffRadiuses[j];

        if (distanceFromCenter > radius) continue;

        float falloffWeight = easeInOutQuadratic(remap1to0(distanceFromCenter, radius));

        float4 weightedTranslate  = mix(zeroVec, translations[j], falloffWeight);
        float4 weightedScale      = mix(oneVec,  scales[j], falloffWeight);
        float4 weightedShear      = mix(zeroVec, shears[j], falloffWeight);
        float4 weightedQuaternion = slerp(zeroVec, quaternions[j], falloffWeight);

        Mat4 transMat = composeMatrix(weightedTranslate, weightedQuaternion, weightedScale, weightedShear);
        float4 localPt = clPt - center;
        localPt.w = 1.0f;
        float4 newPtBase = mat4MulVec4(transMat, localPt);
        addPt += (newPtBase - localPt) * localEnvelopes[j] * localWeights[j * numVertices + vid];
    }

    float weight = (vertWeights ? vertWeights[vid] * deformerEnvelope : deformerEnvelope);
    float3 outPt = pt + weight * (float3)(addPt.x, addPt.y, addPt.z);

    vstore3(outPt, vid, deformedPositions);
}