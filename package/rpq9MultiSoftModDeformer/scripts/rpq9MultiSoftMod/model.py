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
import dataclasses
from typing import Any, get_args, get_origin, get_type_hints, Union
from enum import IntEnum
from logging import getLogger

logger = getLogger(__name__)


def isValueOfType(value: Any, expectedType: Any) -> bool:
    origin = get_origin(expectedType)
    args = get_args(expectedType)

    # Any
    if expectedType is Any:
        return True

    # Union / Optional
    if origin is Union:
        return any(isValueOfType(value, arg) for arg in args)

    # list[T]
    if origin is list:
        if not isinstance(value, list):
            return False
        if not args:
            return True
        return all(isValueOfType(v, args[0]) for v in value)

    # dict[K, V]
    if origin is dict:
        if not isinstance(value, dict):
            return False
        if len(args) != 2:
            return True
        key_type, val_type = args
        return all(
            isValueOfType(k, key_type) and isValueOfType(v, val_type)
            for k, v in value.items()
        )

    # tuple[T1, T2, ...] / tuple[T, ...]
    if origin is tuple:
        if not isinstance(value, tuple):
            return False
        if not args:
            return True

        # tuple[T, ...]
        if len(args) == 2 and args[1] is Ellipsis:
            return all(isValueOfType(v, args[0]) for v in value)

        # tuple[T1, T2, ...]
        if len(value) != len(args):
            return False
        return all(isValueOfType(v, t) for v, t in zip(value, args))

    # set[T]
    if origin is set:
        if not isinstance(value, set):
            return False
        if not args:
            return True
        return all(isValueOfType(v, args[0]) for v in value)

    # nested dataclass
    if isinstance(expectedType, type) and dataclasses.is_dataclass(expectedType):
        return isinstance(value, expectedType)

    # normal type
    if isinstance(expectedType, type):
        return isinstance(value, expectedType)

    # Types that cannot be determined are set to False here, not True.
    return False


def isDataclassAllFieldsValid(obj: Any, verbose:bool=False) -> bool:
    if not dataclasses.is_dataclass(obj) or isinstance(obj, type):
        raise TypeError(f'obj is not dataclass instance.')

    hints = get_type_hints(type(obj))
    res = True

    for f in dataclasses.fields(obj):
        value = getattr(obj, f.name)
        expectedType = hints.get(f.name, f.type)
        ok = isValueOfType(value, expectedType)

        if not ok:
            if verbose:
                logger.info(f'"{f.name}" expected type is "{expectedType}". but actual type is {type(value)}')
                res = False
            else:
                return False

    return res


class FalloffMode(IntEnum):
    kNone = 0,
    kLinear = 1
    kSmooth = 2
    kSpline = 3
    kSmoothStep = 4
    kEaseInOut = 5


@dataclasses.dataclass
class Matrix4x4:
    m11: float
    m12: float
    m13: float
    m14: float
    m21: float
    m22: float
    m23: float
    m24: float
    m31: float
    m32: float
    m33: float
    m34: float
    m41: float
    m42: float
    m43: float
    m44: float

    def __getitem__(self, index: tuple[int, int]) -> float:
        row, col = index
        if not (0 <= row < 4 and 0 <= col < 4):
            raise IndexError("row and col must be in range 0..3")

        names = (
            ("m11", "m12", "m13", "m14"),
            ("m21", "m22", "m23", "m24"),
            ("m31", "m32", "m33", "m34"),
            ("m41", "m42", "m43", "m44"),
        )
        return getattr(self, names[row][col])

    def __setitem__(self, index: tuple[int, int], value: float) -> None:
        row, col = index
        if not (0 <= row < 4 and 0 <= col < 4):
            raise IndexError("row and col must be in range 0..3")

        names = (
            ("m11", "m12", "m13", "m14"),
            ("m21", "m22", "m23", "m24"),
            ("m31", "m32", "m33", "m34"),
            ("m41", "m42", "m43", "m44"),
        )
        setattr(self, names[row][col], value)

    def toList(self) -> list[float]:
        return [
            self.m11, self.m12, self.m13, self.m14,
            self.m21, self.m22, self.m23, self.m24,
            self.m31, self.m32, self.m33, self.m34,
            self.m41, self.m42, self.m43, self.m44,
        ]


class DeformerWeights(dict):
    def __init__(self, *args, **kwargs):
        super().__init__()
        temp = dict(*args, **kwargs)
        for key, value in temp.items():
            self[key] = value

    def __setitem__(self, key:str, value:float):
        if type(key) is not str:
            raise TypeError('The key type must be str.') 

        if type(value) is not float:
            raise TypeError('The value type must be float.')

        super().__setitem__(key, value)



class DeformerWeightList(dict):
    def __init__(self, *args, **kwargs):
        super().__init__()
        temp = dict(*args, **kwargs)
        for key, value in temp.items():
            self[key] = value

    def __setitem__(self, key:str, value:DeformerWeights):
        if type(key) is not str:
            raise TypeError('The key type must be str.') 

        if type(value) is not DeformerWeight:
            raise TypeError('The value type must be DeformerWeight.')

        super().__setitem__(key, value)



@dataclasses.dataclass
class Rpq9MultiSoftModData:
    name: str
    geometries: list[str]
    geometryIndices: list[int]
    envelope: float
    falloffMode: int
    localEnvelopes: dict[str, float]
    centerMatrices: dict[str, Matrix4x4]
    modifyMatrices: dict[str, Matrix4x4]
    falloffRadii: dict[str, float]
    weightList: DeformerWeightList
    localWeightList: dict[str, DeformerWeightList]


    def toDict(self) -> dict:
        return dataclasses.asdict(self)


    def fromJsonDict(cls, data:dict) -> 'Rpq9MultiSoftModData':
        kwargs = {
            'name': data['name'],
            'geometries': data['geometries'],
            'geometryIndices': data['geometryIndices'],
            'envelope': data['envelope'],
            'falloffMode': data['falloffMode'],
            'localEnvelopes': data['localEnvelopes'],
            'centerMatrices': {key: Matrix4x4(**value) for key, value in data['centerMatrices'].items()},
            'modifyMatrices': {key: Matrix4x4(**value) for key, value in data['modifyMatrices'].items()},
            'falloffRadii': data['falloffRadii'],
            'weightList': DeformerWeightList(**data['weightList']),
            'localWeightList': {key: DeformerWeightList({index: DeformerWeights(**value) for index, value in weightListValue.items()}) for key, weightListValue in data['localWeightList'].items()},
        }

        return cls(**kwargs)


    @staticmethod
    def haveSameKeys(*dicts) -> bool:
        if not dicts:
            return True
        first_keys = set(dicts[0].keys())
        return all(set(d.keys()) == first_keys for d in dicts[1:])


    def isValid(self, verbose:bool=False) -> bool:
        res = True
        if len(self.geometries) != len(self.geometryIndices):
            if verbose:
                logger.info('The number of elements in geometries and geometryIndices do not match.')
                res = False
            else:
                return False

        if not(len(self.localEnvelopes) == len(self.centerMatrices) == len(self.modifyMatrices) == len(self.falloffRadii), len(self.localWeightList)):
            if verbose:
                logger.info('The number of elements in the input data is different.')
                res = False
            else:
                return False

        if not Rpq9MultiSoftModData.haveSameKeys(self.localEnvelopes, self.centerMatrices, self.modifyMatrices, self.falloffRadii, self.localWeightList):
            if verbose:
                logger.info('The element indices in inputData are different.')
                res = False
            else:
                return False

        geometryIndicesStr = {str(index) for index in self.geometryIndices}
        for index in self.weightList.keys():
            if index not in geometryIndicesStr:
                if verbose:
                    logger.info(f'Using a non-existent geometry index in weightList. (geoIndex: {index})')
                    res = False
                else:
                    return False

        for matrixIndex, geoWeightData in localWeightList.items():
            for geoIndex in geoWeightData.key():
                if geoIndex not in geometryIndicesStr:
                    if verbose:
                        logger.info(f'Using a non-existent geometry index in localWeightList. (geoIndex: {geoIndex}, inputDataIndex: {matrixIndex})')
                        res = False
                    else:
                        return False

        if not isDataclassAllFieldsValid(self, verbose):
            if verbose:
                logger.info('The field data type is incorrect.')
                res = False

        return res


    def reindexingGeometry(self) -> None:
        geoNum = len(self.geometryIndices)

        new_data = {
            'geometryIndices': [i for i in range(geoNum)],
            'weightList': {},
            'localWeightList': {}
        }

        shiftedIndices = {}
        shiftValueCache = 0
        for i, geoIndex in enumerate(self.geometryIndices):
            geoIndexint = int(geoIndex)
            shiftedIndex = geoIndexint + shiftValueCache
            if shiftedIndex == i:
                shiftedIndices[geoIndex] = str(shiftedIndex)
            else:
                shiftValueCache -= shiftedIndex-i
                shiftedIndices[geoIndex] = str(geoIndexint + shiftValueCache)

        for geoIndex, weights in self.weightList.items():
            new_data['weightList'][shiftedIndices[geoIndex]] = weights

        for matrixIndex, geoWeightData in self.localWeightList.items():
            new_data['localWeightList'][matrixIndex] = {}
            for geoIndex, weights in geoWeightData.items():
                new_data['localWeightList'][matrixIndex][shiftedIndices[geoIndex]] = weights


        self.geometryIndices = new_data['geometryIndices']
        self.weightList = new_data['weightList']
        self.localWeightList = new_data['localWeightList']


    def removeGeometry(self, removeGeo:str) -> int:
        geoIndex = self.geometryIndices[self.geometries.index(removeGeo)]

        del self.weightList[str(geoIndex)]
        for matrixIndex, geoWeightData in self.localWeightList.items():
            del geoWeightData[str(geoIndex)]

        return geoIndex