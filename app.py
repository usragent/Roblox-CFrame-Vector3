import math
import numpy as np

"""


📚 [Vector3 and CFrame Classes in Python] 📚

This module contains two classes: `Vector3` and `CFrame`.

    🚀 Vector3 Class:

        The Vector3 class represents a 3D vector in a Cartesian coordinate system. It is independent
        and can run on its own. It supports basic vector operations such as addition, subtraction, 
        multiplication, and dot product. It also provides methods to calculate the magnitude and unit
        vector, and perform cross product operation.

    ----------------------------------------------------------------------------------

    🧩 CFrame Class:

        The CFrame class represents a coordinate frame, which is a combination of a a Vector3
        and a rotation (a 3x3 matrix). It requires the Vector3 class to function properly.
        It supports multiplication with another CFrame or a Vector3, and provides methods to
        convert between world space and object space, calculate the inverse, and perform operations
        related to Euler angles and axis-angle rotation.

🏆 Credits:

    - (x86) for the creation of Vector3 and CFrame


"""

class Vector3:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"{int(self.x) if self.x % 1 == 0.0 else self.x}, {int(self.y) if self.y % 1 == 0.0 else self.y}, {int(self.z) if self.z % 1 == 0.0 else self.z}"

    def __mul__(self, other):
        if isinstance(other, CFrame):
            return other * self
        return Vector3(self.x * other, self.y * other, self.z * other)

    def __add__(self, other):
        return Vector3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Vector3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __neg__(self):
        return Vector3(-self.x, -self.y, -self.z)

    @property
    def magnitude(self):
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    @property
    def unit(self):
        mag = self.magnitude
        return Vector3(self.x / mag, self.y / mag, self.z / mag)

    def cross(self, b):
        return Vector3(self.y * b.z - self.z * b.y, self.z * b.x - self.x * b.z, self.x * b.y - self.y * b.x)

    def dot(self, b):
        return self.x * b.x + self.y * b.y + self.z * b.z


class CFrame:
    def __init__(self, *args):
            if len(args) == 1 and isinstance(args[0], Vector3):
                self.position = args[0]
                self.rotation = np.identity(3)
            elif len(args) == 3 and all(isinstance(arg, float) for arg in args):
                self.position = Vector3(*args)
                self.rotation = np.identity(3)
            elif len(args) == 12:
                self.position = Vector3(args[0], args[1], args[2])
                self.rotation = np.array(args[3:]).reshape((3, 3))
            else:
                raise ValueError("[ERROR] Invalid args for CFrame constructor...")
                
    def __repr__(self):
        rotation = ', '.join(str(int(i) if i % 1 == 0.0 else i) for i in self.rotation.flatten())
        return f"{int(self.position.x) if self.position.x % 1 == 0.0 else self.position.x}, {int(self.position.y) if self.position.y % 1 == 0.0 else self.position.y}, {int(self.position.z) if self.position.z % 1 == 0.0 else self.position.z}, {rotation}"
    
    def __mul__(self, other):
        if isinstance(other, Vector3):
            rotPos = np.dot(self.rotation, np.array([other.x, other.y, other.z]))
            return Vector3(rotPos[0], rotPos[1], rotPos[2]) + self.position
        elif isinstance(other, CFrame):
            newPos = self * other.position
            newRotation = np.dot(self.rotation, other.rotation)
            return CFrame(newPos.x, newPos.y, newPos.z, *newRotation.flatten())
        else:
            raise ValueError("[ERROR] CFrame can only be multiplied by Vector3 or CFrame...")
    
    def toWorldSpace(self, cf2):
        return self * cf2

    def toObjectSpace(self, cf2):
        invSelf = self.inverse()
        return invSelf * cf2

    def inverse(self):
        invRotation = np.linalg.inv(self.rotation)
        invPosition = np.dot(-invRotation, np.array([self.position.x, self.position.y, self.position.z]))
        return CFrame(invPosition[0], invPosition[1], invPosition[2], *invRotation.flatten())

    def pointToWorldSpace(self, v):
        return self * v

    def pointToObjectSpace(self, v):
        return self.inverse() * v
        
    @staticmethod
    def fromEulerAnglesXYZ(x, y, z):
        cfx = CFrame.fromAxisAngle(Vector3(1, 0, 0), x)
        cfy = CFrame.fromAxisAngle(Vector3(0, 1, 0), y)
        cfz = CFrame.fromAxisAngle(Vector3(0, 0, 1), z)
        return cfx * cfy * cfz

    @staticmethod
    def fromEulerAnglesYXZ(y, x, z):
        cfx = CFrame.fromAxisAngle(Vector3(1, 0, 0), x)
        cfy = CFrame.fromAxisAngle(Vector3(0, 1, 0), y)
        cfz = CFrame.fromAxisAngle(Vector3(0, 0, 1), z)
        return cfy * cfx * cfz

    def toEulerAnglesXYZ(self):
        sy = math.sqrt(self.rotation[0,0] * self.rotation[0,0] +  self.rotation[1,0] * self.rotation[1,0])
        singular = sy < 1e-6

        if not singular:
            x = math.atan2(self.rotation[2,1] , self.rotation[2,2])
            y = math.atan2(-self.rotation[2,0], sy)
            z = math.atan2(self.rotation[1,0], self.rotation[0,0])
        else:
            x = math.atan2(-self.rotation[1,2], self.rotation[1,1])
            y = math.atan2(-self.rotation[2,0], sy)
            z = 0
        return np.array([x, y, z])
    
    @staticmethod
    def Angles(x, y, z):
        return CFrame.fromEulerAnglesXYZ(x, y, z)

    @staticmethod
    def fromAxisAngle(axis, theta):
        axis = axis.unit
        cosTheta = math.cos(theta)
        sinTheta = math.sin(theta)
        oneMinusCosTheta = 1 - cosTheta

        m11 = cosTheta + axis.x**2 * oneMinusCosTheta
        m12 = axis.x * axis.y * oneMinusCosTheta - axis.z * sinTheta
        m13 = axis.x * axis.z * oneMinusCosTheta + axis.y * sinTheta

        m21 = axis.y * axis.x * oneMinusCosTheta + axis.z * sinTheta
        m22 = cosTheta + axis.y**2 * oneMinusCosTheta
        m23 = axis.y * axis.z * oneMinusCosTheta - axis.x * sinTheta

        m31 = axis.z * axis.x * oneMinusCosTheta - axis.y * sinTheta
        m32 = axis.z * axis.y * oneMinusCosTheta + axis.x * sinTheta
        m33 = cosTheta + axis.z**2 * oneMinusCosTheta

        return CFrame(0, 0, 0, m11, m12, m13, m21, m22, m23, m31, m32, m33)

    @staticmethod
    def lookAt(at, look):
        zAxis = (at - look).unit
        xAxis = Vector3.cross(Vector3(0, 1, 0), zAxis).unit
        yAxis = Vector3.cross(zAxis, xAxis).unit

        if xAxis.magnitude == 0:
            xAxis = Vector3(0, 0, zAxis.y).unit
            yAxis = Vector3(zAxis.y, 0, 0).unit
            zAxis = Vector3(0, -1, 0) if zAxis.y < 0 else Vector3(0, 1, 0)

        return CFrame.fromMatrix(at, xAxis, yAxis, zAxis)
        
    @staticmethod
    def fromMatrix(pos, vX, vY, vZ):
        rotation = np.array([vX.x, vY.x, vZ.x, vX.y, vY.y, vZ.y, vX.z, vY.z, vZ.z]).reshape((3, 3))
        return CFrame(pos.x, pos.y, pos.z, *rotation.flatten())
