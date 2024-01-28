# Roblox CFrame & Vector3 Functions

This repository contains Python implementations of Roblox's Vector3 and CFrame classes. These classes are used for 3D mathematical operations such as vector addition, subtraction, multiplication, and dot product, as well as operations related to Euler angles, axis-angle rotation and point to screen space.

## Classes

### Vector3 Class

The Vector3 class represents a 3D vector in a Cartesian coordinate system. It supports basic vector operations such as addition, subtraction, multiplication, and dot product. It also provides methods to calculate the magnitude and unit vector, and perform cross product operation.

### CFrame Class

The CFrame class represents a coordinate frame, which is a combination of a Vector3 and a rotation (a 3x3 matrix). It supports multiplication with another CFrame or a Vector3, and provides methods to convert between world space and object space, calculate the inverse, and perform operations related to Euler angles and axis-angle rotation.

## Usage

You can use these classes in your Python projects that require 3D mathematical operations. Simply import the classes from the module and create instances of them as needed.

# Vector3 usage
```python
v1 = Vector3(1, 2, 3)
v2 = Vector3(4, 5, 6)

dotResult = v1.dot(v2)
crossResult = v1.cross(v2)
magnitude = v1.magnitude
unitVector = v1.unit
```

# CFrame usage
```python
cf1 = CFrame(v1)
cf2 = CFrame(v2)

toWorldSpace = cf1.toWorldSpace(cf2)
toObjectSpace = cf1.toObjectSpace(cf2)
inverse = cf1.inverse()
pointToWorldSpace = cf1.pointToWorldSpace(v1)
pointToObjectSpace = cf1.pointToObjectSpace(v1)
fromEulerAngles = CFrame.fromEulerAnglesXYZ(1, 2, 3)
toEulerAnglesXYZ = cf1.toEulerAnglesXYZ()
fromAxisAngle = CFrame.fromAxisAngle(v1, 1)
lookAt = CFrame.lookAt(v1, v2)
fromMatrix = CFrame.fromMatrix(v1, v2, v1, v2)
```

## Credits

This project was created by (x86).

## License

This project is licensed under the MIT License.
