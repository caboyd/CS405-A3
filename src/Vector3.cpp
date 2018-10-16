#include "Vector3.h"
#include <cassert>


const Vector3 Vector3::ZERO(0.0, 0.0, 0.0);


Vector3::Vector3(): x(0), y(0), z(0)
{
}

Vector3::Vector3(float value):x(value),y(value),z(value)
{
}

Vector3::Vector3(float x, float y, float z): x(x), y(y), z(z)
{
}

Vector3::Vector3(const Vector3& other) = default;
Vector3& Vector3::operator=(const Vector3& other) = default;

bool Vector3::operator==(const Vector3& other) const
{
	if (x != other.x) return false;
	if (y != other.y) return false;
	if (z != other.z) return false;
	return true;
}

void Vector3::set(float x, float y, float z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

void Vector3::normalize()
{
	if (isZero()) return;

	const float inv_length = 1.0f / this->length();
	x *= inv_length;
	y *= inv_length;
	z *= inv_length;
}

void Vector3::negate()
{
	x = -x;
	y = -y;
	z = -z;
}

bool Vector3::isZero() const
{
	if (x != 0.0) return false;;
	if (y != 0.0) return false;
	if (z != 0.0) return false;
	return true;
}

float Vector3::length() const
{
	return sqrtf(x * x + y * y + z * z);
}

float Vector3::getSquaredLength() const
{
	return (x * x + y * y + z * z);
}

Vector3 Vector3::getNormalized() const
{
	if (isZero())
		return Vector3();

	const float inv_length = 1.0f / length();
	return {x * inv_length, y * inv_length, z * inv_length};
}

float Vector3::distance(const Vector3& v) const
{
	return (*this - v).length();
}

float Vector3::distanceSquared(const Vector3& v) const
{
	return (*this - v).getSquaredLength();
}

float Vector3::getLargestComponentIndex() const
{
	float a_x = abs(x);
	float a_y = abs(y);
	float a_z = abs(z);

	if(a_x >= a_y && a_x >= a_z) return 0;
	if(a_y >= a_x && a_y >= a_z) return 1;
	if(a_z >= a_y && a_z >= a_x) return 2;
	assert(false);
	return 0;
}

Vector3 Vector3::operator-() const
{
	return {-x, -y, -z};
}

Vector3 Vector3::operator+(const Vector3& v) const
{
	return {this->x + v.x, this->y + v.y, this->z + v.z};
}

Vector3 Vector3::operator-(const Vector3& v) const
{
	return {this->x - v.x, this->y - v.y, this->z - v.z};
}

Vector3 Vector3::operator*(const Vector3& v) const
{
	return {this->x * v.x, this->y * v.y, this->z * v.z};
}

Vector3 Vector3::operator/(const Vector3& v) const
{
	return {this->x / v.x, this->y / v.y, this->z / v.z};
}

Vector3& Vector3::operator+=(const Vector3& v)
{
	this->x += v.x;
	this->y += v.y;
	this->z += v.z;
	return *this;
}

Vector3& Vector3::operator-=(const Vector3& v)
{
	this->x -= v.x;
	this->y -= v.y;
	this->z -= v.z;
	return *this;
}

Vector3& Vector3::operator*=(const Vector3& v)
{
	this->x *= v.x;
	this->y *= v.y;
	this->z *= v.z;
	return *this;
}

Vector3& Vector3::operator/=(const Vector3& v)
{
	this->x /= v.x;
	this->y /= v.y;
	this->z /= v.z;
	return *this;
}

Vector3 Vector3::operator-(float f) const
{
	return {x - f, y - f, z - f};
}

Vector3 Vector3::operator+(float f) const
{
	return {x + f, y + f, z + f};
}

Vector3 Vector3::operator/(float scalar) const
{
	return {x / scalar, y / scalar, z / scalar};
}

Vector3 Vector3::operator*(float scalar) const
{
	return {x * scalar, y * scalar, z * scalar};
}

Vector3& Vector3::operator+=(float f)
{
	this->x += f;
	this->y += f;
	this->z += f;
	return *this;
}

Vector3& Vector3::operator-=(float f)
{
	this->x -= f;
	this->y -= f;
	this->z -= f;
	return *this;
}

Vector3& Vector3::operator*=(float f)
{
	this->x *= f;
	this->y *= f;
	this->z *= f;
	return *this;
}

Vector3& Vector3::operator/=(float f)
{
	this->x /= f;
	this->y /= f;
	this->z /= f;
	return *this;
}

float Vector3::operator[](int i) const
{
	return data[i];
}

float& Vector3::operator[](int i)
{
	return data[i];
}

Vector3 Vector3::cross(const Vector3& v) const
{
	return Vector3
	(
		y * v.z - z * v.y,
		z * v.x - x * v.z,
		x * v.y - y * v.x
	);
}

float Vector3::dot(const Vector3& v) const
{
	return x * v.x + y * v.y + z * v.z;
}

std::ostream& operator<<(std::ostream& os, const Vector3& v)
{
	if (!os.good())
		return os;
	os  << std::fixed;
	os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	os.unsetf(std::ios_base::fixed);
	return os;
}
