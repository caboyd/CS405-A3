#include "Mat4.h"
#include "Vector3.h"

#include <iomanip>
#include <algorithm>


Mat4::Mat4()
{
	memset(m, 0, sizeof(m));
}

Mat4::Mat4(float value)
{
	for (float & i : m)
		i = value;
}

Mat4::Mat4(float m11, float m12, float m13, float m14,
           float m21, float m22, float m23, float m24,
           float m31, float m32, float m33, float m34,
           float m41, float m42, float m43, float m44)
{
	m[0] = m11;
	m[1] = m12;
	m[2] = m13;
	m[3] = m14;
	m[4] = m21;
	m[5] = m22;
	m[6] = m23;
	m[7] = m24;
	m[8] = m31;
	m[9] = m32;
	m[10] = m33;
	m[11] = m34;
	m[12] = m41;
	m[13] = m42;
	m[14] = m43;
	m[15] = m44;
}

Mat4::Mat4(const Mat4& other)
{
	for (int i = 0; i < 16; i++)
		m[i] = other.m[i];
}


Mat4& Mat4::operator=(const Mat4& other)
{
	for (int i = 0; i < 16; i++)
		m[i] = other.m[i];
	return *this;
}

bool Mat4::operator==(const Mat4& other) const
{
	for (int i = 0; i < 16; i++)
		if (m[i] != other.m[i]) return false;

	return true;
}

bool Mat4::operator!=(const Mat4& other) const
{
	return !(*this == other);
}

void Mat4::set(float m11, float m12, float m13, float m14,
               float m21, float m22, float m23, float m24,
               float m31, float m32, float m33, float m34,
               float m41, float m42, float m43, float m44)
{
	m[0] = m11;
	m[1] = m12;
	m[2] = m13;
	m[3] = m14;
	m[4] = m21;
	m[5] = m22;
	m[6] = m23;
	m[7] = m24;
	m[8] = m31;
	m[9] = m32;
	m[10] = m33;
	m[11] = m34;
	m[12] = m41;
	m[13] = m42;
	m[14] = m43;
	m[15] = m44;
}

void Mat4::transpose()
{
	std::swap(m[1], m[4]);
	std::swap(m[2], m[8]);
	std::swap(m[3], m[12]);
	std::swap(m[6], m[9]);
	std::swap(m[7], m[13]);
	std::swap(m[11], m[14]);
}

void Mat4::inverse()
{
	const float a00 = m[0], a01 = m[1], a02 = m[2], a03 = m[3];
	const float a10 = m[4], a11 = m[5], a12 = m[6], a13 = m[7];
	const float a20 = m[8], a21 = m[9], a22 = m[10], a23 = m[11];
	const float a30 = m[12], a31 = m[13], a32 = m[14], a33 = m[15];

	const float b00 = a00 * a11 - a01 * a10;
	const float b01 = a00 * a12 - a02 * a10;
	const float b02 = a00 * a13 - a03 * a10;
	const float b03 = a01 * a12 - a02 * a11;
	const float b04 = a01 * a13 - a03 * a11;
	const float b05 = a02 * a13 - a03 * a12;
	const float b06 = a20 * a31 - a21 * a30;
	const float b07 = a20 * a32 - a22 * a30;
	const float b08 = a20 * a33 - a23 * a30;
	const float b09 = a21 * a32 - a22 * a31;
	const float b10 = a21 * a33 - a23 * a31;
	const float b11 = a22 * a33 - a23 * a32;

	// Calculate the determinant
	const float det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
	if (det == 0.0) return;

	const float inv_det = 1.0f / det;

	m[0] = (a11 * b11 - a12 * b10 + a13 * b09) * inv_det;
	m[1] = (a02 * b10 - a01 * b11 - a03 * b09) * inv_det;
	m[2] = (a31 * b05 - a32 * b04 + a33 * b03) * inv_det;
	m[3] = (a22 * b04 - a21 * b05 - a23 * b03) * inv_det;
	m[4] = (a12 * b08 - a10 * b11 - a13 * b07) * inv_det;
	m[5] = (a00 * b11 - a02 * b08 + a03 * b07) * inv_det;
	m[6] = (a32 * b02 - a30 * b05 - a33 * b01) * inv_det;
	m[7] = (a20 * b05 - a22 * b02 + a23 * b01) * inv_det;
	m[8] = (a10 * b10 - a11 * b08 + a13 * b06) * inv_det;
	m[9] = (a01 * b08 - a00 * b10 - a03 * b06) * inv_det;
	m[10] = (a30 * b04 - a31 * b02 + a33 * b00) * inv_det;
	m[11] = (a21 * b02 - a20 * b04 - a23 * b00) * inv_det;
	m[12] = (a11 * b07 - a10 * b09 - a12 * b06) * inv_det;
	m[13] = (a00 * b09 - a01 * b07 + a02 * b06) * inv_det;
	m[14] = (a31 * b01 - a30 * b03 - a32 * b00) * inv_det;
	m[15] = (a20 * b03 - a21 * b01 + a22 * b00) * inv_det;
}

void Mat4::identity()
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (j == i)
				m[i * 4 + j] = 1;
			else
				m[i * 4 + j] = 0;
		}
	}
}

Mat4 Mat4::getTransposed() const
{
	Mat4 m(*this);
	m.transpose();
	return m;
}

bool Mat4::isIdentity() const
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = i; j < 4; j++)
		{
			if (m[i * 4 + j] != static_cast<float>(i == j) ||
				m[j * 4 + i] != static_cast<float>(i == j))
				return false;
		}
	}
	return true;
}

bool Mat4::isZero() const
{
	for (float i : m)
		if (i != 0.0f) return false;
	return true;
}

Mat4 Mat4::operator+(const Mat4& rhs) const
{
	Mat4 mat;
	for (int i = 0; i < 16; i++)
		mat[i] = m[i] + rhs[i];
	return mat;
}

Mat4 Mat4::operator-(const Mat4& rhs) const
{
	Mat4 mat;
	for (int i = 0; i < 16; i++)
		mat[i] = m[i] - rhs[i];
	return mat;
}

Mat4 Mat4::fromIdentity()
{
	return Mat4(1.f, 0.f, 0.f, 0.f,
	            0.f, 1.f, 0.f, 0.f,
	            0.f, 0.f, 1.f, 0.f,
	            0.f, 0.f, 0.f, 1.f);
}

Mat4 Mat4::fromTranslation(const Vector3& v)
{
	return Mat4(1.f, 0.f, 0.f, v.x,
	            0.f, 1.f, 0.f, v.y,
	            0.f, 0.f, 1.f, v.z,
	            0.f, 0.f, 0.f, 1.f);
}

Mat4 Mat4::fromScaling(const Vector3& v)
{
	return Mat4(v.x, 0.f, 0.f, 0.f,
	            0.f, v.y, 0.f, 0.f,
	            0.f, 0.f, v.z, 0.f,
	            0.f, 0.f, 0.f, 1.f);
}

Mat4 Mat4::fromXRotation(float rad)
{
	const float s = sin(rad);
	const float c = cos(rad);

	// @formatter:off
	return Mat4(1.f, 0.f, 0.f, 0.f,
	            0.f,   c,  -s, 0.f,
	            0.f,   s,   c, 0.f,
	            0.f, 0.f, 0.f, 1.f);
}

Mat4 Mat4::fromYRotation(float rad)
{
	const float s = sin(rad);
	const float c = cos(rad);

	// @formatter:off
	return Mat4(  c, 0.f,   s, 0.f,
	            0.f, 1.f, 0.f, 0.f,
	             -s, 0.f,   c, 0.f,
	            0.f, 0.f, 0.f, 1.f);
}

Mat4 Mat4::fromZRotation(float rad)
{
	const float s = sin(rad);
	const float c = cos(rad);

	// @formatter:off
	return Mat4(  c,  -s, 0.f, 0.f,
	              s,   c, 0.f, 0.f,
	            0.f, 0.f, 1.f, 0.f,
	            0.f, 0.f, 0.f, 1.f);
}

Mat4 Mat4::operator*(const Mat4& rhs) const
{
	Mat4 mat;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				mat[i * 4 + j] += m[i * 4 + k] * rhs[k * 4 + j];
			}
		}
	}
	return mat;
}


Mat4 Mat4::operator/(const Mat4& rhs) const
{
	Mat4 mat;
	//TODO: Invert rhs then multiply it
	return mat;
}

Mat4 Mat4::operator*(float scalar) const
{
	Mat4 mat;
	for (int i = 0; i < 16; i++)
		mat[i] = m[i] * scalar;
	return mat;
}

Mat4 Mat4::operator/(float scalar) const
{
	Mat4 mat;
	for (int i = 0; i < 16; i++)
		mat[i] = m[i] / scalar;
	return mat;
}

Mat4& Mat4::operator*=(float scalar)
{
	for (int i = 0; i < 16; i++)
		m[i] *= scalar;
	return *this;
}

Mat4& Mat4::operator/=(float scalar)
{
	for (int i = 0; i < 16; i++)
		m[i] /= scalar;
	return *this;
}

Mat4& Mat4::operator+=(const Mat4& rhs)
{
	for (int i = 0; i < 16; i++)
		m[i] += rhs[i];
	return *this;
}

Mat4& Mat4::operator-=(const Mat4& rhs)
{
	for (int i = 0; i < 16; i++)
		m[i] -= rhs[i];
	return *this;
}

Mat4& Mat4::operator*=(const Mat4& rhs)
{
	Mat4 mat(*this);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			m[i * 4 + j] = 0;
			for (int k = 0; k < 4; k++)
			{
				m[i * 4 + j] += mat[i * 4 + k] * rhs[k * 4 + j];
			}
		}
	}
	return *this;
}

Mat4& Mat4::operator/=(const Mat4& rhs)
{
	Mat4 m(rhs);
	m.inverse();
	return *this *= m;
}

std::ostream& operator<<(std::ostream& os, const Mat4& m)
{

	for(int i = 0; i < 4; i++)
	{
		if( i == 0)os << "["; else os << " ";
		for(int j = 0; j < 4; j++)
		{
			os << std::setw(9) << std::fixed;
			os << m[i*4 + j];
			if(i != 3 || j !=3) os << ", ";
		}
		if (i == 3) os << "]";
		os << std::endl;
	}
	return os;
}

Mat4 Mat4::translate(const Vector3& v) const

{
	Mat4 mat(*this);
	mat[3] += m[0] * v.x + m[1] * v.y + m[2] * v.z;
	mat[7] += m[4] * v.x + m[5] * v.y + m[6] * v.z;
	mat[11] += m[8] * v.x + m[9] * v.y + m[10] * v.z;
	return mat;
}

Mat4& Mat4::translate(const Vector3& v)
{
	m[3] += m[0] * v.x + m[1] * v.y + m[2] * v.z;
	m[7] += m[4] * v.x + m[5] * v.y + m[6] * v.z;
	m[11] += m[8] * v.x + m[9] * v.y + m[10] * v.z;
	return *this;
}

Mat4 Mat4::scale(const Vector3& v) const
{
	return Mat4(m[0] * v.x, m[1] * v.y, m[2] * v.z, m[3],
	            m[4] * v.x, m[5] * v.y, m[6] * v.z, m[7],
	            m[8] * v.x, m[9] * v.y, m[10] * v.z, m[11],
	            m[12] * v.x, m[13] * v.y, m[14] * v.z, m[15]);
}

Mat4& Mat4::scale(const Vector3& v)
{
	set(m[0] * v.x, m[1] * v.y, m[2] * v.z, m[3],
	    m[4] * v.x, m[5] * v.y, m[6] * v.z, m[7],
	    m[8] * v.x, m[9] * v.y, m[10] * v.z, m[11],
	    m[12] * v.x, m[13] * v.y, m[14] * v.z, m[15]);
	return *this;
}

float Mat4::operator[](int i) const
{
	return m[i];
}

float& Mat4::operator[](int i)
{
	return m[i];
}

float Mat4::determinant() const
{
	const float b00 = m[0] * m[5] - m[1] * m[4];
	const float b01 = m[0] * m[6] - m[2] * m[4];
	const float b02 = m[0] * m[7] - m[3] * m[4];
	const float b03 = m[1] * m[6] - m[2] * m[5];
	const float b04 = m[1] * m[7] - m[3] * m[5];
	const float b05 = m[2] * m[7] - m[3] * m[6];
	const float b06 = m[8] * m[13] - m[9] * m[12];
	const float b07 = m[8] * m[14] - m[10] * m[12];
	const float b08 = m[8] * m[15] - m[11] * m[12];
	const float b09 = m[9] * m[14] - m[10] * m[13];
	const float b10 = m[9] * m[15] - m[11] * m[13];
	const float b11 = m[10] * m[15] - m[11] * m[14];

	// Calculate the determinant
	return b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
}
