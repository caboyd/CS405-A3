#pragma once
class Vector3;
#include <ostream>

class Mat4
{
public:
	float m[16]{};

	//Constructors
	Mat4();
	explicit Mat4(float value);
	Mat4(float m11, float m12, float m13, float m14,
	     float m21, float m22, float m23, float m24,
	     float m31, float m32, float m33, float m34,
	     float m41, float m42, float m43, float m44);

	//Copy Constructor
	Mat4(const Mat4& other);

	static Mat4 fromIdentity();
	static Mat4 fromTranslation(const Vector3& v);
	static Mat4 fromScaling(const Vector3& v);

	static Mat4 fromXRotation(float rad);
	static Mat4 fromYRotation(float rad);
	static Mat4 fromZRotation(float rad);

	//Copy Assignment 
	Mat4& operator=(const Mat4& other);

	//Comparison
	bool operator==(const Mat4& other) const;
	bool operator!=(const Mat4& other) const;

	//Mutators
	void set(float m11, float m12, float m13, float m14,
	         float m21, float m22, float m23, float m24,
	         float m31, float m32, float m33, float m34,
	         float m41, float m42, float m43, float m44);

	void transpose();
	void inverse();
	void identity();

	//Getters
	Mat4 getTransposed() const;
	bool isIdentity() const;
	bool isZero() const;

	//Mat Mat Arithmentic Operations
	Mat4 operator+(const Mat4& rhs) const;
	Mat4 operator-(const Mat4& rhs) const;
	Mat4 operator*(const Mat4& rhs) const;
	Mat4 operator/(const Mat4& rhs) const;

	Mat4 operator*(float scalar) const;
	Mat4 operator/(float scalar) const;
	Mat4& operator*=(float scalar);
	Mat4& operator/=(float scalar);

	Mat4& operator+=(const Mat4& rhs);
	Mat4& operator-=(const Mat4& rhs);
	Mat4& operator*=(const Mat4& rhs);
	Mat4& operator/=(const Mat4& rhs);

	//Translate
	Mat4 translate(const Vector3& v) const;
	Mat4& translate(const Vector3& v);

	//Scale
	Mat4 scale(const Vector3& v) const;
	Mat4& scale(const Vector3& v);

	float operator [](int i) const;
	float& operator [](int i);

private:
	float determinant() const;
};

//cout operator
std::ostream& operator<<(std::ostream& os, const Mat4& m);