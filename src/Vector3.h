#pragma once
#include <ostream>


class Vector3
{
public:
	union
	{
		float data[3];
		struct{	float x, y, z;};
	};



	static const Vector3 ZERO;

	//Constructors
	Vector3();
	explicit Vector3(float value);
	Vector3(float x, float y, float z);

	//Copy Constructor
	Vector3(const Vector3& other);

	//Copy Assignment
	Vector3& operator=(const Vector3& other);

	//Comparison
	bool operator==(const Vector3& other) const;

	//Mutators
	void set(float x, float y, float z);
	void normalize();
	void negate();


	//Getters
	bool isZero() const;
	float length() const;
	float getSquaredLength() const;
	Vector3 getNormalized() const;
	float distance(const Vector3& v) const;
	float distanceSquared(const Vector3& v) const;
	int getLargestComponentIndex()const;

	//Negation
	Vector3 operator-() const;

	//Vector Vector Arithmetic operations
	Vector3 operator+(const Vector3& v) const;
	Vector3 operator-(const Vector3& v) const;
	Vector3 operator*(const Vector3& v) const;
	Vector3 operator/(const Vector3& v) const;

	//Vector Vector Assignment Arithmetic operations
	Vector3& operator+=(const Vector3& v);
	Vector3& operator-=(const Vector3& v);
	Vector3& operator*=(const Vector3& v);
	Vector3& operator/=(const Vector3& v);

	//Vector Float Arithmetic operations
	Vector3 operator-(float f) const;
	Vector3 operator+(float f) const;
	Vector3 operator/(float scalar) const;
	Vector3 operator*(float scalar) const;

	//Vector Float Assignment Arithmetic operations
	Vector3& operator+=(float f);
	Vector3& operator-=(float f);
	Vector3& operator*=(float f);
	Vector3& operator/=(float f);

	float operator[](int i) const;
	float& operator[](int i);

	//Cross Product
	Vector3 cross(const Vector3& v) const;

	//Dot Product
	float dot(const Vector3& v) const;
};

//cout operator
std::ostream& operator<<(std::ostream& os, const Vector3& v);

//Float Vector Arithmetic operations
inline Vector3 operator+(float f, const Vector3& v) { return v + f; };
inline Vector3 operator-(float f, const Vector3& v) { return Vector3(f - v.x, f - v.y, f - v.z); };
inline Vector3 operator*(float f, const Vector3& v) { return v * f; };
inline Vector3 operator/(float f, const Vector3& v) { return Vector3(f / v.x, f / v.y, f / v.z); };
