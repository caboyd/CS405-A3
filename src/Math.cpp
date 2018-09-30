#include "Vector3.h"
#include "Mat4.h"
#include <iostream>


void orthoNormalVectors(const Vector3& v1, const Vector3& v2, Vector3& u, Vector3& v, Vector3& n)
{
	//v2 is the preferred up vector, but not necessarily the exact up vector
	const Vector3& preferred_up = v2;

	//The normal is the direction of v1 and will be along
	//the Z axis in our new coordinate system.
	const Vector3 z_forward = v1.getNormalized();
	//n is called the forward vector because it is the direction the normal is facing
	n = z_forward;


	//Y cross Z gives +X because we are in right handed coordinate system
	Vector3 x_right = preferred_up.cross(z_forward);
	x_right.normalize();
	//u is called the right vector because positive X is to the right
	u = x_right;

	//Z cross X gives +Y because we are in right handed coordinate system
	const Vector3 y_up = z_forward.cross(x_right);
	//v is called the up vector because positive Y is up.
	v = y_up;
}

void coordinateSystemTransformationMatricesFromPositionNormalUp(
	const Vector3& VRP, const Vector3& VPN,
	const Vector3& VUP,
	Mat4& transformation_matrix,
	Mat4& inverse_transformation_matrix)
{
	Vector3 u, v, n;
	orthoNormalVectors(VPN, VUP, u, v, n);

	Mat4 R = Mat4(u.x, u.y, u.z, 0.f,
	              v.x, v.y, v.z, 0.f,
	              n.x, n.y, n.z, 0.f,
	              0.f, 0.f, 0.f, 1.f);

	Mat4 T = Mat4(1.f, 0.f, 0.f, -VRP.x,
	              0.f, 1.f, 0.f, -VRP.y,
	              0.f, 0.f, 1.f, -VRP.z,
	              0.f, 0.f, 0.f, 1.f);

	//4x4 Matrix Multiplication
	transformation_matrix = R * T;

	Mat4 R_inv = Mat4(u.x, v.x, n.x, 0.f,
	                  u.y, v.y, n.y, 0.f,
	                  u.z, v.z, n.z, 0.f,
	                  0.f, 0.f, 0.f, 1.f);

	Mat4 T_inv = Mat4(1.f, 0.f, 0.f, VRP.x,
	                  0.f, 1.f, 0.f, VRP.y,
	                  0.f, 0.f, 1.f, VRP.z,
	                  0.f, 0.f, 0.f, 1.f);

	//4x4 Matrix Multiplication
	inverse_transformation_matrix = T_inv * R_inv;
}
