#include "Vector3.h"
#include "Mat4.h"
#include <iostream>

void orthoNormalVectors(const Vector3& v1, const Vector3& v2, Vector3& u, Vector3& v, Vector3& n)
{
	const Vector3& preferred_up = v2;

	const Vector3 forward = v1.getNormalized();
	n = forward;

	Vector3 right = forward.cross(preferred_up);
	right.normalize();
	u = right;

	const Vector3 up = right.cross(forward);
	v = up;
}


void worldToCameraMatricesFromPositionNormalUp(const Vector3& VRP, const Vector3& VPN, const Vector3& VUP,
                                               Mat4& world_to_camera,
                                               Mat4& camera_to_world)
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
	world_to_camera = R * T;
	
	Mat4 R_inv = Mat4(u.x, v.x, n.x, 0.f,
	                  u.y, v.y, n.y, 0.f,
	                  u.z, v.z, n.z, 0.f,
	                  0.f, 0.f, 0.f, 1.f);

	Mat4 T_inv = Mat4(1.f, 0.f, 0.f, VRP.x,
	                  0.f, 1.f, 0.f, VRP.y,
	                  0.f, 0.f, 1.f, VRP.z,
	                  0.f, 0.f, 0.f, 1.f);

	//4x4 Matrix Multiplication
	camera_to_world = T_inv * R_inv;
}

void WorldToLightMatricesFromPositionNormalUp(const Vector3& LRP, const Vector3& LPN, const Vector3& LUP,
                                              Mat4& world_to_light,
                                              Mat4& light_to_world)
{
	Vector3 u, v, n;
	orthoNormalVectors(LPN, LUP, u, v, n);

	Mat4 R = Mat4(u.x, u.y, u.z, 0.f,
	              v.x, v.y, v.z, 0.f,
	              n.x, n.y, n.z, 0.f,
	              0.f, 0.f, 0.f, 1.f);

	Mat4 T = Mat4(1.f, 0.f, 0.f, -LRP.x,
	              0.f, 1.f, 0.f, -LRP.y,
	              0.f, 0.f, 1.f, -LRP.z,
	              0.f, 0.f, 0.f, 1.f);

	//4x4 Matrix Multiplication
	world_to_light = R * T;

	Mat4 R_inv = Mat4(u.x, v.x, n.x, 0.f,
	                  u.y, v.y, n.y, 0.f,
	                  u.z, v.z, n.z, 0.f,
	                  0.f, 0.f, 0.f, 1.f);

	Mat4 T_inv = Mat4(1.f, 0.f, 0.f, LRP.x,
	                  0.f, 1.f, 0.f, LRP.y,
	                  0.f, 0.f, 1.f, LRP.z,
	                  0.f, 0.f, 0.f, 1.f);

	//4x4 Matrix Multiplication
	light_to_world = T_inv * R_inv;
}
