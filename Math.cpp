#include "Vector3.h"

inline void OrthoNormalVectors(const Vector3& v1, const Vector3& v2, Vector3& u, Vector3& v, Vector3& n)
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
