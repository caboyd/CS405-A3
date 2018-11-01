#pragma once

#include "Vector3.h"

class Ray
{
public:
	Vector3 origin{};
	Vector3 direction{};

	Ray() = default;

	Ray(const Vector3& origin, const Vector3& direction)
	{
		this->origin = origin;
		this->direction = direction;
	}

	Vector3 point_at_parameter(float t) const
	{
		return origin + t * direction;
	}

	float parameter_at_plane_x(const float x) const
	{
		return ((x - origin.x) / direction.x);
	}

	float parameter_at_plane_y(const float y) const
	{
		return ((y - origin.y) / direction.y);
	}


	float parameter_at_plane_z(const float z) const
	{
		return ((z - origin.z) / direction.z);
	}

	Vector3 point_at_plane_x(const float x) const
	{
		return origin + ((x - origin.x) / direction.x) * direction;
	}

	Vector3 point_at_plane_y(const float y) const
	{
		return origin + ((y - origin.y) / direction.y) * direction;
	}


	Vector3 point_at_plane_z(const float z) const
	{
		return origin + ((z - origin.z) / direction.z) * direction;
	}
};
