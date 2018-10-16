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
};
