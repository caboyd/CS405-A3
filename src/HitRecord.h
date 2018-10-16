#pragma once

#include "LambertianMaterial.h"

class HitRecord
{
public:

	float t{};
	float kd;
	Vector3 position{};
	Vector3 normal{};
	HitRecord() = default;
};


