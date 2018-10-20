#pragma once
class HitRecord
{
public:

	//Default to max float, indicating no object has been hit yet
	float t{};
	float kd{};
	Vector3 position{};
	Vector3 normal{};
	HitRecord() = default;
};


