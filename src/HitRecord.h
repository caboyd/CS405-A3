#pragma once

/**
 * Contains information about a ray hit
 * 
 */
class HitRecord
{
public:

	//Default to max float, indicating no object has been hit yet
	float t{};			//where along the ray the intersection happened
	float kd{};			//The diffused reflection coefficient of the object that was hit
	Vector3 position{};	//The position in world coordinates of the intersection
	Vector3 normal{};	//The normal of the object at the intersection point
	HitRecord() = default;
};


