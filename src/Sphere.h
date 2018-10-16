#pragma once
#include "Hittable.h"
#include "Vector3.h"
#include "Ray.h"
#include "HitRecord.h"

class Sphere : Hittable
{
public:
	Vector3 center;
	float radius;
	float kd;

	Sphere(Vector3 center, float radius, float kd)
	{
		this->center = center;
		this->radius = radius;
		this->kd = kd;
	}

	bool hit(const Ray& ray, HitRecord& hit_record) override;
};

inline bool Sphere::hit(const Ray& ray, HitRecord& hit_record)
{
	Vector3 D =ray.direction;

	//sphere_center_to_ray_origin
	Vector3 L = this->center - ray.origin;

	//Project L onto D
	float tca = L.dot(D);

	//distance_from_ray_to_sphere
	float distance = sqrt(L.getSquaredLength() - (tca * tca));

	if (distance > this->radius)
	{
		//No Intersection
	}
	else if (distance <= this->radius)
	{
		//One Or Two Intersections
		float thc = sqrt((radius * radius) - (distance * distance));

		//Intersection points
		float t0 = tca - thc;
		float t1 = tca + thc;
		if(t0 > 0.005)
		{
			hit_record.t = t0;
			hit_record.position = ray.point_at_parameter(t0);
			hit_record.normal = ((ray.origin + t0 * D) - center).getNormalized();
			hit_record.kd = kd;
			return true;
		}

		if(t1 > 0.005)
		{
			hit_record.t = t1;
			hit_record.position = (ray.origin + t1 * D);
			hit_record.normal = ((ray.origin + t1 * D) - center).getNormalized();
			hit_record.kd = kd;
		//	return true;
		}
	}
	return false;
}
