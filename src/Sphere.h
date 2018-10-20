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

	bool hit(const Ray& ray, float t_min, float t_max, HitRecord& hit_record) override;

private:
	bool sphereIntersectionMethod1(const Ray& ray, float t_min, float t_max, HitRecord& hit_record);
	bool sphereIntersectionMethod2(const Ray& ray, float t_min, float t_max, HitRecord& hit_record);
};

inline bool Sphere::hit(const Ray& ray, float t_min, float t_max, HitRecord& hit_record)
{
	return sphereIntersectionMethod2(ray, t_min, t_max, hit_record);
}


inline bool Sphere::sphereIntersectionMethod1(const Ray& ray, float t_min, float t_max, HitRecord& hit_record)
{
	Vector3 D = ray.direction;

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
		if (t0 > 0.005)
		{
			hit_record.t = t0;
			hit_record.position = ray.point_at_parameter(t0);
			hit_record.normal = ((ray.origin + t0 * D) - center).getNormalized();
			hit_record.kd = kd;
			return true;
		}

		if (t1 > 0.005)
		{
			hit_record.t = t1;
			hit_record.position = (ray.origin + t1 * D);
			hit_record.normal = ((ray.origin + t1 * D) - center).getNormalized();
			hit_record.kd = kd;
			return true;
		}
	}
	return false;
}



inline bool Sphere::sphereIntersectionMethod2(const Ray& ray, float t_min, float t_max, HitRecord& hit_record)
{
	//https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
	//Equation of Sphere
	// (X - X0)^2 + (Y - Y0)^2 + (Z - Z0)^2 = R^2
	// Where center of sphere is (X0, Y0, Z0)

	//Equation of Line
	// ray = origin + t * direction

	//Substitute ray(X,Y,Z) into Equation of Sphere for X,Y,Z
	// (origin + t * direction - center)^2 = radius^2

	// (origin.x + t * direction.x - X0)^2 + 
	// (origin.y + t * direction.y - Y0)^2 +
	// (origin.z + t * direction.z - Z0)^2 = radius^2

	//Expand
	// t^2 (direction * direction ) + 2t( direction * (origin - center)) + (origin - center) * (origin - center) = radius^2

	// Which is the quadratic equation
	// at^2 + bt + c = 0
	// Where
	// a = (direction * direction) dot product
	// b = 2 (direction * (origin-center)) dot product
	// c = (origin-center) * (origin-center) - radius^2

	float a = ray.direction.dot(ray.direction);
	float b = 2 * ray.direction.dot(ray.origin - center);
	float c = (ray.origin - center).dot(ray.origin - center) - radius * radius;

	//Use quadratic formula
	// (-b +- sqrt(b^2 - 4ac)) / 2a

	//find discriminant first to determine type of roots
	// discriminant = (b^2 - 4ac)
	float discriminant = b * b - 4 * a * c;

	if (discriminant < 0)
	{
		//No Roots	
	}
	else if (discriminant > 0)
	{
		//One or Two Roots
		float t0 = (-b - sqrt(discriminant)) / 2 * a;
		if (t0 > t_min && t0 < t_max)
		{
			hit_record.t = t0;
			hit_record.position = (ray.origin + t0 * ray.direction);
			hit_record.normal = (hit_record.position - center).getNormalized();
			hit_record.kd = kd;
			return true;
		}

		float t1 = (-b + sqrt(discriminant)) / 2 * a;
		if (t1 > t_min && t1 < t_max)
		{
			hit_record.t = t1;
			hit_record.position = (ray.origin + t1 * ray.direction);
			hit_record.normal = (hit_record.position - center).getNormalized();
			hit_record.kd = kd;
			return true;
		}
	}
	return false;
}
