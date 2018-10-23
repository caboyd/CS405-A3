#pragma once

class HitRecord;
class Ray;

//An abstract hittable object that is inherited so hittable objects
//can be in a list together
class Hittable
{
public:
	virtual ~Hittable() = default;

	//t_min is the minimum value along the ray that is allowed for an intersection to count.
	//negative t_min means the ray is behind the camera
	//this is the z near plane

	//t_max is the maximum value along the ray this is allowed for an intersection to count.
	//this is the z far plane
	virtual  bool hit(const Ray& ray, float t_min, float t_max, HitRecord& hit_record) = 0;
};

