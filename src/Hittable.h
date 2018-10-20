#pragma once

class HitRecord;
class Ray;

class Hittable
{
public:
	virtual ~Hittable() = default;
	virtual  bool hit(const Ray& ray, float t_min, float t_max, HitRecord& hit_record) = 0;
};

