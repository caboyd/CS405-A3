#pragma once

class HitRecord;
class Ray;

class Hittable
{
public:
	virtual ~Hittable() = default;
	virtual  bool hit(const Ray& ray, HitRecord& hit_record) = 0;
};

