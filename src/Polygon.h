#pragma once
#include "Hittable.h"
#include "Vector3.h"
#include "Ray.h"
#include "HitRecord.h"
#include <vector>

int pnpoly(std::vector<Vector3> vertices, Vector3 point, int x, int y);

class Polygon : Hittable
{
public:
	std::vector<Vector3> vertices;
	Vector3 normal;
	float kd;

	Polygon(std::vector<Vector3> verts, const Vector3& normal, float kd)
	{
		this->vertices = verts;
		this->normal = normal;
		this->kd = kd;
	}

	bool hit(const Ray& ray, HitRecord& hit_record) override;
};

inline bool Polygon::hit(const Ray& ray, HitRecord& hit_record)
{
	float D = -normal.dot(vertices[0]);

	float NdotO = normal.dot(ray.origin);
	float NdotV = normal.dot(ray.direction);

	if (NdotV == 0.0) return false;

	float t = - (NdotO + D) / NdotV;

	Vector3 point = ray.point_at_parameter(t);

	//check if intersection point is in polygon
	int index = normal.getLargestComponentIndex();
	int x = 0, y = 1;
	if(index == 0) x = 1, y = 2;
	if(index == 1) x = 0, y = 2;
	if(index == 2) x = 0, y = 1;
	if (pnpoly(vertices, point, x, y))
	{
		hit_record.normal = normal;
		hit_record.kd = kd;
		hit_record.position = point;
		hit_record.t = t;
		return  true;
	}
	return  false;
}


inline int pnpoly(std::vector<Vector3> vertices, Vector3 point, int x, int y)
{
	int c = 0;
	for (int i = 0, j = vertices.size() - 1; i < vertices.size(); j = i++)
	{
		if ((vertices[i][y] > point[y]) != (vertices[j][y] > point[y]) &&
			(point[x] < (vertices[j][x] - vertices[i][x]) * (point[y] - vertices[i][y]) / (vertices[j][y] - vertices[i][y]) +
				vertices[i][x]))
		{
			c = !c;
		}
	}
	return c;
}
