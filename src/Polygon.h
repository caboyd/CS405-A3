#pragma once
#include "Hittable.h"
#include "Vector3.h"
#include "Ray.h"
#include "HitRecord.h"
#include <vector>



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

	bool hit(const Ray& ray, float t_min, float t_max, HitRecord& hit_record) override;
private:
	static bool pointInPolygon(std::vector<Vector3> vertices, Vector3 point, int x, int y);
};

inline bool Polygon::hit(const Ray& ray, float t_min, float t_max, HitRecord& hit_record)
{
	float D = -normal.dot(vertices[0]);

	float NdotO = normal.dot(ray.origin);
	float NdotV = normal.dot(ray.direction);

	if (NdotV == 0.0) return false;

	float t = - (NdotO + D) / NdotV;
	if (t < t_min || t > t_max) return false;

	Vector3 point = ray.point_at_parameter(t);

	//check if intersection point is in polygon
	int index = normal.getLargestComponentIndex();
	int x = 0, y = 1;
	//Project on Z plane
	if (index == 0) x = 1, y = 2;
	//Project on Y plane
	if (index == 1) x = 0, y = 2;
	//Project on Z plane
	if (index == 2) x = 0, y = 1;

	//Point in polygon test using 2D projection
	if (pointInPolygon(vertices, point, x, y))
	{
		hit_record.normal = normal;
		hit_record.kd = kd;
		hit_record.position = point;
		hit_record.t = t;
		return true;
	}
	return false;
}


//Algorithm from https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
inline bool Polygon::pointInPolygon(std::vector<Vector3> vertices, Vector3 point, int x, int y)
{
	bool intersections_is_odd = false;
	for (size_t i = 0, j = vertices.size() - 1; i < vertices.size(); j = i++)
	{
		//One point above point, one point below point
		//implies intersection with x axis
		if ((vertices[i][y] > point[y]) != (vertices[j][y] > point[y]))
		{
			//x_intersection is the x along the line created by the two vertices when y is 0
			const float x_intersection = (vertices[j][x] - vertices[i][x]) * (point[y] - vertices[i][y]) / (vertices[j][
				y] - vertices[i][y]) + vertices[i][x];

			//If intersection on x-axis is to the right of the point then we count it as an intersection
			if (point[x] < x_intersection)
			{
				intersections_is_odd = !intersections_is_odd;
			}
		}
	}
	return intersections_is_odd;
}
