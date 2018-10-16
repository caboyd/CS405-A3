#include "Vector3.h"
#include "Mat4.h"
#include <iostream>
#include <string>
#include <vector>

#include "Math.h"
#include "tgaimage.h"

#include "Camera.h"
#include "Sphere.h"
#include "Polygon.h"
using std::cout;
using std::endl;

Ray rayConstruction(int i, int j);
int shading(const Vector3& pos, const Vector3& normal, float kd);
int rayTracing(const Ray& ray);
int rayObjectIntersection(const Ray& ray, HitRecord& hit_record);


constexpr unsigned ROWS = 512;
constexpr unsigned COLS = 512;
constexpr float focal = 50.0f / 1000.f; /* focal length simulating 50 mm lens */
constexpr float film_size = 35.0f / 1000.f;
constexpr float aspect_ratio = float(ROWS) / float(COLS);

Vector3 LRP{-10.0, 10.0, 2.0};
constexpr float Ip = 200.0;
constexpr float ambient = 20.0f;

TGAImage image(ROWS, COLS, TGAImage::RGB);

Vector3 VRP{1.0, 2.0, 3.5};
Vector3 VPN{0.0, -1.0, -2.5};
Vector3 VUP{0.0, 1.0, 0.0};
Camera camera(VRP, VPN, VUP, focal, film_size, 1.5);
Mat4 camera_to_world = camera.getViewMatrixInverse();


Sphere s({1, 1, 1}, 1.0, 0.75 );
Polygon poly({{0,0,0},{0,0,2},{2,0,2},{2,0,0}}, {0,1,0}, 0.8);
HitRecord g_hit_record;

int main()
{
	//Height
	for (int i = 0; i < ROWS; i++)
	{
		//Width
		for (int j = 0; j < COLS; j++)
		{
			//Construct th ray starting from center of projection, p0,
			// and passing through the pixel (x, y);
			Ray ray = rayConstruction(i, j);

			float c;
			if ((c = rayTracing(ray)) != -1)
				image.set(j, i, TGAColor((g_hit_record.normal.x+1)*128, (g_hit_record.normal.y+1)*128, (g_hit_record.normal.z+1)*255));
		}
	}

	image.write_tga_file("output.tga");
}


int rayTracing(const Ray& ray)
{
	// If the ray intersects with any object, this function call will
	// return a true value, and the data of the intersection are return
	// in P, N, and kd.
	HitRecord hit_record;
	int C;

	int found = rayObjectIntersection(ray, hit_record);

	if (found)
	{
			g_hit_record = hit_record;
		C = shading(hit_record.position, hit_record.normal, hit_record.kd);
		return (C);
	}

	return (-1);
}

int rayObjectIntersection(const Ray& ray, HitRecord& hit_record)
{
	if (s.hit(ray,hit_record))
	{
		return 1;
	}
	if( poly.hit(ray,hit_record))
	{
		return 1;
	}
	return 0;
}


Ray rayConstruction(int i, int j)
{
	// Step 1:
	// Map (j, i) in the screen coordinates to (xc, yc) in the
	// camera coordinates.
	float xmin = film_size / 2 * aspect_ratio;
	float ymin = film_size / 2;
	float xmax = -film_size / 2 * aspect_ratio;
	float ymax = -film_size / 2;

	float x = (xmax - xmin) * float(j) / float(COLS - 1) + xmin;
	float y = (ymax - ymin) * float(i) / float(ROWS - 1) + ymin;

	// Step 2:
	// Transform the origin (0.0, 0.0, 0.0) of the camera
	// coordinates to P0 in the world coordinates using the
	// transformation matrix Mcw. Note that the transformed result
	// should be VRP.
	Vector3 VRP = camera_to_world * Vector3(0);

	// Step 3:
	// Transform the point (xc, yc, f) on the image plane in
	// the camera coordinates to P1 in the world coordinates using
	// the transformation matrix Mcw.
	// V0 = P1 – P0;
	Vector3 V0 = (camera_to_world * Vector3(x, y, focal)) - VRP;
	// normalize V0 into unit length;
	V0.normalize();

	return Ray(VRP, V0);
};

int shading(const Vector3& pos, const Vector3& normal, float kd)
{
	Vector3 L = LRP - pos;
	L.normalize();

	int C = Ip * kd * normal.dot(L);
	if (normal.dot(L) < 0) C = 0;
	return C;
}
