#include <iostream>

#include "mymodel.h"
#include "Math.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

Ray rayConstruction(int i, int j);
int shading(const Vector3& pos, const Vector3& normal, float kd);
int rayTracing(const Ray& ray);
int rayObjectIntersection(const Ray& ray, HitRecord& hit_record);

int main()
{
	//Initialization of Camera Matrices
	coordinateSystemTransformationMatricesFromPositionNormalUp(VRP,VPN,VUP,Mwc,Mcw);

	//Height
	for (unsigned int i = 0; i < ROWS; i++)
	{
		//Width
		for (unsigned int j = 0; j < COLS; j++)
		{
			//Construct th ray starting from center of projection, p0,
			// and passing through the pixel (x, y);
			Ray ray = rayConstruction(i, j);

			int c;
			//Color the image using the shading value generated
			if ((c = rayTracing(ray)) != -1)
				img[i][j] = c;
		}
	}

	stbi_write_bmp("output.bmp",ROWS,COLS,1,img);
}


int rayTracing(const Ray& ray)
{
	// If the ray intersects with any object, this function call will
	// return a true value, and the data of the intersection are return
	// in P, N, and kd.

	//HitRecord is the container for P, N, and kd
	// P is intersection point
	// N is normal of object at that point
	// kd is the diffuse reflection coefficient
	HitRecord hit_record;

	//If the ray intersects with any objects in the scene
	//return the color using the shading value generated
	if (rayObjectIntersection(ray, hit_record))
		return shading(hit_record.position, hit_record.normal, hit_record.kd);

	//No intersection
	return -1;
}


// Construction of ray
// Input: pixel index (i, j) in the screen coordinates
// Output: P0, V0 (for parametric ray equation P = P0 + V0*t)
// in the world coordinates.
Ray rayConstruction(int i, int j)
{
	// Step 1:
	// Map (j, i) in the screen coordinates to (xc, yc) in the
	// camera coordinates.

	//(0,0) is upper left corner
	float x = (xmax - xmin) * float(j) / float(COLS - 1) + xmin;
	float y = (ymax - ymin) * float(i) / float(ROWS - 1) + ymin;

	// Step 2:
	// Transform the origin (0.0, 0.0, 0.0) of the camera
	// coordinates to P0 in the world coordinates using the
	// transformation matrix Mcw. Note that the transformed result
	// should be VRP.
	Vector3 P0 = Mcw * Vector3(0);

	// Step 3:
	// Transform the point (xc, yc, f) on the image plane in
	// the camera coordinates to P1 in the world coordinates using
	// the transformation matrix Mcw.
	// V0 = P1 – P0;
	Vector3 P1 = Mcw * Vector3(x,y,focal);
	Vector3 V0 = P1 - P0;
	// normalize V0 into unit length;
	V0.normalize();

	return Ray(P0, V0);
};


// Ray-Object Intersection
// Input: ray – ray.origin, ray.direction
// Output: hit_record - the nearest intersection point P[3] if found, 
// along with N[3], the surface normal at that point, and
// kd, the diffuse reflection coefficient of the surface.
// Note: In a general system, the objects should be stored a list
// structure. A loop will scan through each object in the list. The
// nearest intersection point is found. In our case, we will have only
// two hard-coded objects: a sphere and a polygon.
int rayObjectIntersection(const Ray& ray, HitRecord& hit_record)
{
	bool hit_anything = false;

	//The closest intersection point found.
	//The t value in the equation Point = Ray.origin + t * Ray.direction
	float closest_so_far_t = FLT_MAX;

	//Compute sphere intersection
	if (sphere_obj.hit(ray, 0.0001f, closest_so_far_t, hit_record))
	{
		hit_anything = true;
		closest_so_far_t = hit_record.t;
	}

	//Compute polygon intersection
	//Will only be true if polygon is hit and is closer than sphere
	if (polygon_obj.hit(ray,  0.0001f, closest_so_far_t, hit_record))
	{
		hit_anything = true;
	}

	return hit_anything;
}


int shading(const Vector3& pos, const Vector3& normal, float kd)
{
	//Normalized light direction vector
	Vector3 L = LRP - pos;
	L.normalize();

	//Phong lighting model
	//The amount of light is equal to 
	//light intenstity * diffuse coefficient * cos angle between normal and light
	int color = int(Ip * kd * normal.dot(L));

	//Color may be negative if angle between normal and light >90 degrees
	if (color < 0) 
		color = 0;
	return color;
}
