/**
 *	
 *	Name: Chris Boyd
 *	Assignment #2
 *	October 24,2018
 *	
 *	Simple greyscale ray tracer with a sphere and quad in scene
 *	Outputs a 512x512 BMP file using stb_image_write.h (https://github.com/nothings/stb)
 *
 */

#include "mymodel.h"
#include "Math.h"
#include "Ray.h"
#include "HitRecord.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <fstream>
#include <iostream>

Ray rayConstruction(int i, int j);
void write_raw(const char* filename, int bytes, const void* data);
inline int clamp(int value, int min, int max);

void readDensityData(const char* filename, int rows, int cols, int layers, void* data);
void computeShadingVolume();
int rayBoxIntersection(const Ray& ray, float ts[2]);
int volumeRayTracing(const Ray& ray, float ts[2], float min_density, float max_density);
void volumeRayTraceSpecificDensities(const char* filename, float min_density, float max_density);
float trilinearInterpolation(unsigned char buffer[DEN_LAYERS][DEN_ROWS][DEN_COLS], const Vector3& p);
inline float bilinearInterpolation(float v00, float v10, float v01, float v11, float tx, float ty);
inline float lerp(float v0, float v1, float t);

int main()
{
	//Initialization of Camera Matrices
	coordinateSystemTransformationMatricesFromPositionNormalUp(VRP, VPN, VUP, Mwc, Mcw);

	readDensityData("smallHead.den", DEN_ROWS, DEN_COLS, DEN_LAYERS, density_buffer);
	computeShadingVolume();

	//RIGHT SIDE
	VRP = {286, 64, 64};
	VPN = {-1, 0, 0};
	VUP = {0, 0, -1};
	volumeRayTraceSpecificDensities("right_soft_tissue", 30.f, 55.f);
	volumeRayTraceSpecificDensities("right_skin", 55.f, 90.f);
	volumeRayTraceSpecificDensities("right_bone", 110.f, 255.f);

	VRP = {64, -286, 64};
	VPN = {0, 1, 0};
	VUP = {0, 0, 1};
	volumeRayTraceSpecificDensities("front_soft_tissue", 30.f, 55.f);
	volumeRayTraceSpecificDensities("front_skin", 55.f, 90.f);
	volumeRayTraceSpecificDensities("front_bone", 110.f, 255.f);

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
	float x = (xmax - xmin) * float(j) / float(IMG_COLS - 1) + xmin;
	float y = (ymax - ymin) * float(i) / float(IMG_ROWS - 1) + ymin;

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
	Vector3 P1 = Mcw * Vector3(x, y, focal);
	Vector3 V0 = P1 - P0;
	// normalize V0 into unit length;
	V0.normalize();

	return Ray(P0, V0);
};


//Writes the raw bytes to a file
void write_raw(const char* filename, int bytes, const void* data)
{
	//open outfile in out mode, overrite mode, and binary mode
	std::ofstream os(filename, std::ofstream::out | std::ofstream::trunc | std::ios::binary);
	if (os.is_open())
	{
		for (int i = 0; i < bytes; i++)
			os << static_cast<const unsigned char*>(data)[i];
	}
	os.close();
}

//Clamps a value between min and max
inline int clamp(int value, int min, int max)
{
	const int result = value < min ? min : value;
	return result > max ? max : result;
}


//Copies the density data as a char array into the data blob
void readDensityData(const char* filename, int rows, int cols, int layers, void* data)
{
	//char pointer to data
	char* buffer = static_cast<char*>(data);

	std::ifstream is(filename, std::ifstream::binary);
	if (is)
	{
		//Assert length of file is same as expected length of file
		is.seekg(0, is.end);
		int length = is.tellg();
		is.seekg(0, is.beg);
		assert(length == rows*cols*layers);

		//Copy the file into the buffer as char array
		is.read(buffer, length);
	}
}

/*
 * Computes the shading buffer using partial derivatives generated fron density data
 */
void computeShadingVolume()
{
	Vector3 normal;
	int color = 0;
	// boundary layers are skipped.
	for (unsigned int z = 1; z < (DEN_LAYERS - 1); z++)
		for (unsigned int y = 1; y < (DEN_ROWS - 1); y++)
			for (unsigned int x = 1; x < (DEN_COLS - 1); x++)
			{
				// compute the partial derivative at [z, y, x]
				float dfdx = (float(density_buffer[z][y][x + 1]) - float(density_buffer[z][y][x - 1])) / 2.0f;
				float dfdy = (float(density_buffer[z][y + 1][x]) - float(density_buffer[z][y - 1][x])) / 2.0f;
				float dfdz = (float(density_buffer[z + 1][y][x]) - float(density_buffer[z - 1][y][x])) / 2.0f;

				//move values into Vector
				//The partial derivative is the normal at that point
				normal.set(dfdx, dfdy, dfdz);

				// if the magnitude of the gradient is less than a
				// pre-specified epsilon value, set shading to 0.
				// otherwise, compute the diffuse shading.
				// save the result in the shading volume.
				if (normal.getSquaredLength() > EPSILON * EPSILON)
				{
					normal.normalize();
					//color = Light intensity * NdotL
					color = int(Ip * normal.dot(light_dir));
				}

				//colors must be clamped from 0 to 255
				color = clamp(color, 0, 255);

				shading_buffer[z][y][x] = color;
			}
}

/*
 * Returns the t values in P(t) = P0 + t*V of
 * the intersection points of the ray intersecting the box
 */
int rayBoxIntersection(const Ray& ray, float ts[2])
{
	//index for ts
	int i = 0;

	//6 Bounding Box sides
	float TOP = DEN_ROWS-1;
	float BOTTOM = 0;

	float FRONT = DEN_LAYERS-1;
	float BACK = 0;

	float LEFT = 0;
	float RIGHT = DEN_COLS-1;

	//Calculates the position the ray would intersect with each plane
	//If there is an intersection the t value is stored in ts[];
	//When there are two intersection we can early return


	//LEFT SIDE
	Vector3 p = ray.point_at_plane_x(LEFT); //Point of intersection
	//If intersection is within the quad save t
	if (p.y > BOTTOM && p.y < TOP && p.z > BACK && p.z < FRONT)
		ts[i] = ray.parameter_at_plane_x(LEFT), i++;

	//RIGHT SIDE
	p = ray.point_at_plane_x(RIGHT);
	if (p.y > BOTTOM && p.y < TOP && p.z > BACK && p.z < FRONT)
		ts[i] = ray.parameter_at_plane_x(RIGHT), i++;
	if (i == 2) return i;

	//BOTTOM SIDE
	p = ray.point_at_plane_y(BOTTOM);
	if (p.x > LEFT && p.x < RIGHT && p.z > BACK && p.z < FRONT)
		ts[i] = ray.parameter_at_plane_y(BOTTOM), i++;
	if (i == 2) return i;

	//TOP SIDE
	p = ray.point_at_plane_y(TOP);
	if (p.x > LEFT && p.x < RIGHT && p.z > BACK && p.z < FRONT)
		ts[i] = ray.parameter_at_plane_y(TOP), i++;
	if (i == 2) return i;

	//BACK SIDE
	p = ray.point_at_plane_z(BACK);
	if (p.x > LEFT && p.x < RIGHT && p.y > BOTTOM && p.y < TOP)
		ts[i] = ray.parameter_at_plane_z(BACK), i++;
	if (i == 2) return i;

	//FRONT SIDE
	p = ray.point_at_plane_z(FRONT);
	if (p.x > LEFT && p.x < RIGHT && p.y > BOTTOM && p.y < TOP)
		ts[i] = ray.parameter_at_plane_z(FRONT), i++;

	return i;
}

int volumeRayTracing(const Ray& ray, float ts[2], float min_density, float max_density)
{
	float Dt = 1.f; // the interval for sampling along the ray
	float C = 0.0; // for accumulating the shading value
	float T = 1.0; // for accumulating the transparency

	/* Marching through the CT volume from t1 to t2
	*  by step size Dt.
	*/
	float t1 = ts[0];
	float t2 = ts[1];

	//We want to march in front-to-back order
	//so swap if t1 is larger
	if (t1 > t2) std::swap(t1, t2);


	//March front-to-back
	for (float t = t1; t <= t2; t += Dt)
	{
		/* Compute the 3D coordinates of the current
		 * sample position in the volume:
		 * x = x(t);
		 * y = y(t);
		 * z = z(t);
		 */
		Vector3 point = ray.point_at_parameter(t);

		/* Obtain the shading value C and opacity value A
		 * from the shading volume and CT volume, respectively,
		 * by using tri-linear interpolation. 
		 */
		float color = trilinearInterpolation(shading_buffer, point);
		float den = trilinearInterpolation(density_buffer, point);
		float alpha = den / 255.0f;

		/*
		 * Accumulate Color and Transparency onl when density
		 * value is between these bounds.
		 * This allows to only view specific densities
		 * 
		 * 15 - 60 : Soft tissue
		 * 60 - 110 : Skin
		 * 110 - 255 : Bone
		 * 220 - 255 : Hard Bone (teeth)
		 */
		if (den > min_density && den < max_density)
		{
			C += T * alpha * color;
			T *= (1.0f - alpha);
		}
		
	}

	int color = clamp((int)C, 0, 255);

	return color;
}

void volumeRayTraceSpecificDensities(const char* filename, float min_density, float max_density)
{
	//Height
	for (unsigned int i = 0; i < IMG_ROWS; i++)
	{
		//Width
		for (unsigned int j = 0; j < IMG_COLS; j++)
		{
			//Construct th ray starting from center of projection, p0,
			// and passing through the pixel (x, y);
			Ray ray = rayConstruction(i, j);
			float ts[2] = {};

			int n = rayBoxIntersection(ray, ts);
			if (n == 2)
				img[i][j] = volumeRayTracing(ray, ts, min_density, max_density);
		}
	}

		//Output file as bmp picture
	char str[80];
	strcpy_s(str, filename);
	strcat_s(str,".bmp");

	stbi_write_bmp(str, IMG_COLS, IMG_ROWS, 1, img);

	strcpy_s(str, filename);
	strcat_s(str,".raw");
	//Output file as raw buffer data
	//write_raw(str, sizeof(img), img);

}

float trilinearInterpolation(unsigned char buffer[128][128][128], const Vector3& p)
{
	//INDICES
	auto posx = (unsigned char)ceil(p.x);
	auto negx = (unsigned char)floor(p.x);

	auto posy = (unsigned char)ceil(p.y);
	auto negy = (unsigned char)floor(p.y);

	auto posz = (unsigned char)ceil(p.z);
	auto negz = (unsigned char)floor(p.z);

	//TOP 4 corners
	unsigned char posx_posz = buffer[posz][posy][posx];
	unsigned char negx_posz = buffer[posz][posy][negx];
	unsigned char posx_negz = buffer[negz][posy][posx];
	unsigned char negx_negz = buffer[negz][posy][negx];

	float top = bilinearInterpolation(negx_posz, posx_posz, negx_negz, posx_negz, p.x - negx, 1.0f - (p.z - negz));

	//BOTTOM 4 corners
	posx_posz = buffer[posz][negy][posx];
	negx_posz = buffer[posz][negy][negx];
	posx_negz = buffer[negz][negy][posx];
	negx_negz = buffer[negz][negy][negx];

	float bottom = bilinearInterpolation(negx_posz, posx_posz, negx_negz, posx_negz, p.x - negx, 1.0f - (p.z - negz));

	return lerp(bottom, top, p.y - negy);
}

inline float bilinearInterpolation(float v00, float v10, float v01, float v11, float tx, float ty)
{
	return lerp(lerp(v00, v10, tx), lerp(v01, v11, tx), ty);
}


inline float lerp(float v0, float v1, float t)
{
	return (1 - t) * v0 + t * v1;
}
