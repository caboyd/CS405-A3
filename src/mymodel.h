#pragma once
#include "Mat4.h"
#include "Sphere.h"
#include "Polygon.h"

Sphere sphere_obj(
	{1, 1, 1},	//center
	1.0f,		// radius
	0.75f);		// diffuse coefficient

Polygon polygon_obj(
	{{0.5,0,2},  {2.5, 0, 2}, {2.5, 0, 0},{0.5, 0,0}}, // 4 Vertices, CCW
	{0, 1, 0},	// normal vector		
	0.8f);		// diffuse coefficient

/* definition of the image buffer */
constexpr unsigned ROWS = 512;
constexpr unsigned COLS = 512;

unsigned char img[ROWS][COLS];

/* definition of window on the image plane in the camera coordinates */
/* They are used in mapping (j, i) in the screen coordinates into */
/* (x, y) on the image plane in the camera coordinates */
/* The window size used here simulates the 35 mm film. */

constexpr float focal = 50.0f / 1000.f; /* focal length simulating 50 mm lens */
constexpr float film_size = 35.f / 1000.f; /* film size simulating 35 mm film */
constexpr float aspect_ratio = float(ROWS) / float(COLS);

constexpr float xmin = film_size / 2 * aspect_ratio;
constexpr float ymin = film_size / 2;
constexpr float xmax = -film_size / 2 * aspect_ratio;
constexpr float ymax = -film_size / 2;

/* definition of the camera parameters */
Vector3 VRP{-0.5f, 0.8f, 4.5f};
Vector3 VPN{1.2f, 0, -2.5f};
Vector3 VUP{0.0f, 1.0f, 0.0f};

/* definition of light source */
Vector3 LRP{-10.0, 10.0, 2.0};	/* light position */
constexpr float Ip = 200.0;		/* intensity of the point light source */

Mat4 Mwc = Mat4::fromIdentity();
Mat4 Mcw = Mat4::fromIdentity();


