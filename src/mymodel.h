#pragma once
#include "Mat4.h"
#include "Vector3.h"

/* definition of the image buffer */
constexpr unsigned IMG_ROWS = 512;
constexpr unsigned IMG_COLS = 512;

unsigned char img[IMG_ROWS][IMG_COLS];

constexpr unsigned DEN_ROWS = 128;
constexpr unsigned DEN_COLS = 128;
constexpr unsigned DEN_LAYERS = 128;

unsigned char density_buffer[DEN_LAYERS][DEN_ROWS][DEN_COLS];
unsigned char shading_buffer[DEN_LAYERS][DEN_ROWS][DEN_COLS];

const float EPSILON = 18.f;

/* definition of window on the image plane in the camera coordinates */
/* They are used in mapping (j, i) in the screen coordinates into */
/* (x, y) on the image plane in the camera coordinates */
/* The window size used here simulates the 35 mm film. */

constexpr float focal = 50.0f / 1000.f; /* focal length simulating 50 mm lens */
constexpr float film_size = 35.f / 1000.f; /* film size simulating 35 mm film */
constexpr float aspect_ratio = float(IMG_ROWS) / float(IMG_COLS);

//Coordinates of the screen
constexpr float xmin = film_size / 2 * aspect_ratio;
constexpr float ymin = film_size / 2;
constexpr float xmax = -film_size / 2 * aspect_ratio;
constexpr float ymax = -film_size / 2;

/* definition of the camera parameters */
Vector3 VRP{286, 64, 64};
Vector3 VPN{-1, 0, 0};
Vector3 VUP{0, 0, -1};

/* definition of light source */
Vector3 light_dir{-0.577f, 0.577f, 0.577f};	/* light position */

constexpr float Ip = 255.0;		/* intensity of the point light source */

//Matrices for world to camera and camera to world
Mat4 Mwc = Mat4::fromIdentity();
Mat4 Mcw = Mat4::fromIdentity();


