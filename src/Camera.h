#pragma once
#include "Ray.h"
#include "Vector3.h"
#include "Mat4.h"

class Camera
{
public:
	Vector3 position;
	Vector3 normal;

	Vector3 screen_horizontal;
	Vector3 screen_vertical;
	Vector3 lower_left_corner;

	Camera(const Vector3& pos, const Vector3& normal, const Vector3& up, float vFov, float aspect_ratio);
	Camera(const Vector3& pos, const Vector3& normal, const Vector3& up, float focal_length, float film_size, float aspect_ratio);

	Ray getRay(float x, float y) const;

	Mat4 getViewMatrix() const;
	Mat4 getViewMatrixInverse() const;
};
