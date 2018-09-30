#include "Vector3.h"
#include "Mat4.h"
#include <iostream>
#include <string>
#include <vector>

#include "Math.h"
using std::cout;
using std::endl;

/**
 *	Creates transformation and inverse transformation matrices and prints them
 *	Uses the given test points and uses the transformation matrix to transform them to the new coordinate system.
 *	Prints the test input point and output point after is has been transformed.
 *
 */
void print_test(const Vector3& VRP, const Vector3& VPN, const Vector3& VUP, std::string message,
                std::vector<Vector3> test_points);
void print_test_camera_to_light(const Vector3& VRP, const Vector3& VPN, const Vector3& VUP,
                                const Vector3& LRP, const Vector3& LPN, const Vector3& LUP,
                                std::string message, std::vector<Vector3> test_points);

int main()
{
	Vector3 VRP(5.f, 5.f, 5.f), VPN(0.f, 0.f, -1.f), VUP(0.f, 1.f, 0.1f);

	cout << "---- TEST #1 ----\n";
	print_test(Vector3(0.f, 0.f, 0.f), Vector3(0.f, 0.f, 1.f), Vector3(0.f, 1.f, 0.f),
	           "The generated matrix should be an identity matrix",
	           {});

	cout << "---- TEST #2 ----\n";
	print_test(Vector3(0.f, 0.f, 0.f), Vector3(1.f, 0.f, 0.f), Vector3(0.f, 1.f, 0.f),
	           "The generated matrix should be a rotation matrix around y-axis for -90 degrees.",
	           {{1.f, 0.5f, 0.f}, {0.f, 1.f, -2.f}, {0.f, .707f, .707f}});

	cout << "---- TEST #3 ----\n";
	print_test(Vector3(0.f, 0.f, 0.f), Vector3(0.f, 1.f, 1.f), Vector3(0.f, 1.f, 0.f),
	           "The generated matrix should be a rotation matrix around x-axis for 45 degrees.",
	           {{0.f, 1.f, 0.f}, {1.f, 0.707f, 0.707f}});

	cout << "---- TEST #4 ----\n";
	print_test(Vector3(1.f, 1.f, 1.f), Vector3(0.f, 1.f, 1.f), Vector3(0.f, 1.f, 0.f),
	           "The generated matrix should be a combination of translation of (-1.0, -1.0, -1.0) and rotation matrix around x-axis for 45 degrees.",
	           {{0.f, 1.f, 0.f}, {1.f, 0.707f, 0.707f}});

	cout << "---- TEST #6 ----\n";
	print_test_camera_to_light(Vector3(0.f, 2.f, 0.f), Vector3(0.f, 0.f, 1.f), Vector3(0.f, 1.f, 0.f),
	                           Vector3(0.f, -2.f, 0.f), Vector3(0.f, 0.f, 1.f), Vector3(0.f, 1.f, 0.f),
	                           "The generated matrix should be a translation of (0.0, -4.0, 0.0)",
	                           {{0.f, 1.f, 0.f}, {1.f, 2.f, 1.f}});

	cout << "---- TEST #7 ----\n";
	print_test_camera_to_light(Vector3(0.f, 0.f, 0.f), Vector3(0.f, 1.f, 1.f), Vector3(0.f, 1.f, 0.f),
	                           Vector3(0.f, 0.f, 0.f), Vector3(1.f, 0.f, 0.f), Vector3(0.f, 1.f, 0.f),
	                           "The generated matrix should be a rotation around x-axis for 45 degrees and the y-axis for 90 degrees",
	                           {{0.f, .707107f, .707107f}, {1.f, 0.707107f, 0.707107f}});
}

void print_test(const Vector3& VRP, const Vector3& VPN, const Vector3& VUP, std::string message,
                std::vector<Vector3> test_points)
{
	cout << message << endl;
	cout << "VRP = " << VRP << endl;
	cout << "VPN = " << VPN << endl;
	cout << "VUP = " << VUP << endl << endl;

	Vector3 n, u, v;
	orthoNormalVectors(VPN, VUP, u, v, n);

	// cout << "n_z_forward: " << n << std::endl;
	// cout << "v_y_up: " << v << std::endl;
	// cout << "u_x_right: " << u << std::endl;
	// cout << "n dot v: " << n.dot(v) << std::endl;

	Mat4 matrix, inverse_matrix;
	coordinateSystemTransformationMatricesFromPositionNormalUp(VRP, VPN, VUP, matrix, inverse_matrix);

	cout << "The generated transformation matrix is:\n";
	cout << matrix << endl;

	cout << "The generated inverse transformation matrix is:\n";
	cout << inverse_matrix << endl;

	if (test_points.empty()) return;

	cout << "Testing Points:\n";
	cout << "These points will be transformed using the transformation matrix\n";
	for (Vector3& point : test_points)
		cout << point << " -> " << matrix * point << endl;

	cout << endl;
}

void print_test_camera_to_light(const Vector3& VRP, const Vector3& VPN, const Vector3& VUP,
                                const Vector3& LRP, const Vector3& LPN, const Vector3& LUP,
                                std::string message, std::vector<Vector3> test_points)
{
	cout << "Testing camera to world to light transformation\n";
	cout << message << endl;

	cout << "Camera:" << endl;
	cout << "VRP = " << VRP << endl;
	cout << "VPN = " << VPN << endl;
	cout << "VUP = " << VUP << endl << endl;

	cout << "Light:" << endl;
	cout << "LRP = " << LRP << endl;
	cout << "LPN = " << LPN << endl;
	cout << "LUP = " << LUP << endl << endl;

	Vector3 n, u, v;
	orthoNormalVectors(VPN, VUP, u, v, n);
	Mat4 camera_to_world, world_to_camera;
	coordinateSystemTransformationMatricesFromPositionNormalUp(VRP, VPN, VUP, world_to_camera, camera_to_world);

	orthoNormalVectors(LPN, LUP, u, v, n);
	Mat4 light_to_world, world_to_light;
	coordinateSystemTransformationMatricesFromPositionNormalUp(LRP, LPN, LUP, world_to_light, light_to_world);

	Mat4 camera_to_light = world_to_light * camera_to_world;
	Mat4 light_to_camera = world_to_camera * light_to_world;
	cout << "The generated camera_to_light transformation matrix is:\n";
	cout << camera_to_light << endl;

	cout << "The generated light_to_camera transformation matrix is:\n";
	cout << light_to_camera << endl;

	if (test_points.empty()) return;

	cout << "Testing Points:\n";
	cout << "These points will be transformed using the camera_to_light matrix\n";
	for (Vector3& point : test_points)
		cout << point << " -> " << camera_to_light * point << endl;

	cout << endl;
}
