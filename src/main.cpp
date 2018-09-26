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
	           {{1.f, 0.5f, 0.f}, {0.f, 1.f, -2.f}});

	cout << "---- TEST #3 ----\n";
	print_test(Vector3(0.f, 0.f, 0.f), Vector3(0.f, 1.f, 1.f), Vector3(0.f, 1.f, 0.f),
	           "The generated matrix should be a rotation matrix around x-axis for 45 degrees.",
	           {{0.f, 1.f, 0.f}, {1.f, 0.707f, 0.707f}});

	cout << "---- TEST #4 ----\n";
	print_test(Vector3(1.f, 1.f, 1.f), Vector3(0.f, 1.f, 1.f), Vector3(0.f, 1.f, 0.f),
	           "The generated matrix should be a combination of translation of (-1.0, -1.0, -1.0) and rotation matrix around x-axis for 45 degrees.",
	           {{0.f, 1.f, 0.f}, {1.f, 0.707f, 0.707f}});
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

	Mat4 matrix, inverse_matrix;
	coordinateSystemTransformationMatricesFromPositionNormalUp(VRP, VPN, VUP, matrix, inverse_matrix);

	cout << "The generated transformation matrix is:\n";
	cout << matrix << endl;

	cout << "The generated inverse transformation matrix is:\n";
	cout << inverse_matrix << endl;

	if (test_points.empty()) return;

	cout << "Testing Points:\n";
	for (Vector3& point : test_points)
		cout << point << " -> " << matrix * point << endl;

	cout << endl;
}
