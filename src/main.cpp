
#include "Vector3.h"
#include "Mat4.h"
#include <iostream>

#include "Math.h"

using std::cout;


int main()
{
	Vector3 forward(1.f,0.f,0.f);
	Vector3 weird_up(1.0f,0.6f,0.1f);

	Vector3 n,u,v;

	orthoNormalVectors(forward,weird_up,u,v,n);

	cout << "forward: " << n << std::endl;
	cout << "up: " << v << std::endl;
	cout << "right: " << u << std::endl;

	Mat4 m = Mat4::fromIdentity();
	m = Mat4::fromXRotation(3.14159f/2.f);

	cout << std::endl << std::right << m;

}
