#include "Vector3.h"
#include "Mat4.h"
#include <iostream>

#include "Math.h"

using std::cout;


int main()
{
	Vector3 VRP(5.f,5.f,5.f), VPN(0.f,0.f,-1.f), VUP(0.f,1.f,0.1f); 
	Vector3 LRP(0.f,10.f,0.f), LPN(-1.f,-1.f,0.f), LUP(0.1f,1.f,0.f); 

	Vector3 n, u, v;

	orthoNormalVectors(VPN, VUP, u, v, n);

	cout << "forward: " << n << std::endl;
	cout << "up: " << v << std::endl;
	cout << "right: " << u << std::endl;

	Mat4 world_to_camera, camera_to_world;
	worldToCameraMatricesFromPositionNormalUp({5,5,5},{0,0,-1},{0,1,0},world_to_camera,camera_to_world);

	Mat4 world_to_light, light_to_world;
	WorldToLightMatricesFromPositionNormalUp({0,10,0},{-1,-1,0},{0,1,0},world_to_light,light_to_world);

	Mat4 camera_to_light, light_to_camera;

	camera_to_light = world_to_light * camera_to_world;
	light_to_camera = world_to_camera * light_to_world;

	cout << std::endl << world_to_camera ;
	cout << std::endl << camera_to_world ;
	cout << std::endl << camera_to_world * world_to_camera;

	cout << std::endl << world_to_light ;
	cout << std::endl << light_to_world ;
	cout << std::endl << world_to_light * light_to_world;

	cout << std::endl << camera_to_light;
	cout << std::endl << light_to_camera;
	cout << std::endl << light_to_camera * camera_to_light;
}
