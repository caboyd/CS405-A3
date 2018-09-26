#pragma once
class Vector3;

void orthoNormalVectors(const Vector3& v1, const Vector3& v2, Vector3& u, Vector3& v, Vector3& n);
void coordinateSystemTransformationMatricesFromPositionNormalUp(const Vector3& VRP, const Vector3& VPN, const Vector3& VUP,
                                               Mat4& world_to_camera,
                                               Mat4& camera_to_world);

       