#pragma once
class Vector3;

/**
 * Given a direction vector v1, and non-collinear direction v2.
 * This function creates 3 orthonormal vectors where 
 * n (+Z) is the normalized direction of v1 (v1 normalized),
 * u (+X) is the normalized orthogonal vector to v2 and n (v2 cross n).
 * v (+Y) is the normalized up vector calculated using the direction of v2 (n cross u)
 */
void orthoNormalVectors(const Vector3& v1, const Vector3& v2, Vector3& u, Vector3& v, Vector3& n);


/**
 * Given 3 vectors, a position, normal direction, and up vector
 * This function calculates new orthonormal vectors to produce a transformation matrix and inverse transformation matrix
 * transformation_matrix = R * T
 * where R is the rotation matrix created from the orthonormal vectors
 * and T is the translation matrix created from translating from VRP to the origin.
 * 
 * inverse_transformation_matrix = T_inverse * R_inverse
 * where R_inverse is the inverse rotation matrix created from the orthonormal vectors
 * and T_inverse is the inverse translation matrix created from the translating from the origin to VRP.
 */
void coordinateSystemTransformationMatricesFromPositionNormalUp(const Vector3& VRP, const Vector3& VPN, const Vector3& VUP,
                                               Mat4& transformation_matrix,
                                               Mat4& inverse_transformation_matrix);

       