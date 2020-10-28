#pragma once

#include <cmath>

/**
 * A class representing a 3-component vector (x, y, z) of real numbers
 */
class Vec3D
{
private:
	double m_x;
	double m_y;
	double m_z;

public:
	/**
	 * The default constructor of the class, which initializes all components to 0
	 */
	Vec3D();

	/**
	 * A constructor of the class allowing to choose the initialization values of all the components
	 *
	 * @param x The value of the first component
	 *
	 * @param y The value of the second component
	 *
	 * @param z The value of the third component
	 */
	Vec3D(double x, double y, double z);

	/**
	 * The destructor of the class
	 */
	~Vec3D();

	/**
	 * A method to read the value of the first component of the instance
	 *
	 * @returns The value of the first component
	 */
	double get_x() const;

	/**
	 * A method to read the value of the second component of the instance
	 *
	 * @returns The value of the second component
	 */
	double get_y() const;

	/**
	 * A method to read the value of the third component of the instance
	 *
	 * @returns The value of the third component
	 */
	double get_z() const;

	/**
	 * A method to set the value of the first component of the instance
	 *
	 * @param x The new value of the first component
	 *
	 * @warning The previous value will be lost when calling this method
	 */
	void set_x(double x);
	
	/**
	 * A method to set the value of the second component of the instance
	 *
	 * @param y The new value of the second component
	 *
	 * @warning The previous value will be lost when calling this method
	 */
	void set_y(double y);

	/**
	 * A method to set the value of the third component of the instance
	 *
	 * @param z The new value of the third component
	 *
	 * @warning The previous value will be lost when calling this method
	 */
	void set_z(double z);

	/**
	 * A method to set the values of all 3 components of the instance
	 *
	 * @param x The new value of the first component
	 *
	 * @param y The new value of the second component
	 *
	 * @param z The new value of the third component
	 *
	 * @warning The previous values will be lost when calling this method
	 */
	void set(double x, double y, double z);

	/**
	 * A method to compute the square of the euclidian norm of the vector
	 *
	 * @returns The squared norm of the instance
	 */
	double compute_squared_norm() const;

	/**
	 * A method to compute the euclidian norm of the vector
	 *
	 * @returns The norm of the instance
	 *
	 * @warning Computing a squared norm with this method is much heavier than with the compute_squared_norm() method
	 */
	double compute_norm() const;

	/**
	 * A method transforming the vector into a unitary normed one with the same orientation
	 *
	 * @warning The values of all the components will be irreversibly modified when calling this method
	 */
	void normalize();

	/**
	 * A boolean operator checking the equality of two vectors
	 *
	 * @returns true if all the components are equal, false otherwise
	 */
	bool operator==(Vec3D const& v);

	/**
	 * A boolean operator ckecking the non-equality of two vectors
	 *
	 * @returns true if at least one component is different, false otherwise
	 */
	bool operator!=(Vec3D const& v);

	/**
	 * A summation operator computing the result of the addition with another vector
	 *
	 * @param v The other vector
	 *
	 * @returns The sum of the two vectors
	 */
	Vec3D operator+(Vec3D const& v) const;
	
	/**
	 * A difference operator computing the result of the difference with another vector
	 *
	 * @param v The other vector
	 *
	 * @returns The difference of the two vectors
	 */
	Vec3D operator-(Vec3D const& v) const;

	/**
	 * An external product operator compution de result of the multiplication by a scalar
	 *
	 * @param a The scalar
	 *
	 * @returns The product of the vector by the scalar
	 *
	 * @warning The product is only defined by the right
	 */
	Vec3D operator*(double a) const;
	
	/**
	 * An external division operator compution de result of the division by a scalar
	 *
	 * @param a The scalar
	 *
	 * @returns The division of the vector by the scalar
	 *
	 * @warning The division is only defined by the right
	 */
	Vec3D operator/(double a) const;

	/**
	 * An incrementation operator based on the summation operator
	 *
	 * @param v The vector to add to the instance
	 *
	 * @warning Using this operator changes all the components of the instance
	 */
	void operator+=(Vec3D const& v);

	/**
	 * A decrementation operator based on the difference operator
	 *
	 * @param v The vector to substract to the instance
	 *
	 * @warning Using this operator changes all the components of the instance
	 */
	void operator-=(Vec3D const& v);

	/**
	 * A homotecy operator based on the external product operator
	 *
	 * @param a The scalar to multiply the instance by
	 *
	 * @warning Using this operator changes all the components of the instance
	 */
	void operator*=(double a);
	
	/**
	 * A homotecy operator based on the external division operator
	 *
	 * @param a The scalar to divide the instance by
	 *
	 * @warning Using this operator changes all the components of the instance
	 */
	void operator/=(double a);

	/**
	 * A static method computing the dot product of two vectors
	 *
	 * @param v1 The first vector
	 *
	 * @param v2 The second vector
	 *
	 * @returns The dot product v1.v2
	 */
	static double dot_product(Vec3D const& v1, Vec3D const& v2);

	/**
	 * A static method computing the cross product of two vectors
	 *
	 * @param v1 The first vector
	 *
	 * @param v2 The second vector
	 *
	 * @returns The cross product v1^v2
	 */
	static Vec3D cross_product(Vec3D const& v1, Vec3D const& v2);
};