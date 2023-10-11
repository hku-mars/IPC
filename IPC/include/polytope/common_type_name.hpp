#ifndef _UTILS_COMMON_TYPE_NAME_
#define _UTILS_COMMON_TYPE_NAME_

#include "Eigen/Dense"
#include "iostream"
#include "vector"
#include "fmt/color.h"
using namespace fmt;
using namespace std;

#include <std_msgs/ColorRGBA.h>

/*
 * @\brief Rename the float type used in lib

    Default is set to be double, but user can change it to float.
*/
typedef double decimal_t;

///Pre-allocated std::vector for Eigen using vec_E
template<typename T>
using vec_E = std::vector<T, Eigen::aligned_allocator<T>>;
///Eigen 1D float vector
template<int N>
using Vecf = Eigen::Matrix<decimal_t, N, 1>;
///Eigen 1D int vector
template<int N>
using Veci = Eigen::Matrix<int, N, 1>;
///MxN Eigen matrix
template<int M, int N>
using Matf = Eigen::Matrix<decimal_t, M, N>;
///MxN Eigen matrix with M unknown
template<int N>
using MatDNf = Eigen::Matrix<decimal_t, Eigen::Dynamic, N>;

///MxN Eigen matrix with N unknown
template<int M>
using MatMDf = Eigen::Matrix<decimal_t, M, Eigen::Dynamic>;

///Vector of Eigen 1D float vector
template<int N>
using vec_Vecf = vec_E<Vecf<N>>;
///Vector of Eigen 1D int vector
template<int N>
using vec_Veci = vec_E<Veci<N>>;

///Eigen 1D float vector of size 2
typedef Vecf<2> Vec2f;
///Eigen 1D int vector of size 2
typedef Veci<2> Vec2i;
///Eigen 1D float vector of size 3
typedef Vecf<3> Vec3f;
///Eigen 1D int vector of size 3
typedef Veci<3> Vec3i;
///Eigen 1D float vector of size 4
typedef Vecf<4> Vec4f;
///Column vector in float of size 6
typedef Vecf<6> Vec6f;

///Vector of type Vec2f.
typedef vec_E<Vec2f> vec_Vec2f;
///Vector of type Vec2i.
typedef vec_E<Vec2i> vec_Vec2i;
///Vector of type Vec3f.
typedef vec_E<Vec3f> vec_Vec3f;
///Vector of type Vec3i.
typedef vec_E<Vec3i> vec_Vec3i;

///2x2 Matrix in float
typedef Matf<2, 2> Mat2f;
///3x3 Matrix in float
typedef Matf<3, 3> Mat3f;
///4x4 Matrix in float
typedef Matf<4, 4> Mat4f;
///6x6 Matrix in float
typedef Matf<6, 6> Mat6f;

///Dynamic Nx1 Eigen float vector
typedef Vecf<Eigen::Dynamic> VecDf;
///Nx2 Eigen float matrix
typedef MatDNf<2> MatD2f;
///Nx3 Eigen float matrix
typedef MatDNf<3> MatD3f;
///Nx4 Eigen float matrix
typedef MatDNf<4> MatD4f;
///4xM Eigen float matrix
typedef MatMDf<4> Mat4Df;
typedef MatD4f MatPlanes;
///3xM Eigen float matrix
typedef MatMDf<3> Mat3Df;
typedef Mat3Df MatPoints;

///Dynamic MxN Eigen float matrix
typedef Matf<Eigen::Dynamic, Eigen::Dynamic> MatDf;

///Allias of Eigen::Affine2d
typedef Eigen::Transform<decimal_t, 2, Eigen::Affine> Aff2f;
///Allias of Eigen::Affine3d
typedef Eigen::Transform<decimal_t, 3, Eigen::Affine> Aff3f;
#endif

#ifndef EIGEN_QUAT
#define EIGEN_QUAT
///Allias of Eigen::Quaterniond
typedef Eigen::Quaternion<decimal_t> Quatf;
#endif

#ifndef EIGEN_EPSILON
#define EIGEN_EPSILON
///Compensate for numerical error
constexpr decimal_t epsilon_ = 1e-10; // numerical calculation error


/// For std color
class Color : public std_msgs::ColorRGBA {
 public:
  Color() : std_msgs::ColorRGBA() {}

  Color(double red, double green, double blue) : Color(red, green, blue, 1.0) {}

  Color(double red, double green, double blue, double alpha) : Color() {
	r = red;
	g = green;
	b = blue;
	a = alpha;
  }

  static const Color White() { return Color(1.0, 1.0, 1.0); }

  static const Color Black() { return Color(0.0, 0.0, 0.0); }

  static const Color Gray() { return Color(0.5, 0.5, 0.5); }

  static const Color Red() { return Color(1.0, 0.0, 0.0); }

  static const Color Green() { return Color(0.0, 0.96, 0.0); }

  static const Color Blue() { return Color(0.0, 0.0, 1.0); }

  static const Color SteelBlue() { return Color(0.4, 0.7, 1.0); }

  static const Color Yellow() { return Color(1.0, 1.0, 0.0); }

  static const Color Orange() { return Color(1.0, 0.5, 0.0); }

  static const Color Purple() { return Color(0.5, 0.0, 1.0); }

  static const Color Chartreuse() { return Color(0.5, 1.0, 0.0); }

  static const Color Teal() { return Color(0.0, 1.0, 1.0); }

  static const Color Pink() { return Color(1.0, 0.0, 0.5); }
};

typedef Eigen::Matrix<double, 3, 3> StatePVA;
typedef Eigen::Matrix<double, 3, 4> StatePVAJ;
typedef std::pair<double, Vec3f> TimePosPair;
typedef std::pair<Vec3f, Vec3f> Line;

struct QuadState {
  Vec3f position, velocity, acceleration, jerk;
  double yaw;
  double callback_time;
  bool rcv{false};
  Quatf q;
};

#endif