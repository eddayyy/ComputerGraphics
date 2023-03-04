//
//  quaternion_T.h
//  linear_algebra
//
//  Created by William McCarthy on 172//21.
//

#ifndef quaternion_T_h
#define quaternion_T_h

#include <cmath>
#include "vector3d_T.h"
#include "matrix3d_T.h"

template <typename T> class quaternion;
template <typename T> using quat = class quaternion<T>;
typedef quat<double> quatD;

template <typename T>
class quaternion {
public:
  quaternion(T w_=T(), T x_=T(), T y_=T(), T z_=T())
  : w(w_), x(x_), y(y_), z(z_) { }

  static quaternion i(){ return quaternion(0.0, 1.0, 0.0, 0.0); }
  static quaternion j(){ return quaternion(0.0, 0.0, 1.0, 0.0); }
  static quaternion k(){ return quaternion(0.0, 0.0, 0.0, 1.0); }

  static double ii(){ return -1; }
  static double jj(){ return -1; }
  static double kk(){ return -1; }
  static double ijk(){ return -1; }

  static quaternion ij(){ return quaternion::k();}
  static quaternion jk(){ return quaternion::i();}
  static quaternion ki(){ return quaternion::j();}

  static quaternion ji(){ return -quaternion::k(); }
  static quaternion kj(){ return -quaternion::i();}
  static quaternion ik(){ return -quaternion::j();}

  friend quaternion operator+(const quaternion& a, const quaternion& b){
    if (typeid(T) == typeid(double)) {
        return quaternion(b.w, a.x + b.x, a.y + b.y, a.z + b.z);
    }
    return quaternion(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z);
  };

  friend quaternion operator-(const quaternion& a, const quaternion& b){ return quaternion(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z); }

  friend quaternion operator*(const quaternion& a, const quaternion& b){
    double w1 = a.w, x1 = a.x, y1 = a.y, z1 = a.z;
    double w2 = b.w, x2 = b.x, y2 = b.y, z2 = b.z;
    return quaternion(w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
                       w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                       w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
                       w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2);
  }

  friend quaternion operator+(const quaternion& q, T k){ return quaternion(q.w + k, q.x, q.y, q.z); }

  friend quaternion operator+(T k, const quaternion& q){ return quaternion(q.w + k, q.x, q.y, q.z); }


  friend quaternion operator-(const quaternion& q, T k){ return quaternion(q.w + k, q.x, q.y, q.z); }
  friend quaternion operator-(T k, const quaternion& q){ return q + k; }

  friend quaternion operator*(const quaternion& q, T k){ return quaternion<T>(q.w * k, q.x * k, q.y * k, q.z * k); }
  friend quaternion operator*(T k, const quaternion& q){ return quaternion<T>(k * q.w, k * q.x, k * q.y, k * q.z); }
  friend quaternion operator/(const quaternion& q, T k){ return quaternion(q.w / k, q.x / k, q.y / k, q.z / k); }


  quaternion operator-() const{ return quaternion<T>(-w, -x, -y, -z); }

  friend bool operator==(const quaternion& q, const quaternion& r){ return (q.w == r.w && q.x == r.x && q.y == r.y && q.z == r.z); }
  friend bool operator!=(const quaternion& q, const quaternion& r){ return !(q == r); }
  vector3d<T> vector() const{ return vector3d<T>("Vector", 3, {x, y, z}); }
  T scalar() const{ return w; }

  quaternion unit_scalar() const{
    T mag = std::sqrt(w * w + x * x + y * y + z * z);
    return quaternion<T>(w / mag, x / mag, y / mag, z / mag);
  }

  quaternion conjugate() const { return quaternion(w, -x, -y, -z); }

  quaternion inverse() const{
    double n = norm();
    if (n != 0) {
      return conjugate() / n;
    }
    return quaternion();
  }

  quaternion unit() const{
    double n = norm();
    if (n != 0) {
      return (*this) / n;
    }
    return quaternion();
  }

  double norm() const{ return w * w + x * x + y * y + z * z; }
  double magnitude() const{ return std::sqrt(norm()); }

  double dot(const quaternion& v) const{ return w * v.w + x * v.x + y * v.y + z * v.z; }

  double angle(const quaternion& v) const{
    double dp = dot(v);
    double mag = magnitude() * v.magnitude();
    if (mag != 0) {
      return std::acos(dp / mag);
    }
    return 0;
  }

  matrix3d<T> rot_matrix() const{
    matrix3d<T> m;
    T xx = x * x, xy = x * y, xz = x * z, xw = x * w;
    T yy = y * y, yz = y * z, yw = y * w;
    T zz = z * z, zw = z * w;

    m(0, 0) = 1 - 2 * (yy + zz);
    m(0, 1) = 2 * (xy - zw);
    m(0, 2) = 2 * (xz + yw);

    m(1, 0) = 2 * (xy + zw);
    m(1, 1) = 1 - 2 * (xx + zz);
    m(1, 2) = 2 * (yz - xw);

    m(2, 0) = 2 * (xz - yw);
    m(2, 1) = 2 * (yz + xw);
    m(2, 2) = 1 - 2 * (xx + yy);

    return m;
  }

 // rotates point pt (pt.x, pt.y, pt.z) about (axis.x, axis.y, axis.z) by theta
 static vec3 rotate(const vector3D& pt, const vector3D& axis, double theta){
     // normalize the axis vector
    vector3D axis_norm = axis.mag();

    // calculate sin and cos of theta/2
    double sintheta2 = sin(theta/2);
    double costheta2 = cos(theta/2);

    // create quaternion representing the rotation
    quaternion<T> q(costheta2, axis_norm.x() * sintheta2, axis_norm.y() * sintheta2, axis_norm.z() * sintheta2);

    // calculate the inverse of the quaternion
    quaternion<T> q_inv = q.inverse();

    // create a quaternion representing the point to be rotated
    quaternion<T> p(0, pt.x(), pt.y(), pt.z());

    // apply the rotation
    quaternion<T> q_result = q * p * q_inv;

    // return the rotated vector
    return vec3(q_result.x(), q_result.y(), q_result.z());
 }

 friend std::ostream& operator<<(std::ostream& os, const quaternion& q) {
   os << "Quat(";
   if (q ==  quaternion::i())  { return os <<  "i)"; }
   if (q == -quaternion::i())  { return os << "-i)"; }
   if (q ==  quaternion::j())  { return os <<  "j)"; }
   if (q == -quaternion::j())  { return os << "-j)"; }
   if (q ==  quaternion::k())  { return os <<  "k)"; }
   if (q == -quaternion::k())  { return os << "-k)"; }

   if (q.magnitude() == 0.0 && q.w == 0)   { return os << "0)"; }
   if (q.magnitude() == 0.0 && q.w == 0)   { return os << "0)"; }
   if (q.magnitude() == 1.0 && q.w == 1)   { return os << "1)"; }
   if (q.vector().mag() == 0.0)      { return os << q.w << ")"; }
   else { return os << q.w << q.vector() << ")"; }
 }

 static void run_tests();

private:
 T w, x, y, z;
};

void plane_rotation(const std::string& msg, const quatD& plane, const std::initializer_list<double>& li) {
 matrix3D rotate = matrix3D("rot_matrix", 3, li);
 assert(plane.rot_matrix() == rotate);
 std::cout << msg << " is: " << plane << plane.rot_matrix() << "\n";
}


std::string yes_or_no(bool condition) { return condition ? "YES" : "no"; }

template <typename T>
void quaternion<T>::run_tests() {
 quatD a(1, 2, 3, 4), b(4, 0, 0, 7), c(0, 1, 1, 0), d(0, 0, 1, 0);
 quatD e(0, 0, 0, 1), f(0, 0, 0, 0), g(1, 0, 0, 0), h(3, 0, 0, 0);
 std::cout << "a = " << a << ")\nb = " << b << ")\nc = " << c << ")\nd = " << d
           << ")\ne = " << e << ")\nf = " << f << ")\ng = " << g << ")\nh = " <<  h << "\n";

 std::cout << "c + d = " <<  c + d << "\nc + d + e = " << c + d + e << "\n";
 std::cout << "5 * h = " << 5 * h << "\nh * 5 = " << h * 5 << "\nh / 3.0 = " << h / 3.0 << "\n\n";

 std::cout << "h.magnitude() is " << h.magnitude() << "\nh.unit() is " << h.unit();
 std::cout << "g.unit() is " << g.unit() << "\na.unit() is " << a.unit() << ")\n\n";

 std::cout << "a.vector() is " << a.vector() << "\na.scalar() is " << a.scalar() << "\n";
 std::cout << "a.conjugate() is " << a.conjugate() << "\na.inverse() is " << a.inverse()
           << "\na * a.inverse() is " << a * a.inverse() << "\n\n";

 std::cout << "c == d is " << yes_or_no(c == d) << "\nc != d is " << yes_or_no(c != d);
 std::cout << "\ne == e is " << yes_or_no(e == e) << "\ne != e is " << yes_or_no(e != e) << "\n";

 std::cout << "\n\nquat.ij is: " << quatD::ij() << "\nquat.jk is: " << quatD::jk()
           << "\nquat.ki is: " << quatD::ki() << "\n";
 assert(quatD::ij() == quatD::k());
 assert(quatD::jk() == quatD::i());
 assert(quatD::ki() == quatD::j());

 std::cout << "\nquat.ji is: " << quatD::ji() << "\nquat.kj is: " << quatD::kj()
           << "\nquat.ik is: " << quatD::ik() << "\nquat.ijk is: " << quatD::ijk() << "\n";
 assert(quatD::ji() == -quatD::k());
 assert(quatD::kj() == -quatD::i());
 assert(quatD::ik() == -quatD::j());

 std::cout << "\nquat.ii is: " << quatD::ii() << "\nquat.jj is: " << quatD::jj()
           << "\nquat.kk is: " << quatD::kk() << "\n";
 assert(quatD::ii()  == -1);
 assert(quatD::jj()  == -1);
 assert(quatD::kk()  == -1);
 assert(quatD::ijk() == -1);

 std::cout << "\nangle (deg) between c and d is: " << c.angle(d) << "\n";
 quatD c_minus_d = c - d;
 std::cout << "c_minus_d is: " << c_minus_d;
 matrix3D rot_matrix = c_minus_d.rot_matrix();
 std::cout << "rot_matrix of c_minus_d is: " << c_minus_d.rot_matrix() << "\n";

 double rad2_2 = sqrt(2)/2.0;
 std::cout << "// -------------- LEVEL FLIGHT -------------------')\n";
 plane_rotation("levelFlight(E)", quatD(1),                     {  1,  0,  0,   0,  1,  0,   0,  0,  1 });
 plane_rotation("levelFlight(N)", quatD(rad2_2, 0, rad2_2,  0), {  0,  0,  1,   0,  1,  0,  -1,  0,  0 });
 plane_rotation("levelFlight(W)", quatD(0,      0,  1,      0), { -1,  0,  0,   0,  1,  0,   0,  0, -1 });
 plane_rotation("levelFlight(S)", quatD(rad2_2, 0, -rad2_2, 0), {  0,  0, -1,   0,  1,  0,   1,  0,  0} );
 std::cout << "LEVEL FLIGHT assertions passed ..............................................\n";
 std::cout << "// --------- end LEVEL FLIGHT ------------------------)\n";

 std::cout << "// -------------- STRAIGHT UP -------------------')\n";
 plane_rotation("straightUp(E)", quatD(rad2_2, 0, 0, rad2_2),   {  0, -1,  0,   1,  0,  0,   0,  0,  1 } );
 plane_rotation("straightUp(N)", quatD(0.5, 0.5, 0.5, 0.5),     {  0,  0,  1,   1,  0,  0,   0,  1,  0 } );
 plane_rotation("straightUp(W)", quatD(0, rad2_2, rad2_2, 0),   {  0,  1,  0,   1,  0,  0,   0,  0, -1 } );
 plane_rotation("straightUp(S)", quatD(0.5, -0.5, -0.5, 0.5),   {  0,  0, -1,   1,  0,  0,   0, -1,  0 } );
 std::cout << "STRAIGHT UP assertions passed..............................................\n";
 std::cout << "// --------- end STRAIGHT UP ------------------------)\n\n";


 std::cout << "// -------------- STRAIGHT DOWN ------------------')\n";
 plane_rotation("straightDown(E)", quatD(rad2_2, 0, 0, -rad2_2), {  0,  1,  0,  -1,  0,  0,   0,  0,  1 } );
 plane_rotation("straightDown(E)", quatD(0.5, -0.5, 0.5, -0.5),  {  0,  0,  1,  -1,  0,  0,   0, -1,  0 });
 plane_rotation("straightDown(E)", quatD(0, -rad2_2, rad2_2, 0), {  0, -1,  0,  -1,  0,  0,   0,  0, -1 } );
 plane_rotation("straightDown(E)", quatD(0.5, 0.5, -0.5, -0.5),  {  0,  0, -1,  -1,  0,  0,   0,  1,  0 });
 std::cout << "STRAIGHT DOWN assertions passed..............................................\n";
 std::cout << "// --------- end STRAIGHT DOWN ----------------------)\n\n";


 std::cout << "\n\n -------- BANK/ROLL ----------------\n";
 std::cout << "\nBanking/Rolling 90 degrees left...\n";
 plane_rotation("plane_E_bankLeft90", quatD(rad2_2, rad2_2, 0, 0),  {  1,  0,  0,   0,  0, -1,   0,  1,  0 } );
 plane_rotation("plane_N_bankLeft90", quatD(0.5, 0.5, 0.5, -0.5),   {  0,  1,  0,   0,  0, -1,  -1,  0,  0 } );
 plane_rotation("plane_W_bankLeft90", quatD(0, 0, rad2_2, -rad2_2), { -1,  0,  0,   0,  0, -1,   0, -1,  0 }  );
 plane_rotation("plane_W_bankLeft90", quatD(0.5, 0.5, -0.5, 0.5),   {  0, -1,  0,   0,  0, -1,   1,  0,  0 } );
 std::cout << "ROLL 90 deg left assertions passed..............................................\n";

 std::cout << "\n\nBanking/Rolling 180 degrees...\n";
 plane_rotation("plane_E_bankLeft180", quatD(0, 1, 0, 0),            {  1,  0,  0,   0, -1,  0,   0,  0, -1 });
 plane_rotation("plane_N_bankLeft180", quatD(0, rad2_2, 0, -rad2_2), {  0,  0, -1,   0, -1,  0,  -1,  0,  0 });
 plane_rotation("plane_W_bankLeft180", quatD(0, 0, 0, 1),            { -1,  0,  0,   0, -1,  0,   0,  0,  1 });
 plane_rotation("plane_S_bankLeft180", quatD(0, rad2_2, 0, rad2_2),  {  0,  0,  1,   0, -1,  0,   1,  0,  0 });
 std::cout << "ROLL 180 degrees assertions passed..............................................\n";

 std::cout << "\n\nBanking/Rolling 90 degrees right...\n";
 plane_rotation("plane_E_bankRight90", quatD(rad2_2, -rad2_2, 0, 0), {  1,  0,  0,   0,  0,  1,   0, -1,  0 });
 plane_rotation("plane_N_bankRight90", quatD(0.5, -0.5, 0.5, 0.5),   {  0, -1,  0,   0,  0,  1,  -1,  0,  0 });
 plane_rotation("plane_W_bankRight90", quatD(0, 0, rad2_2, rad2_2),  { -1,  0,  0,   0,  0,  1,   0,  1,  0 });
 plane_rotation("plane_S_bankRight90", quatD(0.5, -0.5, -0.5, -0.5), {  0,  1,  0,   0,  0,  1,   1,  0,  0 });
 std::cout << "ROLL 90 deg right assertions passed..............................................\n";
 std::cout << "\n -------- end BANK/ROLL ----------------\n";

 std::cout << "\nALL PLANE ROTATION ASSERTIONS PASSED ............................................\n\n";

 std::cout << "SEE THIS WEBSITE for DETAILED DIAGRAMS on the TESTS of the PLANE's rotations\n";
 std::cout << "https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/examples/index.htm\n";
}


#endif /* quaternion_T_h */

