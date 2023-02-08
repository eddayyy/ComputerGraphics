//============================================================
// FILE: vector3dT.h
//============================================================

#ifndef __vector3dT_H__
#define __vector3dT_H__

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <initializer_list>

template <typename T> class vector3d;

typedef vector3d<double> vector3D;
typedef vector3d<float>  vector3F;
typedef vector3d<int>    vector3I;
typedef vector3d<long>   vector3L;
typedef vector3D vec3;

template <typename T>
class vector3d {  // class that serves as both vectors and points
public:
  vector3d();
  vector3d(const std::string& name, int dims);
  vector3d(const std::string& name, int dims, const std::initializer_list<T>& li);
//---------------------------------------------------------------------
  T  operator[](int i) const;
  T& operator[](int i);
//---------------------------------------------------------------------
  std::string name() const;
  void name(const std::string& name);

  T x() const { return data_[0]; }   // read-only value of x
  T y() const  { return data_[1]; }
  T z() const  { return data_[2]; }

  T& x() { return data_[0]; }        // read-write value of x
  T& y() { return data_[1]; }
  T& z() { return data_[2]; }

//---------------------------------------
  vector3d<T>& operator+=(const vector3d<T>& v);
  vector3d<T>& operator-=(const vector3d<T>& v);
//---------------------------------------------------------------------
  vector3d<T>& operator+=(T k);
  vector3d<T>& operator-=(T k);
  vector3d<T>& operator*=(T k);
  vector3d<T>& operator/=(T k);
//---------------------------------------------------------------------
  vector3d<T> operator-();
//---------------------------------------------------------------------
  friend vector3d<T> operator+(const vector3d<T>& u, const vector3d<T>& v) {
    check_equal_dims(u, v);
    return vector3d<T>(u.name_ + "+" + v.name_, u.dims_,
                    { u[0] + v[0], u[1] + v[1], u[2] + v[2], 0} );
  }  

  friend vector3d<T> operator-(const vector3d<T>& u, const vector3d<T>& v) { /* TODO */ }
//---------------------------------------------------------------------
  friend vector3d<T> operator+(T k, const vector3d<T>& v) {
    return vector3d<T>(std::to_string(k) + "+" + v.name_, v.dims_, { k + v[0], k + v[1], k + v[2], 0 });
  }
  friend vector3d<T> operator+(const vector3d<T>& v, T k) { /* TODO */ }
//---------------------------------------------------------------------
  friend vector3d<T> operator-(T k, const vector3d<T>& v) { /* TODO */ }
  friend vector3d<T> operator-(const vector3d<T>& v, T k) { /* TODO */ }
  //---------------------------------------------------------------------
  friend vector3d<T> operator*(T k, const vector3d<T>& v) {
    return vector3d<T>(std::to_string(k) + v.name_, v.dims_, { k * v[0], k * v[1], k * v[2], 0 });
  }
  friend vector3d<T> operator*(const vector3d<T>& v, T k) { /* TODO */ };
  //---------------------------------------------------------------------
  friend vector3d<T> operator/(const vector3d<T>& v, T k) {
    if (k == 0) { throw new std::invalid_argument("divide by zero"); }
    double kinv = 1.0 / k;
    return kinv * v;
  }
//---------------------------------------------------------------------
  friend bool operator<(const vector3d<T>& u, const vector3d<T>& v) {
    check_equal_dims(u, v);
    return u.mag() < v.mag();
  }
  friend bool operator==(const vector3d<T>& u, const vector3d<T>& v) {
    check_equal_dims(u, v);
    return /* u.name_ == v.name_ && */ u[0] == v[0] && u[1] == v[1] && u[2] == v[2];
  }
  friend bool operator!=(const vector3d<T>& u, const vector3d<T>& v) { return !(u == v); }
//---------------------------------------------------------------------
  double dot(const vector3d<T>& other) const;
  double mag() const;
  double norm() const { return mag(); }  // L2 norm
  double angle(const vector3d<T>& other) const;
  vector3d<T> cross(const vector3d<T>& other) const;
//---------------------------------------------------------------------
  static vector3d<T> zero();
  static double value(double val) { return abs(val) < 1e-5 ? 0 : val; }

//---------------------------------------------------------------------
friend std::ostream& operator<<(std::ostream& os, const vector3d<T>& v) {
  os << "<" << v.name_ << ", ";
  if (v.dims_ == 0) { os << "empty>"; }
  else {
    for (int i = 0; i < v.dims_ + 1; ++i) {
      os << std::setw(3) << std::setprecision(3) << value(v[i]);
      if (i < v.dims_) { os << " "; }
    }
    os << ">";
  }
  return os;
}

void show() { std::cout << *this << "\n"; }
void show(const std::string& msg) { std::cout << msg << *this << "\n"; }

static void run_tests() {
  std::cout << "\n====================  TESTING VECTORS  ========================" << "\n";
  vector3d<double> u("u", 3, {1,  2,  4});
  // std::cout << u.name() << "\n";
  // std::cout << u << "\n";
  std::cout << "u.name_ is: " << u.name() << "\n";
  u.zero();
  u.show();
  vector3D v("v", 3, {8, 16, 32});
  vector3D i("i", 3, {1, 0, 0}), j("j", 3, {0, 1, 0}), k("k", 3, {0, 0, 1});
  vector3D w(3 * i + 4 * j - 2 * k);

  u.show();
  v.show();
  i.show();
  j.show();
  k.show();
  std::cout << "j + k is: " << j + k << "\n";
  w.show();

  std::cout << "*** asserting u == u and u != v" << "\n";
  assert(u == u);
  assert(u != v);
  
  std::cout << "*** asserting u + v == v + u   and  u - v == -(v - u)    and   -(-u) == u" << "\n";
  assert(u + v == v + u);
  assert(u - v == -(v - u));
  assert(-(-u) == u);
  
  std::cout << "*** 3.0 + u == u + 3.0   and   3.0 * u == u * 3.0" << "\n";
  assert(3.0 + u == u + 3.0);
  assert(3.0 * u == u * 3.0);

  // vector3d<double> a = u - 3.0;       // example of tracking down a failed assertion
  // vector3d<double> b = -( 3.0 - u);
  // u.show();
  // a.show();
  // b.show();
  std::cout << "*** asserting u - 3.- == -(3.0 - u)" << "\n";
  assert((u - 3.0) == -(3.0 - u));
  std::cout << "*** asserting 5.0 * u == u * 5.0" << "\n";
  assert((5.0 * u) / 5.0 == u);

  std::cout << "*** asserting u + vector3D::zero() == u" << "\n";
  assert(u + vector3D::zero() == u);

  std::cout << "*** asserting i.dot(j) == j.dot(k) == k.dot(i) == 0" << "\n";
  assert(i.dot(j) == j.dot(k) == k.dot(i) == 0);

  std::cout << "*** asserting i.cross(j) == k   and  j.cross(k) == i   and   k.cross(i) == j" << "\n";
  assert(i.cross(j) == k);
  assert(j.cross(k) == i);
  assert(k.cross(i) == j);
  
  std::cout << "*** asserting u.cross(v) == -v.cross(u)" << "\n";
  assert(u.cross(v) == -v.cross(u));
  assert(u.cross(v + w) == u.cross(v) + u.cross(w));
  assert((u.cross(v)).dot(u) == 0);

  std::cout << "i.angle(j) is: " << i.angle(j) << "\n";
  std::cout << "pi/2 is: " << M_PI/2 << "\n";

  std::cout << "*** asserting i.angle(j) == j.angle(k) == k.angle(i) == M_PI/2" << "\n";
  assert(i.angle(j) == M_PI_2);
  assert(j.angle(k) == M_PI_2);
  assert(k.angle(i) == M_PI_2);
  
  vector3D uhat = u / u.mag();
  u.show();
  uhat.show();
  std::cout << "length of uhat.mag() is... " << uhat.mag() << "\n";
  std::cout << "*** asserting u.hat.mag() - 1.0 < 1.0e-10" << "\n";
  assert(uhat.mag() - 1.0 < 1.0e-10);

  std::cout << "...test vectors assertions passed" << "\n";
  std::cout << "====================  FINISHED testing vectors  ========================\n" << "\n";
}

private:
  friend void check_equal_dims(const vector3d& u, const vector3d& v) {
    if (u.dims_ != v.dims_) { throw new std::invalid_argument("vector3d dims mismatch"); }
  }
  void check_bounds(int i) const;

private:
  std::string name_;
  int dims_;
  T data_[4];
};

//==============================================================================

template <typename T> vector3d<T>::vector3d() : vector3d("no_name", 3) {}  // 3d default dims
template <typename T> vector3d<T>::vector3d(const std::string& name, int dims)
  : name_(name), dims_(dims) {
    memset(data_, 0, dims_ * sizeof(double));
    data_[3] = 0;  // vectors have 0 at end, pts have 1
  }
template <typename T> vector3d<T>::vector3d(const std::string& name, int dims,
                                            const std::initializer_list<T>& li)
  : vector3d(name, dims) {
    int i = 0;
    for (T value : li) {
      if (i > dims_) { break; }
      data_[i++] = value;
    }
    data_[3] = 0;
  }

template <typename T> std::string vector3d<T>::name() const { return name_; }
template <typename T> void vector3d<T>::name(const std::string& name) { name_ = name; }

template <typename T>
vector3d<T>& vector3d<T>::operator+=(const vector3d<T>& v) {
  vector3d<T>& u = *this;
  u[0] += v[0];  u[1] += v[1];  u[2] += v[2];  u[3] = v[3];
  return *this;
}
template <typename T>
vector3d<T>& vector3d<T>::operator-=(const vector3d<T>& v) { /* TODO */ }
//---------------------------------------------------------------------
template <typename T> vector3d<T>& vector3d<T>::operator+=(T k) { /* TODO */ }
template <typename T> vector3d<T>& vector3d<T>::operator-=(T k) { /* TODO */ }
template <typename T> vector3d<T>& vector3d<T>::operator*=(T k) { /* TODO */ }
template <typename T> vector3d<T>& vector3d<T>::operator/=(T k) { /* TODO */ }

//---------------------------------------------------------------------
template <typename T>  /* read only idx */
T  vector3d<T>::operator[](int i) const {  check_bounds(i);  return data_[i]; }

template <typename T> T& vector3d<T>::operator[](int i) { /* TODO */ } // rw idx

//-----------------------
template <typename T>
vector3d<T> vector3d<T>::operator-() { return vector3d("-" + name_, dims_, { -data_[0], -data_[1], -data_[2], 0 }); }
//-----------------------
template <typename T>
double vector3d<T>::dot(const vector3d<T>& v) const { /* TODO */ }
//-----------------------
template <typename T>
double vector3d<T>::mag() const {  return sqrt(dot(*this));  }

template <typename T>
double vector3d<T>::angle(const vector3d<T>& v) const { /* TODO */ }
//-----------------------
template <typename T>
vector3d<T> vector3d<T>::cross(const vector3d<T>& v) const {
  check_equal_dims(*this, v);
  if (v.dims_ != 3) {
    throw new std::invalid_argument("cross_product only implemented for vector3d's");
  }
  return vector3d<T>(name_ + " x " + v.name_, dims_,
                    { data_[1]*v[2] - data_[2]*v[1], -(data_[0]*v[2] - data_[2]*v[0]),
                      data_[0]*v[1] - data_[1]*v[0], 0 });
}
//-----------------------
template <typename T>
vector3d<T> vector3d<T>::zero() { return vector3d<T>("zero", 3, {0, 0, 0, 0}); }
//-----------------------

template <typename T>
void vector3d<T>::check_bounds(int i) const {  // 1 extra dimension for pts/vectors
  if (i > dims_) { throw new std::invalid_argument("out of bounds"); }
}


#endif
