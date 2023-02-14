//============================================================
// FILE: matrix3dT.h
//============================================================
#ifndef __matrix_3d__
#define __matrix_3d__

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>

template <typename T> class matrix3d;
typedef matrix3d<double> matrix3D;
typedef matrix3d<float> matrix3F;
typedef matrix3d<int> matrix3I;
typedef matrix3d<long> matrix3L;
typedef matrix3D mat3;

template <typename T>
class matrix3d {
public:
  matrix3d();
  matrix3d(const std::string& name, int dims);
  matrix3d(const std::string& name, int dims, const std::initializer_list<vector3d<T>>& li);
  matrix3d(const std::string& name, int dims, const std::initializer_list<T>& li);
//=======================================================================
  matrix3d<T>& operator=(T array[9]);
  matrix3d<T>& operator=(T k);
//=======================================================================
// indexing ops...
  vector3d<T> operator[](int i) const;
  vector3d<T>& operator[](int i);
  
  T  operator()(int row, int col) const;
  T& operator()(int row, int col);
  
  T* opengl_memory(int row, int col);
//=======================================================================
  void name(const std::string& name);
  const std::string& name() const;

//============================ LINEAR ALGEBRA =========================
  matrix3d<T>& operator+=(T k);
  matrix3d<T>& operator-=(T k);
  matrix3d<T>& operator*=(T k);
  matrix3d<T>& operator/=(T k);
//=======================================================================
  matrix3d<T>& operator+=(const matrix3d<T>& b);
  matrix3d<T>& operator-=(const matrix3d<T>& b);
//=======================================================================
  matrix3d<T> operator-();
  matrix3d<T> operator+(const matrix3d<T>& b);
  matrix3d<T> operator-(const matrix3d<T>& b);
//=======================================================================
  friend matrix3d operator+(const matrix3d& a, T k) {
    return matrix3d(std::to_string(k) + "+" + a.name(), 3,
                    { a[0] + k, a[1] + k, a[2] + k });
  }
  friend matrix3d operator+(T k, const matrix3d& a) { return  a + k; }
  // Todo 
  friend matrix3d operator-(const matrix3d& a, T k) { 
    matrix3d<T> result(a.rows(), a.cols(), a.depth());
    for (int i = 0; i < a.rows(); ++i) {
      for (int j = 0; j < a.cols(); ++j) {
        for (int k = 0; k < a.depth(); ++k) {
          result(i, j, k) = a(i, j, k) - k;
        }
      }
    }
  return result;
  }
  // Todo
  friend matrix3d operator-(T k, const matrix3d& a) { 
    matrix3d<T> result(a.rows_, a.cols_, a.depth_);
    for (int i = 0; i < a.rows_; i++) {
      for (int j = 0; j < a.cols_; j++) {
        for (int k = 0; k < a.depth_; k++) {
          result(i, j, k) = k - a(i, j, k);
        }
      }
    }
    return result;
  }
  // Todo
  friend matrix3d operator*(const matrix3d& a, T k) { 
    matrix3d<T> result(a);
    for (auto& element : result.elements_) {
      element *= k;
    }
    return result;
   }
  // Todo
  friend matrix3d<T> operator*(T k, const matrix3d& a) { 
  matrix3d<T> result(a.get_rows(), a.get_cols(), T(0));
  for (int i = 0; i < a.get_rows(); i++) {
    for (int j = 0; j < a.get_cols(); j++) {
      result(i, j) = k * a(i, j);
    }
  }
  return result;
  }
  friend matrix3d operator/(const matrix3d& a, T k) { 
    if (abs(k) < epsilon_) { throw new std::invalid_argument("divide by zero error"); }
    return a * (1.0 / k);  
  }
//=======================================================================
  friend matrix3d operator*(const matrix3d& m, const vector3d<T>& v) { 
    matrix3d<T> result(m.name_ + " * " + v.name(), 3, { m(0,0) * v[0] + m(0,1) * v[1] + m(0,2) * v[2],
                                                        m(1,0) * v[0] + m(1,1) * v[1] + m(1,2) * v[2],
                                                        m(2,0) * v[0] + m(2,1) * v[1] + m(2,2) * v[2]});
    return result;
  }
  friend matrix3d operator*(const vector3d<T>& v, const matrix3d& m) { 
    matrix3d<T> result(m.name_ + " * " + v.name(), 3, { m(0,0) * v[0] + m(0,1) * v[1] + m(0,2) * v[2],
                                                        m(1,0) * v[0] + m(1,1) * v[1] + m(1,2) * v[2],
                                                        m(2,0) * v[0] + m(2,1) * v[1] + m(2,2) * v[2]});
    return result;
  }
  matrix3d<T> operator*(const matrix3d<T>& b);
//=======================================================================
  matrix3d<T> transpose() const;// create a new matrix transpose()
  T determinant() const;
  T trace() const;
//=======================================================================
  matrix3d<T> minors() const;     // see defn
  matrix3d<T> cofactor() const;   // (-1)^(i+j)*minors()(i, j)
  matrix3d<T> adjoint() const;    // cofactor.transpose()
  matrix3d<T> inverse() const;    // adjoint()/determinant()
//=======================================================================
  static matrix3d<T> identity(int dims);  // identity matrix
  static matrix3d<T> zero(int dims);      // zero matrix
//=======================================================================
  bool operator==(const matrix3d<T>& b) const;
  bool operator!=(const matrix3d<T>& b) const;
//=======================================================================

void show() { std::cout << *this << "\n"; }

friend std::ostream& operator<<(std::ostream& os, const matrix3d<T>& m) {
  os << "<'" << m.name_ << "', ";
  for (int i = 0; i < 3; ++i) {
    os << m.cols_[i];
  }
  os << "> OR by rows...\n";
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      os << std::setw(3) << m(i, j) << " ";
    }
    os << "\n";
  }
  return os << ">";
}

static void run_tests() {
  std::cout << "\n====================  TESTING MATRICES  ====================================================\n";
  matrix3D a("a", 3, {3, 2, 0,   0, 0, 1,   2, -2, 1});
  a.show();

  matrix3D b("b", 3, {1, 0, 5,   2, 1, 6,   3,  4, 0});
  b.show();

  matrix3D id = matrix3D::identity(3);
  id.show();

  assert(a * id == a);
  (a * id).show();

  // assert(a * b != -b * a);
  (a * b).show();
  (b * a).show();
  (-b * a).show();

  assert(a == a);
  assert(b == b);
  assert(a * b == a * b);
  assert(a * b != -b * a);

  std::cout << "so far so good\n";

  matrix3D aT = a.transpose();
  aT.show();

  matrix3D bT = b.transpose();
  bT.show();

  matrix3D abT = (a * b).transpose();
  abT.show();


  assert((a * b).transpose() == b.transpose() * a.transpose());

  std::cout << "test copy constructor -- program will crash at end if this is incorrect\n";
  matrix3D acopy(a);    // copy constructor
  acopy.show();

  std::cout << "test copy constructor using =  -- program will crash at end if this is incorrect\n";
  matrix3D a2copy = a;  // copy constructor
  a2copy.show();

  std::cout << "test assignment operator  -- program will crash at end if this is incorrect\n";
  matrix3D bcopy;
  bcopy = b;        // assignment operator
  bcopy.show();

  std::cout << "test negative unary operator\n";
  matrix3D aneg = -a;
  matrix3D bneg = -b;
  aneg.show();
  bneg.show();

  std::cout << "test determinant\n";
  printf("|a| = %.2f\n", a.determinant());
  printf("|b| = %.2f\n", b.determinant());

  // a.transpose().show();
  // b.transpose().show();

  std::cout << "test minors\n";
  a.minors().show();
  b.minors().show();

  std::cout << "test cofactor\n";
  a.cofactor().show();
  b.cofactor().show();

  std::cout << "test adjoint\n";
  a.adjoint().show();
  b.adjoint().show();

  std::cout << "test inverse\n";
  matrix3D ainv = a.inverse();
  ainv.show();

  matrix3D binv = b.inverse();
  binv.show();

  std::cout << "test a * ainv, b * binv\n";
  (a * ainv).show();
  (b * binv).show();

  // std::cout << "test a * ainv == id\n";
  matrix3D::identity(3).show();

  assert(a * ainv == matrix3D::identity(3));
  assert(a * ainv == ainv * a);
  assert(b * binv == matrix3D::identity(3));
  assert(b * binv == binv * b);

  std::cout << "test a.transpose().tranpose() == a\n";
  assert(a.transpose().transpose() == a);

  std::cout << "test a.determinant() == a.transpose().determinant()\n";
  assert(a.determinant() == a.transpose().determinant());

  std::cout << "test a + b == b + a\n";
  assert(a + b == b + a);

  std::cout << "test a - b == -(b - a)\n";
  assert(a - b == -(b - a));

  std::cout << "test 3.0 + a == a + 3.0\n";
  assert(3.0 + a == a + 3.0);

  std::cout << "test 3.0 * a == a * 3.0\n";
  assert(3.0 * a == a * 3.0);

  std::cout << "test a + 3.0 - 3.0 == a\n";
  matrix3D a_plus_3 = a + 3.0;
  a_plus_3.show();
  matrix3D a_minus_3 = a - 3.0;
  a_minus_3.show();

  matrix3D a_plus_3_minus_3 = a + 3.0 - 3.0;
  a_plus_3_minus_3.show();
  assert((a + 3.0) - 3.0 == a);

  std::cout << "test a * 3.0 / 3.0 == a\n";
  matrix3D a_times_3 = a * 3.0;
  a_times_3.show();
  matrix3D a_dividedby_3 = a / 3.0;
  a_dividedby_3.show();

  assert((3.0 * a) / 3.0 == a);

  std::cout << "-(-a) == a\n";
  assert(-(-a) == a);

  std::cout << "test matrix(1 2 3   4 5 6   7 8 9) has determinant 0\n";
  matrix3D zerod("zerod", 3, {1, 2, 3,   4, 5, 6,   7, 8, 9});
  assert(zerod.determinant() == 0);

  std::cout << "testing matrix vector multiplication\n";
  vector3D p("p", 2, {1, 2});
  matrix3D m("m", 2, {1, 2,   3, 4});
  p.show();
  m.show();
  std::cout << "asserting that p * m == m * p\n";
  assert(p * m == m * p);

  vector3D q("q", 3, {1, 2, 3});
  matrix3D n("n", 3, {1, 2, 3,   4, 5, 6,   7, 8, 9});
  q.show();
  n.show();
  std::cout << "asserting that q * n == n * q\n";
  assert(q * n == n * q);

  std::cout << "\n\n...test matrices assertions passed" << "\n\n";
  std::cout << "====================  FINISHED testing matrices  ============================================" << "\n\n";

  std::cout << "\n                >>>>>>>>>>  CONGRATULATIONS -- all assertions passed  <<<<<<<<<< \n\n\n";
}


private:
  void check_equal_dims(const matrix3d<T>& v) const; void check_bounds(int i) const;
  void swap(T& x, T& y);

private:
  std::string name_;
  int dims_;
  vector3d<T> cols_[4];
  T data_[16];

  static double epsilon_;
};

template<typename T> double matrix3d<T>::epsilon_ = 1e-10;



//================================================================================================\
// Implementation
//================================================================================================\

template <typename T> matrix3d<T>::matrix3d() : matrix3d("", 3) { }
template <typename T> matrix3d<T>::matrix3d(const std::string& name,
                                            int dims)
: name_(name), dims_(dims) {
  for (int i = 0; i < 4; ++i) {
    cols_[i].name("col" + std::to_string(i));
  }
  std::memset(data_, 0, 16 * sizeof(T));
}
template <typename T> matrix3d<T>::matrix3d(const std::string& name, int dims,
                                            const std::initializer_list<vector3d<T>>& li)
: matrix3d(name, dims) {
  int i = 0;
  for (vector3d<T> value : li) {
    if (i > dims_) { break; }
    cols_[i++] = value;
  }
}

template <typename T> matrix3d<T>::matrix3d(const std::string& name, int dims,
                                            const std::initializer_list<T>& li)
: matrix3d(name, dims) {
  // std::cout << "in constructor for name: " << name << "\n";
  int i = 0;
  matrix3d& m = *this;
  for (const T& value : li) { m(i / 3, i % 3) = value;  ++i; }
}
//=================================================================================================
template <typename T> matrix3d<T>& matrix3d<T>::operator=(T array[9]) {
  matrix3d& m = *this;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++i) {
      m(i, j) = array[i + j];
    }
  }
  return *this;
}
template <typename T> matrix3d<T>& matrix3d<T>::operator=(T k)
{ return operator=({ k, k, k,  k, k, k,  k, k, k }); }

//=================================================================================================
template <typename T> vector3d<T> matrix3d<T>::operator[](int i) const {
  check_bounds(i);
  return cols_[i];
}
template <typename T> vector3d<T>& matrix3d<T>::operator[](int i) {
  check_bounds(i);
  return cols_[i];
}
template <typename T> T matrix3d<T>::operator()(int row, int col) const { return cols_[col][row]; }
// Todo
template <typename T> T& matrix3d<T>::operator()(int row, int col) { 
  if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
      throw std::out_of_range("matrix3d::operator(): index out of range");
  }
  return data_[row * cols_ + col];
}
 // Todo
template <typename T> T* matrix3d<T>::opengl_memory(int row, int col) { 
  if (row >= rows_ || col >= cols_ || row < 0 || col < 0) {
    throw std::out_of_range("Index out of range");
  }
  return &data_[row * cols_ + col];
}
// implement code here
//=================================================================================================
template <typename T> void matrix3d<T>::name(const std::string& name) { name_ = name; }
template <typename T> const std::string& matrix3d<T>::name() const { return name_; }
//=================================== LINEAR ALGEBRA ================================
template <typename T> matrix3d<T>& matrix3d<T>::operator+=(T k) {
  matrix3d<T>& a = *this;
  name_ = std::to_string(k) + "+" + name_;
  for (int i = 0; i < 4; ++i) {
    a[i] += k;
  }
  return *this;
} 
template <typename T> matrix3d<T>& matrix3d<T>::operator-=(T k) { operator+=(-k); }
template <typename T> matrix3d<T>& matrix3d<T>::operator*=(T k) { 
    for (int i = 0; i < size_; i++) {
        data_[i] *= k;
    }
    return *this;
 }
template <typename T> matrix3d<T>& matrix3d<T>::operator/=(T k) { 
  if (abs(k) < epsilon_) { throw new std::invalid_argument("divide by zero\n"); }
  return operator*=(1.0 / k); 
}
//=================================================================================================
template <typename T> matrix3d<T>& matrix3d<T>::operator+=(const matrix3d<T>& b) { 
  check_equal_dims(b);
  matrix3d<T>& a = *this;
  a.name_ = a._name_ + " + " + b.name_;
  for (int i = 0; i < 4; ++i) {
    a[i] += b[i];
  }
  return *this;
}
template <typename T> matrix3d<T>& matrix3d<T>::operator-=(const matrix3d<T>& b) { return operator+=(b * -1); }
//=================================================================================================
template <typename T> matrix3d<T> matrix3d<T>::operator-() {
  const matrix3d<T>& a = *this;
  return matrix3d<T>("-" + name_, 3, { -a[0], -a[1], -a[2] });
}
template <typename T> matrix3d<T> matrix3d<T>::operator+(const matrix3d<T>& b) {
  const matrix3d<T>& a = *this;
  check_equal_dims(b);
  return matrix3d<T>(name_ + "+" + b.name_, dims_,
                    { a[0] + b[0], a[1] + b[1], a[2] + b[2] });
}
template <typename T> matrix3d<T> matrix3d<T>::operator-(const matrix3d<T>& b) { 
  if (dims_ != b.dims_) {
    throw std::invalid_argument("Matrix dimensions do not match");
  }

  matrix3d<T> result(name_, dims_);
  for (int i = 0; i < size_; i++) {
    result.data_[i] = data_[i] - b.data_[i];
  }
  return result;
}
//=================================================================================================
template <typename T> matrix3d<T> matrix3d<T>::operator*(const matrix3d<T>& b) {
  matrix3d res(name_ + "+" + b.name_, 3);
  const matrix3d<T>& a = *this;
  for (int i = 0; i < dims_; ++i) {
    for (int j = 0; j < dims_; ++j) {
      T result = T();
      for (int k = 0; k < dims_; ++k) {
        result += a(i, k) * b(k, j);         // SUM OVER K(a[i][k] * b[k][j]
      }
      res(i, j) = result;
    }
  }
  return res;
}
//=================================================================================================
template <typename T> matrix3d<T> matrix3d<T>::transpose() const {
  const matrix3d<T>& m = *this;
  matrix3d res(name_ + "Transpose", 3, {m(0,0), m(1,0), m(2,0),  m(0,1), m(1,1), m(2,1),  
                                m(0,2), m(1,2), m(2,2),  m(0,3), m(1,3), m(2,3)} );
  return res;
}
// Todo
template <typename T> T matrix3d<T>::determinant() const { 
  T a11 = data_[0], a12 = data_[1], a13 = data_[2];
  T a21 = data_[3], a22 = data_[4], a23 = data_[5];
  T a31 = data_[6], a32 = data_[7], a33 = data_[8];
  
  return a11*(a22*a33 - a23*a32) - a12*(a21*a33 - a23*a31) + a13*(a21*a32 - a22*a31);
 }
template <typename T> T matrix3d<T>::trace() const {
        const matrix3d<T>& m = *this;
        return m(0, 0) + m(1, 1) + m(2, 2);
}
//=================================================================================================
//||ef| |df| |de||
//||hi| |gi| |gh||
// | |
//||bc| |ac| |ab||
//||hi| |gi| |gh||
// | |
//||bc| |ac| |ab||
//||ef| |df| |de|| //---------------------------------------------------------------
template <typename T> matrix3d<T> matrix3d<T>::minors() const {
  const matrix3d<T>& m = *this;
  return matrix3d<T>(name_ + "_Minors", 3, {
    m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1),
    m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0),
    m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0),

    m(0, 1) * m(2, 2) - m(0, 2) * m(2, 1),
    m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0),
    m(0, 0) * m(2, 1) - m(0, 1) * m(2, 0),

    m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1),
    m(0, 0) * m(1, 2) - m(0, 2) * m(1, 0),
    m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0),
  });
}
template <typename T> matrix3d<T> matrix3d<T>::cofactor() const { 
  matrix3d<T> c;

  c(0, 0) =  det_2x2((*this)(1, 1), (*this)(1, 2), (*this)(2, 1), (*this)(2, 2));
  c(0, 1) = -det_2x2((*this)(1, 0), (*this)(1, 2), (*this)(2, 0), (*this)(2, 2));
  c(0, 2) =  det_2x2((*this)(1, 0), (*this)(1, 1), (*this)(2, 0), (*this)(2, 1));
  c(1, 0) = -det_2x2((*this)(0, 1), (*this)(0, 2), (*this)(2, 1), (*this)(2, 2));
  c(1, 1) =  det_2x2((*this)(0, 0), (*this)(0, 2), (*this)(2, 0), (*this)(2, 2));
  c(1, 2) = -det_2x2((*this)(0, 0), (*this)(0, 1), (*this)(2, 0), (*this)(2, 1));
  c(2, 0) =  det_2x2((*this)(0, 1), (*this)(0, 2), (*this)(1, 1), (*this)(1, 2));
  c(2, 1) = -det_2x2((*this)(0, 0), (*this)(0, 2), (*this)(1, 0), (*this)(1, 2));
  c(2, 2) =  det_2x2((*this)(0, 0), (*this)(0, 1), (*this)(1, 0), (*this)(1, 1));
  return c;
}

template <typename T>
T det_2x2(T a, T b, T c, T d) {
    return a * d - b * c;
}

template <typename T> matrix3d<T> matrix3d<T>::adjoint() const {
  
}
template <typename T> matrix3d<T> matrix3d<T>::inverse() const { /* TODO */ }
//=================================================================================================
template <typename T> matrix3d<T> matrix3d<T>::identity(int dims) { /* TODO */ }
template <typename T> matrix3d<T> matrix3d<T>::zero(int dims) { /* TODO */ }
template <typename T> bool matrix3d<T>::operator==(const matrix3d<T>& b) const {
  check_equal_dims(b);
  const matrix3d<T>& a = *this;
  T error = T();
  for (int i = 0; i < dims_; ++i) {
    for (int j = 0; j < dims_; ++j) {
      error += abs((double)(a(i, j) - b(i, j)));
    }
  }
  return error < epsilon_;
}
template <typename T> bool matrix3d<T>::operator!=(const matrix3d<T>& b) const {
  return !(*this == b);
}
//=================================================================================================
// Matrix of minors
//=================================================================================================


//=================================================================================================
template <typename T> void matrix3d<T>::check_equal_dims(const matrix3d<T>& v) const {
  if (dims_ != v.dims_) {
    throw new std::invalid_argument("matrix3d dims mismatch");
  }
}
template <typename T> void matrix3d<T>::check_bounds(int i) const {
  if (i > dims_) {
    throw new std::invalid_argument("out of bounds");
  }
}
template <typename T> void matrix3d<T>::swap(T& x, T& y) {
  T temp = x;
  x = y;
  y = temp;
}


#endif  /* __matrix_3d__ */
