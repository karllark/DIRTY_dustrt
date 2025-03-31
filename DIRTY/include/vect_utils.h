#ifndef _VECT_UTILS_H_INCLUDED
#define _VECT_UTILS_H_INCLUDED

/*******************************************************************************
 * vect_utils:  A collection of classes to implement multi-dimensional arrays. *
 *                                                                             *
 * Includes:                                                                   *
 *   class Matrix                                                              *
 *   class Tensor                                                              *
 *   class FourVector                                                          *
 *   class FiveVector                                                          *
 *                                                                             *
 * All are template classes and so can be assigned any C++ type.  I have only  *
 * tested them for primitive types though (ints, floats, doubles, etc), so be  *
 * careful.                                                                    *
 *                                                                             *
 * These are all sub-classed from the STL vector.  The elements are stored as  *
 * 1d vectors and hence occupy contiguous blocks of memory.  That makes for a  *
 * memory advantage in that no extra arrays of pointers need to allocated; it  *
 * also has speed advantages in that contiguous blocks can be easily copied    *
 * using fast built in functions like memcpy.  For example assume you want to  *
 * copy one matrix to another.  Rather than looping over elements and testing  *
 * conditionals at each iteration, you could simply do something like the      *
 * following:                                                                  *
 *                                                                             *
 *    #include <memory>                                                        *
 *    #include "vect_utils.h"                                                  *
 *                                                                             *
 *    Matrix <double> old_matrix(n1,n2);                                       *
 *    Matrix <double> new_matrix(n1,n2);                                       *
 *                                                                             *
 *    memcpy(old_matrix.begin(),new_matrix.begin(),sizeof(double)*n1*n2);      *
 *                                                                             *
 * Which is very fast.  Two things to note;  no checking is done here to see   *
 * that you've not done something foolish like tried to copy vectors with      *
 * a different number of total elements or different dimensions into each      *
 * other. Ideally, I should write a member function that copies each class and *
 * checks such things.  But that's for a later day.  Next - since the classes  *
 * inherit from vector, all vector member functions are available to them; ie. *
 * .begin(), .size(), .insert(), .reserve() etc.                               *
 *                                                                             *
 * Other important things. I'm assuming the fits standard for data storage.    *
 * I've overloaded operator() to access/store elements within the structure    *
 * where elements are located/placed assuming first dimension changes the most *
 * quickly.  For example, given naxis1, naxis2, ....., naxisn are the number   *
 * elements in each axis of your data structure.  They will be accessed/stored *
 * in the following manner:                                                    *
 *                                                                             *
 * A(0,0,0,...0) A(1,0,0,...0) ..... A(naxis1-1,0,0,....0) --->                *
 * A(0,1,0,...0) A(1,1,0,...0) ..... A(naxsi1-1,1,0,....0) ---> ....           *
 * A(0,naxis2-1,0,...0), A(1,naxis2-1,0,...0),....A(naxis1-1,naxis2-1,0,...0)  *
 * etc....                                                                     *
 * Note that this doesn't look like what you normally think of for mathmatic   *
 * arrays (ie. the first dimension looks more like a column than a row) but is *
 * more like x,y,z... coordinates - image like even.  In any case, accessing   *
 * and storing data should be transparent.  Once you've allocated your array,  *
 * simply use () - if its a Matrix, the something like matrix(i,j) works.  For *
 * a Tensor, Tensor(i,j,k) and for a FourVector, FourVector(i,j,k,l).          *
 * I have put in some rudimentary exception handling to make sure that you     *
 * don't try to access elements beyond what you allocated.  I'm only using the *
 * STL exception out_of_range. For example, doing the following will crash:    *
 *                                                                             *
 * Matrix M(5,5);                                                              *
 * M(5,2) = 10;  // Fails because you're trying to access a sixth i element    *
 *               // of an array you've only allocated five i elements for.     *
 * M(4,2) = 10;  // OK.                                                        *
 *                                                                             *
 * Note that the first case will core dump without any real neat information.  *
 * In order to get more diagnostic information about what the problem is, you  *
 * could do something like the following:                                      *
 *                                                                             *
 * int n1=28,n2=2,n3=500;                                                      *
 * Tensor T(n1,n2,n3);                                                         *
 * for (int i=0;i<n1;i++)                                                      *
 *   for (int j=0;j<n2;j++)                                                    *
 *     for (int k=0;k<=n3;k++)                                                 *
 *       try {                                                                 *
 *         cout << T(i,j,k) << endl;                                           *
 *       }                                                                     *
 *       catch (const out_of_range &excp) {                                    *
 *         cerr << excp.what() << endl;                                        *
 *       }                                                                     *
 *                                                                             *
 * Since we are going beyond the bounds of the initial allocation in the 3rd   *
 * dimension, this code will throw an out_of_bounds exception and you'll get   *
 * the message: "out_of_range Tensor::operator(), element 3" on the crash.     *
 *                                                                             *
 * History:                                                                    *
 *          Added five-vector class.                   07/26/03 KAM            *
 *          Added after instatiation sizing functions. 02/06/02 KAM            *
 *          Written.                                      12/01 KAM            *
 ******************************************************************************/

#include <stdexcept>
#include <vector>

#include <iostream>

using namespace std;

/*******************************************************************************
 * Matrix class.                                                               *
 ******************************************************************************/
template <typename T> class Matrix;

template <typename T> class Matrix : public std::vector<T>
{

public:
  // Constructors/destructors.
  Matrix () : std::vector<T> () {}
  Matrix (int n1, int n2, const T &ival) : std::vector<T> (n1 * n2, ival), _n1 (n1), _n2 (n2) {}
  explicit Matrix (int n1, int n2) : std::vector<T> (n1 * n2), _n1 (n1), _n2 (n2) {}

  ~Matrix () {}

  // Return reference to correct element in vector; All vector type assignments
  // should work.
  T &operator() (int n1_id, int n2_id);
  void MSize (int n1_id, int n2_id);

private:
  // Keep track of how many elements we're allowed to have in each dimension..
  int _n1, _n2;
};

// Overload of operator()
template <typename T>
inline T &
Matrix<T>::operator() (int n1_id, int n2_id)
{
  if (n1_id < 0 || n1_id >= _n1)
    {
      cout << "out_of_range Matrix::operator(), element 1" << endl;
      std::string ExceptionObject = "out_of_range Matrix::operator(), element 1";
      throw std::out_of_range (ExceptionObject);
    }
  if (n2_id < 0 || n2_id >= _n2)
    {
      cout << "out_of_range Matrix::operator(), element 2" << endl;
      std::string ExceptionObject = "out_of_range Matrix::operator(), element 2";
      throw std::out_of_range (ExceptionObject);
    }
  return *(this->begin () + n2_id * _n1 + n1_id);
}

// Size the Matrix after instatiation.
template <typename T>
void
Matrix<T>::MSize (int n1_id, int n2_id)
{
  _n1 = n1_id;
  _n2 = n2_id;
  Matrix::clear ();
  Matrix::resize (n1_id * n2_id);
  // Matrix::clear();
}

/******************************************************************************/

/*******************************************************************************
 * Tensor class.                                                               *
 ******************************************************************************/
template <typename T> class Tensor;

template <typename T> class Tensor : public std::vector<T>
{

public:
  // Constructors/destructors.
  Tensor () : std::vector<T> () {}
  Tensor (int n1, int n2, int n3, const T &ival)
      : std::vector<T> (n1 * n2 * n3, ival), _n1 (n1), _n2 (n2), _n3 (n3)
  {
  }
  explicit Tensor (int n1, int n2, int n3)
      : std::vector<T> (n1 * n2 * n3), _n1 (n1), _n2 (n2), _n3 (n3)
  {
  }

  ~Tensor () {}

  // Return reference to correct element in vector; All vector type assignments
  // should work.
  T &operator() (int n1_id, int n2_id, int n3_id);
  void TSize (int n1_id, int n2_id, int n3_id);

private:
  // Keep track of how many elements we're allowed to have in each dimension.
  int _n1, _n2, _n3;
};

// Overload of operator()
template <typename T>
inline T &
Tensor<T>::operator() (int n1_id, int n2_id, int n3_id)
{
  if (n1_id < 0 || n1_id >= _n1)
    {
      cout << "out_of_range Tensor::operator(), element 1" << endl;
      std::string ExceptionObject = "out_of_range Tensor::operator(), element 1";
      throw std::out_of_range (ExceptionObject);
    }
  if (n2_id < 0 || n2_id >= _n2)
    {
      cout << "out_of_range Tensor::operator(), element 2" << endl;
      std::string ExceptionObject = "out_of_range Tensor::operator(), element 2";
      throw std::out_of_range (ExceptionObject);
    }
  if (n3_id < 0 || n3_id >= _n3)
    {
      cout << "out_of_range Tensor::operator(), element 3" << endl;
      std::string ExceptionObject = "out_of_range Tensor::operator(), element 3";
      throw std::out_of_range (ExceptionObject);
    }
  return *(this->begin () + n3_id * _n1 * _n2 + n2_id * _n1 + n1_id);
}

// Size the Tensor after instatiation.
template <typename T>
void
Tensor<T>::TSize (int n1_id, int n2_id, int n3_id)
{
  _n1 = n1_id;
  _n2 = n2_id;
  _n3 = n3_id;
  // Tensor::clear();
  Tensor::resize (n1_id * n2_id * n3_id);
  // Tensor::clear();
}
/******************************************************************************/

/*******************************************************************************
 * FourVector class.                                                           *
 ******************************************************************************/
template <typename T> class FourVector;

template <typename T> class FourVector : public std::vector<T>
{

public:
  // Constructors/destructors.
  FourVector () : std::vector<T> () {}
  FourVector (int n1, int n2, int n3, int n4, const T &ival)
      : std::vector<T> (n1 * n2 * n3 * n4, ival), _n1 (n1), _n2 (n2), _n3 (n3), _n4 (n4)
  {
  }
  explicit FourVector (int n1, int n2, int n3, int n4)
      : std::vector<T> (n1 * n2 * n3 * n4), _n1 (n1), _n2 (n2), _n3 (n3), _n4 (n4)
  {
  }

  ~FourVector () {}

  // Return reference to correct element in vector; All vector type assignments
  // should work.
  T &operator() (int n1_id, int n2_id, int n3_id, int n4_id);
  void FVSize (int n1_id, int n2_id, int n3_id, int n4_id);

private:
  // Keep track of how many elements we're allowed to have in each dimension.
  int _n1, _n2, _n3, _n4;
};

// Overload of operator()
template <typename T>
inline T &
FourVector<T>::operator() (int n1_id, int n2_id, int n3_id, int n4_id)
{
  // thows not working... do it the "clutz way".
  if (n1_id < 0 || n1_id >= _n1)
    {
      cout << "out_of_range FourVector::operator(), element 1" << endl;
      std::string ExceptionObject = "out_of_range FourVector::operator(), element 1";
      throw std::out_of_range (ExceptionObject);
    }
  if (n2_id < 0 || n2_id >= _n2)
    {
      cout << "out_of_range FourVector::operator(), element 2" << endl;
      std::string ExceptionObject = "out_of_range FourVector::operator(), element 2";
      throw std::out_of_range (ExceptionObject);
    }
  if (n3_id < 0 || n3_id >= _n3)
    {
      cout << "out_of_range FourVector::operator(), element 3" << endl;
      std::string ExceptionObject = "out_of_range FourVector::operator(), element 3";
      throw std::out_of_range (ExceptionObject);
    }
  if (n4_id < 0 || n4_id >= _n4)
    {
      cout << "out_of_range FourVector::operator(), element 4" << endl;
      std::string ExceptionObject = "out_of_range FourVector::operator(), element 4";
      throw std::out_of_range (ExceptionObject);
    }
  return *(this->begin () + n4_id * _n1 * _n2 * _n3 + n3_id * _n1 * _n2 + n2_id * _n1 + n1_id);
}

// Size the FourVector after instatiation.
template <typename T>
void
FourVector<T>::FVSize (int n1_id, int n2_id, int n3_id, int n4_id)
{
  _n1 = n1_id;
  _n2 = n2_id;
  _n3 = n3_id;
  _n4 = n4_id;
  // zero size, resize, zero size!
  FourVector::clear ();
  FourVector::resize (n1_id * n2_id * n3_id * n4_id);
  // FourVector::clear();
}

/******************************************************************************/

/*******************************************************************************
 * Ah, FiveVector class?                                                       *
 ******************************************************************************/
template <typename T> class FiveVector;

template <typename T> class FiveVector : public std::vector<T>
{

public:
  // Constructors/destructors.
  FiveVector () : std::vector<T> () {}
  FiveVector (int n1, int n2, int n3, int n4, int n5, const T &ival)
      : std::vector<T> (n1 * n2 * n3 * n4 * n5, ival), _n1 (n1), _n2 (n2), _n3 (n3), _n4 (n4),
        _n5 (n5)
  {
  }
  explicit FiveVector (int n1, int n2, int n3, int n4, int n5)
      : std::vector<T> (n1 * n2 * n3 * n4 * n5), _n1 (n1), _n2 (n2), _n3 (n3), _n4 (n4), _n5 (n5)
  {
  }

  ~FiveVector () {}

  // Return reference to correct element in vector; All vector type assignments
  // should work.
  T &operator() (int n1_id, int n2_id, int n3_id, int n4_id, int n5_id);
  void FiVSize (int n1_id, int n2_id, int n3_id, int n4_id, int n5_id);

private:
  // Keep track of how many elements we're allowed to have in each dimension.
  int _n1, _n2, _n3, _n4, _n5;
};

// Overload of operator()
template <typename T>
inline T &
FiveVector<T>::operator() (int n1_id, int n2_id, int n3_id, int n4_id, int n5_id)
{
  // thows not working... do it the "clutz way".
  if (n1_id < 0 || n1_id >= _n1)
    {
      cout << "out_of_range FiveVector::operator(), element 1" << endl;
      std::string ExceptionObject = "out_of_range FiveVector::operator(), element 1";
      throw std::out_of_range (ExceptionObject);
    }
  if (n2_id < 0 || n2_id >= _n2)
    {
      cout << "out_of_range FiveVector::operator(), element 2" << endl;
      std::string ExceptionObject = "out_of_range FiveVector::operator(), element 2";
      throw std::out_of_range (ExceptionObject);
    }
  if (n3_id < 0 || n3_id >= _n3)
    {
      cout << "out_of_range FiveVector::operator(), element 3" << endl;
      std::string ExceptionObject = "out_of_range FiveVector::operator(), element 3";
      throw std::out_of_range (ExceptionObject);
    }
  if (n4_id < 0 || n4_id >= _n4)
    {
      cout << "out_of_range FiveVector::operator(), element 4" << endl;
      std::string ExceptionObject = "out_of_range FiveVector::operator(), element 4";
      throw std::out_of_range (ExceptionObject);
    }
  if (n5_id < 0 || n5_id >= _n5)
    {
      cout << "out_of_range FiveVector::operator(), element 5" << endl;
      std::string ExceptionObject = "out_of_range FiveVector::operator(), element 5";
      throw std::out_of_range (ExceptionObject);
    }
  return *(this->begin () + n5_id * _n1 * _n2 * _n3 * _n4 + n4_id * _n1 * _n2 * _n3
           + n3_id * _n1 * _n2 + n2_id * _n1 + n1_id);
}

// Size the FiveVector after instatiation, make zero size.
template <typename T>
void
FiveVector<T>::FiVSize (int n1_id, int n2_id, int n3_id, int n4_id, int n5_id)
{
  _n1 = n1_id;
  _n2 = n2_id;
  _n3 = n3_id;
  _n4 = n4_id;
  _n5 = n5_id;
  FiveVector::clear ();
  FiveVector::resize (n1_id * n2_id * n3_id * n4_id * n5_id);
  // FiveVector::clear();
}

/******************************************************************************/

#endif // already included conditional closed.
