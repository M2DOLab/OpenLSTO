/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//FileName:mpVector.h
//Author:Michael Y. Polyakov
//email:myp@andrew.cmu.eduor  mikepolyakov@hotmail.com
//Website:www.angelfire.com/linux/myp
//Date:7/16/2002
//
//Provides basic vector handling.
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MPVECTOR_H
#define MPVECTOR_H

#include <math.h>
#include <stdio.h>
//#include <wx/string.h>

class mp4Vector;

class mpVector
{
 public:
  float x, y, z;

  mpVector();
  mpVector(float xx, float yy, float zz);
  mpVector(const mpVector& other);

  mpVector& Normalize();
  float Magnitude();
  mpVector Cross(const mpVector& other);
  mpVector operator - (mpVector v);
  mpVector operator + (mpVector v);
  float operator * (mpVector v);
  mpVector operator - (float c);
  mpVector operator + (float c);
  mpVector operator / (float c);
  mpVector operator * (float c);
  mpVector& operator = (const mpVector& other);
  operator mp4Vector() const;

  operator char*()  const;

};

class mp4Vector
{
 public:
  float x, y, z, val;

  mp4Vector();
  mp4Vector(float aa, float bb, float cc, float dd);
  mp4Vector(const mp4Vector& other);
  mp4Vector(const mpVector& v, const float value);

  void operator = (const mp4Vector& v);

  operator mpVector() const;
};
#endif
