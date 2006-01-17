// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file polynom.h
    Datatypes and functions for polynomials.
    Evaluation possible both as Taylor and Chebychev series. Note that the length of the
    double list is equal to the order of the polynomial plus 1, so that Polynom->n does not give
    the order of the polynomial, but one more.
*/
#ifndef POLYNOM_H
#define POLYNOM_H
#include "utils.h"

/** basically, a polynomial is just a list of coefficients */
typedef DoubleList Polynom;

/** evaluate the polynomial interpreted as a Taylor series via the Horner scheme */
MDINLINE double evaluateAsTaylorSeriesAt(Polynom *series, double x)
{
  int cnt   = series->n - 1;
  double *c = series->e;
  double r  = c[cnt];
  while (--cnt >= 0)
    r = r*x + c[cnt];
  return r;
}

/** evaluate the polynomial interpreted as a Chebychev series. Requires a series with at least
    three coefficients, i.e. no linear approximations! */
MDINLINE double evaluateAsChebychevSeriesAt(Polynom *series, double x)
{
  int j;
  double *c = series->e;
  double x2 = 2.0 * x;
  double dd = c[series->n - 1];
  double d  = x2*dd + c[series->n - 2];
  for(j = series->n - 3; j >= 1; j--) {
    double tmp = d;
    d = x2*d - dd + c[j];
    dd = tmp;
  }
  return x*d - dd + 0.5 * c[0];
}

#endif
