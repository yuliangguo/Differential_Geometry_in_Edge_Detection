// This is knot.h
#ifndef knot_h_
#define knot_h_

/*****************************************************************************
 // file: knot.h
 // brief: A data structure for representing knot points
 // Reference: lemsvxl code by Amir Tamrakar
 // author: Xiaoyan Li
 // date: 02/03/2016
******************************************************************************/

#include <vector>
#include <iostream>
#include <deque>
#include <cmath>

class knot
{
public:
    //dbdet_curvelet L_cvlet;
    //dbdet_curvelet R_cvlet;
    
    unsigned nB; ///< # of edgels before anchor
    unsigned nA; ///< # of edgels after anchor
    
    //: default constructor
    knot(unsigned nb=0, unsigned na=0): nB(nb), nA(na) {}
    
    //: destructor
    virtual ~knot(){}
    
};

#endif // dbdet_sel_knot_h
