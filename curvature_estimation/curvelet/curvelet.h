// This is curvelet.h
#ifndef curvelet_h
#define curvelet_h

/*****************************************************************************
// file: curvelet.h
// brief: A class to represent an edgel grouping with an associated curve model
//        (curvelet now called curve bundle)
// author: Xiaoyan Li
// date: 01/26/2015
******************************************************************************/

#include <list>
#include <deque>
#include "edgemap.h"
#include "CC_curve_model_3d.h"
#include "Array.h"
//#include "ref_count.h"


//: The curvelet class: stores the ordered list of edgels defining the curvelet
//  and also the curve model defined by the grouping
//  It also stores a list of higher order curvelets that it forms.
class curvelet
{
public:
    edgel* _ref_edgel;                        ///< ref edgel (the edgel to which it is anchored)
    
    edgel_chain _edgel_chain;          ///< the ordered list of edgels
    CC_curve_model_3d* _curve_model;                ///< associated curve model
    
    bool _forward;                                  ///< is this a forward cvlet or a reverse cvlet
    double _length;                                 ///< length of the curvelet
    double _quality;                                ///< the quality of this grouping (determined by various means)
    
    bool _used;                                     ///< to keep track of whether this curvelet was used in linking
    
    //: default constructor
    curvelet(): _ref_edgel(0), _edgel_chain(0), _curve_model(0), _forward(true), _length(0.0), _quality(0.0){}
    
    //: constructor 1
    curvelet(edgel* e, CC_curve_model_3d* cm, edgel_chain &echain, bool dir=true) :
    _ref_edgel(e), _curve_model(cm), _forward(dir), _length(0.0), _quality(0.0)
    {
        _edgel_chain.insert(_edgel_chain.end(), echain.begin(), echain.end());
    }
    
    //: copy constructor
    curvelet(const curvelet& other);
    
    //: destructor
    ~curvelet();
    
    //: return the order of this grouping
    unsigned order() const { return _edgel_chain.size(); }
    
    //: compute properties of this curvelet once formed
    void compute_properties(double R, double token_len);
    
    //: update length according to the curve model angle
    // if _length<0, the arc is in an opposite direction
    void update_length();
    
    //: record the curvelet map data in array format
    void set_output(arrayi &id_chain, arrayd &curvelet_info, unsigned posy, unsigned posx);
    
    //: print info to file
    void print(std::ostream&);
    
    //: return intersection of two curvelets
    curvelet* intersect(curvelet* c2);

};

typedef std::list<curvelet*> cvlet_list;

#endif // curvelet_h
