// This is CC_curve_model_3d.h
#ifndef CC_curve_model_3d_h
#define CC_curve_model_3d_h

/*****************************************************************************
// file: CC_curve_model_3d.h
// brief: 3d curve mdoels used for linking using gemoetric consistency
// author: Xiaoyan Li
// date: 01/26/2015
******************************************************************************/

#include <vector>
#include <deque>
#include "edgemap.h"
#include "Array.h"

//Note: The grid size of the perturbations in position and orientation should be a parameter
// in the algorithm but for now I'm gonna keep them as constants
#undef _dx_
#define _dx_ 0.1 //pixels

#undef _dt_
#define _dt_ 0.08 //radians

//: The circular arc curve model class (3d curve bundle)
class CC_curve_model_3d
{
public:
    arrayd _Kmin, _Kmax; ///< The min and max curvature surfaces for the bundle (as a function of dx, dt)
    //vcl_vector<double> Dx, Dt;       ///< The values of the Dx and Dt parameters at the grid points
    
    bool _is_poly_arc_formed;  ///< Use poly arc chain to set best fit
    
    // Ref PT are enforced to save K range at last !!!!
    point_2d _ref_pt;   ///< Extrinsic anchoring position
    double _ref_theta;  ///< Extrinsic anchoring orientation
    
    point_2d _pt;       ///< best position estimate
    double _theta;      ///< best orientation estimate
    double _k;          ///< best estimate of curvature
    point_2d _center;   ///< center estimate
    double _angle;      ///< angle estimate
    
    
    //: Constructor: From a pair of edgels
    CC_curve_model_3d(edgel* ref_e, edgel* ne, double dx, double dt, double max_k, bool verbose):
    _is_poly_arc_formed(false),
    _ref_pt(ref_e->_pt), _ref_theta(ref_e->_orientation), _pt(0.0,0.0), _theta(0), _k(0.0),
    _center(0.0,0.0), _angle(0.0)
    {
        //construct the curve bundle from the pair of edgels and the given uncertainty
        compute_curve_bundle(ref_e, ne, dx, dt, max_k, verbose);
    }
    
    //: constructor 3: From the intersection of two curve bundles
    CC_curve_model_3d(CC_curve_model_3d* cm1, CC_curve_model_3d* cm2)
    :
    _is_poly_arc_formed(false),
    _ref_pt(cm1->_ref_pt),
    _ref_theta(cm1->_ref_theta),
    _pt(0.0,0.0),
    _theta(0),
    _k(0.0),
    _center(0.0,0.0),
    _angle(0.0)
    {
        //take the intersection of the two bundles by intersecting the four bounding surfaces
        _Kmin.array_max(cm1->_Kmin, cm2->_Kmin);
        _Kmax.array_min(cm1->_Kmax, cm2->_Kmax);
    }
    
    //: constructor 5: Create default bundle only
    CC_curve_model_3d(edgel* ref_e, double dpos, double dtheta, double max_k);
    
    //: copy constructor
    CC_curve_model_3d(const CC_curve_model_3d& other);
    
    virtual ~CC_curve_model_3d(){};
    //: construct and return a curve model of the same type by intersecting with another curve bundle
    CC_curve_model_3d* intersect(CC_curve_model_3d* cm);
    
    //: get the value of the position perturbation at a grid point
    double Dx(unsigned i){ return _dx_*((double)i-(double(_Kmax.h())-1.0)/2.0); }
    
    //: get the value of the orientation perturbation at a grid point
    double Dt(unsigned j){ return _dt_*((double)j-(double(_Kmax.w())-1.0)/2.0); }
    
    //: is this bundle valid?
    virtual bool bundle_is_valid()
    {
        //if the Kmax surface is completely less than the Kmin surface, there is no bundle
        for (unsigned i=0; i<_Kmax.h(); ++i)
            for (unsigned j=0; j<_Kmax.w(); ++j)
                if (_Kmax.val(i,j) > _Kmin.val(i,j)) return true;
        
        return false;
    }
    
    //: compute the CC curve bundle for an edgel pair given expected errors
    void compute_curve_bundle(edgel* ref_e, edgel* ne, double dpos, double dtheta, double max_k, bool verbose);
    
    //: Compute the best fit curve from the curve bundle
    virtual point_2d compute_best_fit();
    
    //: Set the best fit curve
    virtual void set_best_fit(point_2d dx_dt, double k);
    
    //: Compute the center of the circular arc
//    virtual void compute_center_and_angle(const point_2d &ref_tail);
    
    //: update the best fit curve according to points on the curve
    virtual void update_best_fit(point_2d p1, point_2d p2);
    
    //: function to check if the curve fit is reasonable
    virtual bool curve_fit_is_reasonable(std::deque<edgel*> &edgel_chain, edgel* ref_e, double dpos);
    
    //: are these two curve bundles C^2?
    virtual bool is_C2_with(CC_curve_model_3d* cm);
    
    //: are these two curve bundles C^1?
    virtual bool is_C1_with(CC_curve_model_3d* cm);
    
    //: resize the curve model when the image is resized
    void resize(double scale);
    
    //: report accuracy of measurement
    // virtual void report_accuracy(double *estimates, double *min_estimates, double *max_estimates);
    
    // utility functions
    
    //: get the size of the bundle for computing the saliency heuristic
    //double get_CC_3d_bundle_area();
    
    //: record the curvelet map data in array format
    void set_output(arrayd &curvelet_info, unsigned pos);
    
    //: print bundle info
    virtual void print_bundle(std::ostream&  os);
    
    //: print central info to file
    virtual void print(std::ostream&);
    
    //: read central info from file
    //virtual void read(vcl_istream&);
};

#endif // curve_model_3d_h
