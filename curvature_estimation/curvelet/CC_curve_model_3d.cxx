 #include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm> 

#include "CC_curve_model_3d.h"

extern double To2Pi (double angle);
extern double sq_distance(const point_2d& A, const point_2d& B);
extern double sq_norm(const point_2d& A);
extern const double k_th;

//*****************************************************************************//
// Circular Arc curve model (3d bundle)
//*****************************************************************************//

//: construct and return a curve model of the same type by intersecting with another curve bundle
CC_curve_model_3d* CC_curve_model_3d::intersect(CC_curve_model_3d* cm)
{
    return new CC_curve_model_3d(this, cm);
}

//: copy constructor
CC_curve_model_3d::CC_curve_model_3d(const CC_curve_model_3d& other):_Kmin(other._Kmin),_Kmax(other._Kmax)
{
    //Dx = other.Dx;
    //Dt = other.Dt;
    
    _ref_pt = other._ref_pt;
    _ref_theta = other._ref_theta;
    
    _pt = other._pt;
    _theta = other._theta;
    _k = other._k;
    _center = other._center;
    _angle = other._angle;
    _is_poly_arc_formed = other._is_poly_arc_formed;
}

//: constructor 5: Create default bundle only
CC_curve_model_3d::CC_curve_model_3d(edgel* ref_e, double dpos, double dtheta, double max_k):
_ref_pt(ref_e->_pt), _ref_theta(ref_e->_orientation), _pt(0.0,0.0), _theta(0), _k(0.0), _center(0.0,0.0),
_angle(0.0)
{
    //determine the size of the uncertainty at the reference edgel
    double ref_dx = dpos, ref_dt=dtheta;
    
    //compute the curvature bounds due to the token length and dpos
    //This is to be used for determining the min and max of the curvature
    //values for the default bundle
    // double dk = 2*ref_dx/(token_len*token_len/4.0 + ref_dx*ref_dx);
    
    //1) construct the bundle at the ref edgel and the other edgel
    
    // create the default bundle (omega: the restricted curve bundle space)
    // (currently assuming independence of the parameters)
    //
    // Also determine the size of the grid necessary to represent the curve bundle
    // given its uncertainty parameters
    
    _Kmin.init(int(2*floor(ref_dx/_dx_+0.5)+1), int(2*floor(ref_dt/_dt_+0.5)+1));
    _Kmin.fill(-max_k);
    
    _Kmax.init(int(2*floor(ref_dx/_dx_+0.5)+1), int(2*floor(ref_dt/_dt_+0.5)+1));
    _Kmax.fill(max_k);
    
    _is_poly_arc_formed = false;
}

//#ifdef fffUSE_APPROX_SINE
//
////always wrap input angle to -M_PI..M_PI
//inline double mysin_lowprecision(double x) {
//    if (x < -3.14159265)
//        x += 6.28318531;
//    else
//        if (x >  3.14159265)
//            x -= 6.28318531;
//    
//    //compute sine
//    if (x < 0)
//        return 1.27323954 * x + .405284735 * x * x;
//    else
//        return 1.27323954 * x - 0.405284735 * x * x;
//}
//
//inline double mycos_lowprecision(double x) {
//    if (x < -3.14159265)
//        x += 6.28318531;
//    else
//        if (x >  3.14159265)
//            x -= 6.28318531;
//    //compute cosine: sin(x + M_PI/2) = cos(x)
//    x += 1.57079632;
//    if (x >  3.14159265)
//        x -= 6.28318531;
//    
//    if (x < 0)
//        return 1.27323954 * x + 0.405284735 * x * x;
//    else
//        return 1.27323954 * x - 0.405284735 * x * x;
//}
//
//inline double mysin_precise(double x) {
//    //always wrap input angle to -M_PI..M_PI
//    if (x < -3.14159265)
//        x += 6.28318531;
//    else
//        if (x >  3.14159265)
//            x -= 6.28318531;
//    
//    //compute sine
//    if (x < 0)
//    {
//        double const sin = 1.27323954 * x + .405284735 * x * x;
//        
//        if (sin < 0)
//            return .225 * (sin *-sin - sin) + sin;
//        else
//            return .225 * (sin * sin - sin) + sin;
//    }
//    else
//    {
//        double const sin = 1.27323954 * x - 0.405284735 * x * x;
//        
//        if (sin < 0)
//            return .225 * (sin *-sin - sin) + sin;
//        else
//            return .225 * (sin * sin - sin) + sin;
//    }
//}
//
//inline double mycos_precise(double x) {
//    if (x < -3.14159265)
//        x += 6.28318531;
//    else
//        if (x >  3.14159265)
//            x -= 6.28318531;
//    //compute cosine: sin(x + M_PI/2) = cos(x)
//    x += 1.57079632;
//    if (x >  3.14159265)
//        x -= 6.28318531;
//    
//    if (x < 0)
//    {
//        double const cos = 1.27323954 * x + 0.405284735 * x * x;
//        
//        if (cos < 0)
//            return .225 * (cos *-cos - cos) + cos;
//        else
//            return .225 * (cos * cos - cos) + cos;
//    }
//    else
//    {
//        double const cos = 1.27323954 * x - 0.405284735 * x * x;
//        
//        if (cos < 0)
//            return .225 * (cos *-cos - cos) + cos;
//        else
//            return .225 * (cos * cos - cos) + cos;
//    }
//}
//
//// You can also try mysin_lowprecision for really fast speedups, but this will
//// generate largely different results.
//#define mysin(x) (mysin_precise((x)))
//#define mycos(x) (mycos_precise((x)))
//
//#else

#define mysin(x) (sin((x)))
#define mycos(x) (cos((x)))

//#endif

//: compute the CC curve bundle for an edgel pair at the ref edgel
// (transport default bundle at ne to ref_e)
// \Remarks This is an optimzied version of the implementation.  See the commented block
// after this function for the original implementation, if you think this
// optimized one has become unreadable.
void CC_curve_model_3d::
compute_curve_bundle(edgel* ref_e, edgel* ne, double dx, double dt, double max_k, bool verbose)
{
    //determine the size of the uncertainty at the reference edgel
    double ref_dx = dx, ref_dt=dt;
    
    //1) construct the bundle at the ref edgel and the other edgel
    
    // create the default bundle
    // (currently assuming independence of the parameters)
    _Kmin.init(2*floor(ref_dx/_dx_+0.5)+1, 2*floor(ref_dt/_dt_+0.5)+1);
    _Kmin.fill(-max_k);
    
    _Kmax.init(2*floor(ref_dx/_dx_+0.5)+1, 2*floor(ref_dt/_dt_+0.5)+1);
    _Kmax.fill(max_k);
    
    //now compute the four constraint surfaces due to the neighboring edgel
    const double DX = ref_e->_pt.first - ne->_pt.first;
    const double DY = ref_e->_pt.second - ne->_pt.second;
    const double DX2 = DX*DX;
    const double DY2 = DY*DY;
    
    //if neighboring point is close to the ref edgel, no static constraint can be computed of leave it as the default bundle
    if ( DX2 + DY2 < 4*dx*dx )
        return;
    
    const double T0 = ne->_orientation;
    const double T2 = ref_e->_orientation;
    
    //The uncertainty of the neighboring edgel (default uncertainty = adaptive)
    const double dx_sq = dx*dx;
    const double T0_p_dt = T0 + dt;
    const double T0_m_dt = T0 - dt;
    
    if ( (DX*mycos(T2)+DY*mysin(T2)) < 0 ) {
        for (unsigned i=0; i<_Kmax.h(); ++i){
            for (unsigned j=0; j<_Kmax.w(); ++j){
                const double dx2 = Dx(i);
                const double dt2 = Dt(j);
                
                const double T2t = T2+dt2;
                const double tt = mysin(T2t)*DX - mycos(T2t)*DY;
                
                // changed k_dx_min and k_dx_max according to the thesis
                const double denom = (-2.0*dx2*tt + DX2 + DY2 - dx_sq + dx2*dx2);
                
                const double ttmdx2 = tt - dx2;
                const double k_dx_min = 2.0*(ttmdx2-dx)/denom;
                //compute the intersection of all the surfaces to compute the final bundle
                if (_Kmin.val(i,j) < k_dx_min)
                    _Kmin.set_val(i, j, k_dx_min);
                
                const double k_dx_max = 2.0*(ttmdx2+dx)/denom;
                //compute the intersection of all the surfaces to compute the final bundle
                if (_Kmax.val(i,j) > k_dx_max)
                    _Kmax.set_val(i, j, k_dx_max);
                
                const double sin_t2_p_dt2_m_t0_m_dt = mysin(T2t-T0_m_dt);
                const double sin_t2_p_dt2_m_t0_p_dt = mysin(T2t-T0_p_dt);
                const double sin_T0pdt_times_DY = DY*mysin(T0_p_dt);
                const double cos_T0pdt_times_DX = DX*mycos(T0_p_dt);
                const double sin_T0mdt_times_DY = DY*mysin(T0_m_dt);
                const double cos_T0mdt_times_DX = DX*mycos(T0_m_dt);
                
                const double k_dt_min =sin_t2_p_dt2_m_t0_m_dt/(-sin_t2_p_dt2_m_t0_m_dt*dx2 + cos_T0mdt_times_DX + sin_T0mdt_times_DY);
                //compute the intersection of all the surfaces to compute the final bundle
                if (_Kmin.val(i,j) < k_dt_min)
                    _Kmin.set_val(i,j,k_dt_min);
                
                const double k_dt_max = sin_t2_p_dt2_m_t0_p_dt/(-sin_t2_p_dt2_m_t0_p_dt*dx2 + cos_T0pdt_times_DX + sin_T0pdt_times_DY);
                //compute the intersection of all the surfaces to compute the final bundle
                if (_Kmax.val(i,j) > k_dt_max)
                    _Kmax.set_val(i,j,k_dt_max);
                
                if (verbose)
                    std::cout << std::max(k_dt_min,k_dx_min) << ',' << std::min(k_dt_max,k_dx_max) << '\t';
            }
            if (verbose)
                std::cout << std::endl;
        }
    } else {
        for (unsigned i=0; i<_Kmax.h(); ++i){
            for (unsigned j=0; j<_Kmax.w(); ++j){
                const double dx2 = Dx(i);
                const double dt2 = Dt(j);
                
                const double T2t = T2+dt2;
                const double tt = mysin(T2t)*DX - mycos(T2t)*DY;
                
                const double denom = -(2.0*dx2*tt - DX2 - DY2 + dx_sq - dx2*dx2)/2.0;
                
                const double ttmdx2 = tt - dx2;
                const double k_dx_min = (ttmdx2-dx)/denom;
                //compute the intersection of all the surfaces to compute the final bundle
                if (_Kmin.val(i,j) < k_dx_min)
                    _Kmin.set_val(i,j,k_dx_min);
                
                const double k_dx_max = (ttmdx2+dx)/denom;
                //compute the intersection of all the surfaces to compute the final bundle
                if (_Kmax.val(i,j) > k_dx_max)
                    _Kmax.set_val(i,j,k_dx_max);
                
                const double sin_t2_p_dt2_m_t0_m_dt = mysin(T2t-T0_m_dt);
                const double sin_t2_p_dt2_m_t0_p_dt = mysin(T2t-T0_p_dt);
                const double sin_T0pdt_times_DY = DY*mysin(T0_p_dt);
                const double cos_T0pdt_times_DX = DX*mycos(T0_p_dt);
                const double sin_T0mdt_times_DY = DY*mysin(T0_m_dt);
                const double cos_T0mdt_times_DX = DX*mycos(T0_m_dt);
                
                const double k_dt_max = sin_t2_p_dt2_m_t0_m_dt/(-sin_t2_p_dt2_m_t0_m_dt*dx2 + cos_T0mdt_times_DX + sin_T0mdt_times_DY);
                //compute the intersection of all the surfaces to compute the final bundle
                if (_Kmax.val(i,j) > k_dt_max)
                    _Kmax.set_val(i,j,k_dt_max);
                
                const double k_dt_min = sin_t2_p_dt2_m_t0_p_dt/(-sin_t2_p_dt2_m_t0_p_dt*dx2 + cos_T0pdt_times_DX + sin_T0pdt_times_DY);
                //compute the intersection of all the surfaces to compute the final bundle
                if (_Kmin.val(i,j) < k_dt_min)
                    _Kmin.set_val(i,j,k_dt_min);
                if (verbose)
                    std::cout << std::max(k_dt_min,k_dx_min) << ',' << std::min(k_dt_max,k_dx_max) << '\t';
            }
            if (verbose)
                std::cout << std::endl;
        }
    }
}

#undef mysin
#undef mycos

//: Compute the best fit curve from the curve bundle
point_2d CC_curve_model_3d::compute_best_fit()
{
    if (this->bundle_is_valid() && !_is_poly_arc_formed)
    {
        //find the point of minimum deviation of the edgel
        unsigned mini = _Kmax.h(); unsigned minj = _Kmax.w();
        double mind = 100;
        //unsigned imax = 0, imin = _Kmax.h(), jmax = 0, jmin = _Kmax.w();
        for (unsigned i=0; i<_Kmax.h(); i++){
            for (unsigned j=0; j<_Kmax.w(); j++){
                
                if (_Kmax.val(i,j)>_Kmin.val(i,j)){
                    double d = Dx(i)*Dx(i)+Dt(j)*Dt(j);
//                    imin = imin>i?i:imin;
//                    jmin = jmin>j?j:jmin;
//                    imax = imax<i?i:imax;
//                    jmax = jmax<j?j:jmax;
                    if (d<mind){
                        mind = d;
                        mini=i;
                        minj=j;
                    }
                }
            }
        }
        
        assert(mini<_Kmax.h() && minj<_Kmax.w());
        
        //use this point to set the best fit of the curve bundle
        point_2d best_fit(Dx(mini), Dt(minj));
        
        this->set_best_fit(best_fit, (_Kmax.val(mini, minj)+_Kmin.val(mini, minj))/2);
	//this->set_best_fit(best_fit, _Kmin(mini, minj));        

    	// Ref PT are enforced to save K range at last !!!!
        _ref_pt.first = _Kmax.val(mini, minj);
        _ref_pt.second = _Kmin.val(mini, minj);
        
        //        std::cout << Dx(imin) << ' ' << Dx(imax) << ' ' << Dt(jmin) << ' ' << Dt(jmax) << ' ' << _Kmax.val(mini, minj) << ' ' << _Kmin.val(mini, minj) << ';' << std::endl;
        
        return best_fit;
    }
    
    return point_2d(0.0,0.0); //no optimal, just the ref
}


//: Set the best fit curve
void CC_curve_model_3d::set_best_fit(point_2d dx_dt, double kk)
{
    if(!_is_poly_arc_formed){
        double dx = dx_dt.first;
        double dt = dx_dt.second;
        
        //compute the extrinsic point and tangent of the centroid
        _theta = To2Pi(_ref_theta + dt);
        _pt = point_2d(_ref_pt.first-dx*sin(_theta), _ref_pt.second+dx*cos(_theta));
        _k= kk;
    }
}

//: update the best fit curve according to points on the curve
void CC_curve_model_3d::update_best_fit(point_2d p1, point_2d p2)
{
    
    // p1 is a point on the arc, p2 is the tail of the chain
    if (std::fabs(_k)<k_th) {
        // treat as a line
        _k = 0.0;
        _theta = std::atan2(p2.second-_pt.second, p2.first-_pt.first);
        // std::cout<<"line"<<std::endl;
        return;
    }
    double d1 = sq_norm(_pt), d2 = sq_norm(p1), d3 = sq_norm(p2);
    
    double a11 = 2.0*(_pt.first-p1.first), a12 = 2.0*(_pt.second-p1.second),
           a21 = 2.0*(p2.first-p1.first),  a22 = 2.0*(p2.second-p1.second);
    double b1 = d1 - d2, b2 = d3 - d2;
    double det = a11*a22 - a21*a12;
    
    // update center
    if (std::fabs(det)<k_th) {
        _k = 0.0;
        _theta = std::atan2(p2.second-_pt.second, p2.first-_pt.first);
        // std::cout<<"line"<<std::endl;
        return;
    }
    _center = point_2d((b1*a22-b2*a12)/det, (b2*a11-b1*a21)/det);
    
    point_2d v1(_pt.first-_center.first, _pt.second-_center.second), v2(p2.first-_center.first, p2.second-_center.second);
    
    // update angle
    double c = v1.first*v2.second-v1.second*v2.first, d = v1.first*v2.first+v1.second*v2.second;
    _angle = std::atan2(std::fabs(c), d);
    
    // update _theta
    _theta = std::atan2(c*v1.first, -c*v1.second);
    
    // update curvature
    _k = 1.0/sqrt(sq_norm(v2));
    
    _is_poly_arc_formed = true; // flag to update length

}

//: function to check if the curve fit is reasonable
bool CC_curve_model_3d::
curve_fit_is_reasonable(std::deque<edgel*> & edgel_chain, edgel* /*ref_e*/, double /*dpos*/)
{
    //compute LG-ratio and store as the quality
    compute_best_fit();
    return true;
}

//: are these two curve bundles C^2?
bool CC_curve_model_3d::is_C2_with(CC_curve_model_3d* cm)
{
    if (!cm)
        return true;

//    if (fabs(_k/cm->_k)>1.2 || fabs(_k/cm->_k)<0.8) {
//        return false;
//    }
    
    bool valid = false;
    for (unsigned i=0; i<_Kmax.h(); i++)
        for (unsigned j=0; j<_Kmax.w(); j++)
            valid = valid || (std::min(_Kmax.val(i,j), cm->_Kmax.val(i,j)) > std::max(_Kmin.val(i,j), cm->_Kmin.val(i,j)));
    
    return valid;
}

//: are these two curve bundles C^1?
bool CC_curve_model_3d::is_C1_with(CC_curve_model_3d* cm)
{
    if (!cm)
        return false;
    
    bool valid = false;
    for (unsigned i=0; i<_Kmax.h(); i++)
        for (unsigned j=0; j<_Kmax.w(); j++)
            valid = valid || ((_Kmax.val(i,j)>_Kmin.val(i,j)) && (cm->_Kmax.val(i,j)>cm->_Kmin.val(i,j)));
    
    return valid;
}

//: resize the curve model when the image is resized
void CC_curve_model_3d::resize(double scale)
{
    for (unsigned i=0; i<_Kmax.h(); i++){
        for (unsigned j=0; j<_Kmax.w(); j++){
            _Kmax.set_val(i, j, _Kmax.val(i, j)/scale);
            _Kmin.set_val(i, j, _Kmin.val(i, j)/scale);
        }
    }
    
    _ref_pt.first *= scale; _ref_pt.second *= scale;
    
    _pt.first *= scale; _pt.second *= scale;
    _center.first *= scale; _center.second *= scale;
    
    _k /= scale;
}

// //: report accuracy of measurement
// void CC_curve_model_3d::
// report_accuracy(double *estimates, double *min_estimates, double *max_estimates)
// {
//     double dt = To2Pi(_theta - _ref_theta);
//     double dx = sqrt((_ref_pt.first - _pt.first)*(_ref_pt.first - _pt.first) + (_pt.second - _ref_pt.second)* (_pt.second - _ref_pt.second));
//     _pt = point_2d(_ref_pt.first-dx*sin(_theta), _ref_pt.second+dx*cos(_theta));
//     _k= kk;
// double theta_min=1000.0, theta_max=-1000.0, k_min=1000.0, k_max=-1000.0;
// 
// //FIX ME
// 
// //report these numbers
// estimates[0] = theta;  //theta
// estimates[1] = k;      //curvature
// estimates[2] = 0.0;    //curvature derivative
// 
// min_estimates[0] = dbdet_angle0To2Pi(theta_min+ref_theta);
// min_estimates[1] = k_min;
// min_estimates[2] = 0.0;
// 
// max_estimates[0] = dbdet_angle0To2Pi(theta_max+ref_theta);
// max_estimates[1] = k_max;
// max_estimates[2] = 0.0;
// 
// }


void CC_curve_model_3d::set_output(arrayd &curvelet_info, unsigned pos)
{
    curvelet_info.set_val(pos, 1, _ref_pt.first);
    curvelet_info.set_val(pos, 2, _ref_pt.second);
    curvelet_info.set_val(pos, 3, _ref_theta);
    curvelet_info.set_val(pos, 4, _pt.first);
    curvelet_info.set_val(pos, 5, _pt.second);
    curvelet_info.set_val(pos, 6, _theta);
    curvelet_info.set_val(pos, 7, _k);
//    curvelet_info.set_val(pos, 7, _angle);
    if(_is_poly_arc_formed){
        curvelet_info.set_val(pos, 10, _center.first);
        curvelet_info.set_val(pos, 11, _center.second);
    }
}

//: print central info to file
void CC_curve_model_3d::print(std::ostream&  os)
{
    //FIX ME
    os << "[";
    os << _ref_pt.first<<" ";
    // os << "] [";
    os << _ref_pt.second<<" ";
    //os << "] [";
    os << _ref_theta<<" ";
    //os << "] [";
    os << _pt.first<<" ";
    //os << "] [";
    os << _pt.second<<" ";
    //os << "] [";
    os << _theta<<" ";
    //os << "] [";
    os << _k;
    
    //for (int i=0; i<NkClasses; i++){
    //  os << "<";
    //  if (cv_bundles[i].num_sheets()==1){
    //    for (unsigned p=0; p<cv_bundles[i][0].size(); p++)
    //      os << "(" << cv_bundles[i][0][p].x() << " " << cv_bundles[i][0][p].y() << ")";
    //  }
    //  os << ">";
    //}
    os << "]";
}

//: print bundle info
void CC_curve_model_3d::print_bundle( std::ostream&  os)
{
    //if the Kmax surface is completely less than the Kmin surface, there is no bundle
    // save the current settings
    std::ios::fmtflags old_settings = os.flags(); //save previous format flags
    int old_precision = os.precision(); // save previous precision setting
    os.precision(4);
    os.setf(std::ios::fixed, std::ios::floatfield);
    
    os << '\n';
    for (unsigned i=0; i<_Kmax.h(); ++i){
        if (i==0) {
            os << "\t\t";
            for (unsigned j=0; j<_Kmax.w(); ++j){
                os << Dt(j) << "\t\t\t";
            }
            os << '\n';
        }
        for (unsigned j=0; j<_Kmax.w(); ++j){
            if(j==0){
                os << Dx(i) << "\t\t";
            }
            if (_Kmax.val(i,j) > _Kmin.val(i,j))
                os << _Kmin.val(i,j) << ',' << _Kmax.val(i,j) << "\t\t";
            else
                os << "*************" <<"\t\t";
        }
        os << '\n';
    }
    os<<std::endl;
    os.precision(old_precision);
    os.flags(old_settings);
}
