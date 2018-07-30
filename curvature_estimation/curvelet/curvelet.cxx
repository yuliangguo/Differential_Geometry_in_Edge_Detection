#include <iostream>
#include <cmath>
#include "curvelet.h"

extern double sq_distance(const point_2d& A, const point_2d& B);
extern const double k_th;

// connect two edge chain
edgel_chain connect(const edgel_chain &first, const edgel_chain &second)
{
    edgel_chain new_edgel_chain(first);
    edgel_chain::iterator fit = new_edgel_chain.begin(), fit2 = new_edgel_chain.begin();
    edgel_chain::const_iterator sit;
//    std::cout << "first:";
//    for (; fit2!=new_edgel_chain.end(); fit2++) {
//        std::cout << (*fit2)->_id <<' ';
//    }
//    std::cout << std::endl;
//    std::cout << "second:";
//    for (sit=second.begin(); sit!=second.end(); sit++) {
//        std::cout << (*sit)->_id <<' ';
//    }
//    std::cout << std::endl;
    while (fit!=new_edgel_chain.end()) {
        // find the common element in second
        sit = find(second.begin(),second.end(),*fit);
        if (sit!=second.end()) {
            fit2 = find(new_edgel_chain.begin(),new_edgel_chain.end(),*sit);
            // find the element only in second
            // fit point to the last common element in first
            while (fit2!=new_edgel_chain.end() && sit!=second.end()) {
                sit++;
                fit = fit2;
                fit2 = find(new_edgel_chain.begin(),new_edgel_chain.end(),*sit);
            }
            if (sit!=second.end()) {
                new_edgel_chain.insert(fit+1, *sit);
            }
            else
                break;
        }
        fit++;
    }
//    std::cout << "final:";
//    for (fit2=new_edgel_chain.begin(); fit2!=new_edgel_chain.end(); fit2++) {
//        std::cout << (*fit2)->_id <<' ';
//    }
//    std::cout << std::endl;
    return new_edgel_chain;
}

//: copy constructor
curvelet::curvelet(const curvelet& other)
{
    //the edgels have to copied as links because curvelets are just groupings of the edgels
    _ref_edgel = other._ref_edgel;
    _edgel_chain = other._edgel_chain;
    
    //but the curve model has to be deep copied
    _curve_model = new CC_curve_model_3d(*(other._curve_model));
    
    _forward = other._forward;
    _length = other._length;
    _quality = other._quality;
    _used = other._used;
}

curvelet::~curvelet()
{
    if(_curve_model)
        delete _curve_model;
    _edgel_chain.clear();
}

// weighting constants for the heuristic
#define alpha3 1.0
#define alpha4 1.0

//: compute properties of this curvelet once formed
void curvelet::compute_properties(double R, double token_len)
{
  //find out the # of edgels before and after the reference edgel
  //also find out the length before and after the reference edgel
  int num_before=0, num_after=0;
  double Lm=0, Lp=0;

  bool before_ref = true;
  for (unsigned i=0; i<_edgel_chain.size()-1; i++){
    if (before_ref) { Lm += sqrt(sq_distance(_edgel_chain[i]->_pt, _edgel_chain[i+1]->_pt)); num_before++; }
    else            { Lp += sqrt(sq_distance(_edgel_chain[i]->_pt, _edgel_chain[i+1]->_pt)); num_after++; }

    if (_edgel_chain[i+1]==_ref_edgel)
      before_ref = false;
  }

  //compute the length of the curvelet (extrinsic length)
  _length = Lm+Lp;

  //also compute the LG ratio and store as quality
  //quality = (num_before+num_after)*token_len/length;

  //new quality measure (1/cost of the compatibility heauristic)
  _quality = 2/(alpha3*R/_length + alpha4*_length/token_len/_edgel_chain.size());
}

//: print info to file
void curvelet::print(std::ostream& os)
{
    // anchor point
    os << "[" << _ref_edgel->_id << "] ";
  
  //first output the edgel chain
  os << "[";
  for (unsigned i=0; i< _edgel_chain.size(); i++){
    os << _edgel_chain[i]->_id << " ";
  }
  os << "] ";

  //forward/backward tag
  os << "(";
  if (_forward) os << "F";
  else         os << "B";
  os << ") ";

  //next output the curve model
  _curve_model->print(os);

  //then output the other properties
  os << " " << _length << " " << _quality << std::endl;
}

void curvelet::set_output(arrayi &id_chain, arrayd &curvelet_info, unsigned posy, unsigned posx)
{
    assert(_edgel_chain.size()<id_chain.w());
    for (unsigned i=0; i< _edgel_chain.size(); i++){
        // convert to Matlab index
        id_chain.set_val(posy, i+posx, (_edgel_chain[i]->_id)+1);
    }
    
    _curve_model->set_output(curvelet_info, posy);
    curvelet_info.set_val(posy, 0, _forward?1.0:0.0);
    curvelet_info.set_val(posy, 8, _length);
    curvelet_info.set_val(posy, 9, _quality);
}

curvelet* curvelet::intersect(curvelet* c2)
{
    // NOTE: Assume c2 comes after c1
    // get the intersection of curve model
    CC_curve_model_3d* cur_cm  = _curve_model;
    cur_cm  = cur_cm->intersect(c2->_curve_model);
    
    // get the intersection of edgel chain
    edgel_chain cur_edgel_chain = connect(_edgel_chain,c2->_edgel_chain);
    
    edgel * ref_e = *(_edgel_chain.begin());
    
    return new curvelet(ref_e, cur_cm, cur_edgel_chain, _forward);
    
}

void curvelet::update_length()
{
    if (std::fabs(_curve_model->_angle)>0 && std::fabs(_curve_model->_k)>k_th
        && _curve_model->_is_poly_arc_formed) {
        // if length<0, arc is in the opposite direction
        _length = (_curve_model->_angle)/std::fabs(_curve_model->_k);
    }
}
