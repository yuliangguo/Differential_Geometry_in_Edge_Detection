// This is ployarcmap.h
#ifndef polyarcmap_h_
#define polyarcmap_h_

/*****************************************************************************
 // file: ployarcmap.h
 // brief: polyarc map data structure
 // author: Xiaoyan Li
 // date: 02/05/2016
 ******************************************************************************/

#include <algorithm>

#include "edgemap.h"
#include "curvelet.h"

extern point_2d rotate(const point_2d &pt, double theta);
extern double vPointPoint (point_2d& start, point_2d& end);
extern double dot (double v1, double v2);
extern const double k_th;

//: This class stores the map of poly arcs formed by the SEL edge linker.
class polyarcmap
{
public:
    
    //: constructor
    polyarcmap(edgemap* EM=0): _EM(EM), _samples(NULL){};
    
    //: destructor
    ~polyarcmap(){
        clear(); //delete everything upon exit
        _EM=0; //delete the reference to the edgemap
    };
    
    //: to check if the Curvelet map is valid
    bool is_valid() { return _list.size()>0 && _list.size()==_EM->num_edgels(); };
    
    //: add a curvelet to an edgel
    void add_curvelet(curvelet* cvlet)
    {
        _list.push_back(cvlet);
        
    }
    
    //: does this curvelet exist at this edgel?
    curvelet* does_curvelet_exist(edgel_chain & chain)
    {
        //go over all the curvelets of the current size formed by the current edgel
        cvlet_list::iterator cv_it;
        for ( cv_it = _list.begin(); cv_it != _list.end(); cv_it++){
            curvelet* cvlet = (*cv_it);
            
            if (cvlet->_edgel_chain.size() != chain.size())
                continue;
            
            bool cvlet_exists = true; //reset flag
            for (unsigned k=0; k<chain.size(); k++)
                cvlet_exists = cvlet_exists && (cvlet->_edgel_chain[k]==chain[k]);
            
            //the flag will remain true only if all the edgels match
            if (cvlet_exists)
                return cvlet; //return matching curvelet
        }
        return 0; //curvelet does not exist
    }
    
    //: return the number of curvelets in this curvelet map
    unsigned num_curvelets()
    {
        return _list.size();
    }
    
    //: record the curvelet map data in array format
    void get_output_array(arrayi &id_chain, arrayd &curvelet_info)
    {
        int count = 0;
        for (cvlet_list::iterator it=_list.begin(); it!=_list.end(); it++) {
            // Note: no anchor edgel
            (*it)->set_output(id_chain, curvelet_info, count, 0);
            count++;
        }
    }
    
    //: return the maximum length of edgel chain in this map
    unsigned edgel_chain_max_length()
    {
        size_t max_len = 0;
        cvlet_list::iterator cit;
        for(cit= _list.begin(); cit!=_list.end(); cit++)
        {
            max_len = std::max(max_len, ((*cit)->_edgel_chain).size());
        }
        return unsigned(max_len);
    }
    
    //: end points pose correction
    void compute_best_fit_according_to_curvelets()
    {
        cvlet_list::iterator cit;
        unsigned i;
        for(cit=_list.begin(), i=0; cit!=_list.end(); cit++, i++ ){
            edgel* eA = (*cit)->_edgel_chain.front();
            (*cit)->_curve_model->_pt = eA->_pt;
            
            unsigned ind = (*cit)->_edgel_chain.size()/2;
            edgel* center_pt = ((*cit)->_edgel_chain)[ind];
            (*cit)->_curve_model->update_best_fit(center_pt->_pt, (*cit)->_edgel_chain.back()->_pt);
            (*cit)->update_length();
        }
    }
    
    void sample_from_poly_arc()
    {
        if(!_samples){
            unsigned sampledensity = 2;
            
            _samples = new edgemap(_EM->width(), _EM->height());
            
            cvlet_list::iterator cit;
            unsigned i;
            
            unsigned edgeid = 0;
            
            for(cit=_list.begin(), i=0; cit!=_list.end(); cit++, i++ ){
                //if((*cit)->_curve_model->_is_poly_arc_formed){
                    double k = (*cit)->_curve_model->_k;
                    double theta = (*cit)->_curve_model->_theta;
                    point_2d pt((*cit)->_curve_model->_pt);
                    double length = (*cit)->_length;
                    unsigned numsample = sampledensity*fabs(length);
                    numsample = numsample<1?1:numsample;
                    
                    if(fabs(k)<k_th){
                        double detalen = length/numsample;
                        
                        for(unsigned j=0; j<numsample+1; j++){
                            double x = pt.first + cos(theta)*j*detalen;
                            double y = pt.second + sin(theta)*j*detalen;
                            
                            _samples->insert(x,y,theta,2.0,edgeid);
                            _map[edgeid] = i;
                            //std::cout<<x<<','<<y<<','<<theta<<','<<edgeid<<';'<<std::endl;
                            edgeid++;
                        }
                    }
                    else{
                        double radius = fabs(1.0/k);
                        point_2d center((*cit)->_curve_model->_center);
                        point_2d v1(pt.first - center.first, pt.second - center.second);
                        double angle = (*cit)->_curve_model->_angle;
                        double costheta = cos(theta), sintheta = sin(theta);
                        
                        double detangle = angle/numsample;
                        
                        double c = v1.second*(v1.second*cos(angle)+radius*sintheta*sin(angle))-
                        v1.second*(v1.first*cos(angle)+radius*costheta*sin(angle));
                        
                        for(unsigned j=0; j<numsample+1; j++){
                            double x = v1.first*cos(detangle*j)+radius*costheta*sin(detangle*j)+center.first;
                            double y = v1.second*cos(detangle*j)+radius*sintheta*sin(detangle*j)+center.second;
                            double orientation = atan2(c*x, -c*y);
                            
                            _samples->insert(x,y,orientation,2.0,edgeid);
                            _map[edgeid] = i;
                            //std::cout<<x<<','<<y<<','<<orientation<<','<<edgeid<<';'<<std::endl;
                            edgeid++;
                        }
                    }
                //}
            }
            //_samples->print(std::cout);
        }
    }
    
    //  resize the polyarc map
    void resize(double scale)
    {
        for (cvlet_list::iterator it=_list.begin(); it!=_list.end(); it++) {
            (*it)->_length *= scale;
            (*it)->_curve_model->resize(scale);
        }
        if(_samples){
            delete _samples;
            sample_from_poly_arc();
        }
    }
    
    
    //: print info to file
    void print(std::ostream& os)
    {
        for (cvlet_list::iterator it=_list.begin(); it!=_list.end(); it++) {
            (*it)->print(os);
        }
    }
    
private:
    
    //: The edgemap on which these curvelets have been formed
    //  (due to this smart pointer to the edgemap, the curvelets remain valid even if the edgemap is deleted elsewhere)
    edgemap* _EM;
    
    //: The polyarc map
    cvlet_list _list;
    
    //: edgemap formed by samples of polyarc map
    edgemap* _samples;
    
    //: map edge id to polyarc id
    std::map<int,int> _map;
    
    //: clear the graph
    void clear()
    {
        std::set<curvelet*> deleted;
        
        cvlet_list::iterator p_it;
        for (p_it = _list.begin(); p_it != _list.end(); p_it++) {
            delete (*p_it);
        }
        _list.clear();
        
        if(_samples) delete _samples;
        _map.clear();
    }
    
    friend class edgemap;
};

#endif // curveletmap_h_
