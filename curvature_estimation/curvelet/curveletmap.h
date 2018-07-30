// This is curveletmap.h
#ifndef curveletmap_h_
#define curveletmap_h_

/*****************************************************************************
 // file: curveletmap.h
 // brief: Curvelet Map data structure
 // author: Xiaoyan Li
 // date: 01/26/2016
 ******************************************************************************/

#include "edgemap.h"
#include "curvelet.h"

//: This class stores the map of curvelets formed by the SEL edge linker.
//  This is an intermediate structure before actual edge linking occurs.
class curveletmap
{
public:
    
    //: constructor
    curveletmap(edgemap* EM=0);
    
    //: destructor
    ~curveletmap();
    
    //: to check if the Curvelet map is valid
    bool is_valid() { return _map.size()>0 && _map.size()==_EM->num_edgels(); };
    
    //: set the edgemap
    //void set_edgemap(edgemap* EM) { _EM = EM; resize(EM->num_edgels()); }
    
    //: access the curvelets for an edge using its id
    const cvlet_list& getcurvelets(unsigned id) const { return _map[id]; };
    cvlet_list& getcurvelets(unsigned id) { return _map[id]; };
    
    //: add a curvelet to an edgel
    void add_curvelet(curvelet* curvelet);
    
    //: remove a curvelet from this edgel
    // void remove_curvelet(dbdet_curvelet* curvelet);
    
    //: delete all the curvelets formed by this edgel
    // void delete_all_curvelets(dbdet_edgel* e);
    
    //: does this curvelet exist at this edgel?
    curvelet* does_curvelet_exist(edgel* e, std::deque<edgel*> & chain);
    
    //: return the number of curvelets in this curvelet map
    unsigned num_curvelets();
    
    //: record the curvelet map data in array format
    void get_output_array(arrayi &id_chain, arrayd &curvelet_info);
    
    // forms a curvelet map which is a mapping from each edgel to all curvelets it participates in
    void form_full_cvlet_map();
    
    //: print info to file
    void print(std::ostream&);
    
private:
    
    //: The edgemap on which these curvelets have been formed
    //  (due to this smart pointer to the edgemap, the curvelets remain valid even if the edgemap is deleted elsewhere)
    edgemap* _EM;
    
    //: The curvelet map, indexed by edgel IDs
    std::vector<cvlet_list > _map;
    
    bool _fullcvmap;
    
    //: resize the graph
    void resize(unsigned size);
    
    //: clear the graph
    void clear();//!!!!!!!!
    
    friend class edgemap;
};

#endif // curveletmap_h_
