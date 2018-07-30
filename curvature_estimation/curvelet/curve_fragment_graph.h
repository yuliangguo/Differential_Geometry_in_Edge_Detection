// This is curve_fragment_graph.h
#ifndef curve_fragment_graph_h
#define curve_fragment_graph_h

/*****************************************************************************
 // file: curveletmap.h
 // brief: Curve Fragment graph structure
 // Reference: lemsvxl code by Amir Tamrakar
 // author: Xiaoyan Li
 // date: 02/03/2016
 ******************************************************************************/

#include <vector>
#include <list>
#include <set>
#include "edgemap.h"
#include "curve_frag_topo_graph.h"

//: This class represents the curve fragment graph formed from the edgels
//  The links are curve fragments represented by edgel chains.
class curve_fragment_graph
{
public:
    std::vector<edgel_chain_list> _cFrags; ///< child curve fragments
    std::vector<edgel_chain_list> _pFrags; ///< parent curve fragments
    
    edgel_chain_list _frags; ///< redundant single list of all fragments
    
    //this is a hack (need to move this out of here into the storage class)
    curve_frag_topo_graph _CFTG; ///< The Curve Fragment Topology Graph (CFTG)
    
    std::set<int> _participate_edge_id;  // include edges in unambiguous fragments and edges participate in hypothesis tree;
    
    //: constructor
    curve_fragment_graph(int size=0): _cFrags(size), _pFrags(size){}
    
    //: destructor
    ~curve_fragment_graph()
    {
        clear(); //delete everything upon exit
        _CFTG.clear();
    }
    
    //Access functions
    unsigned size() { return _cFrags.size(); }//should be the same as edgels.size()
    
    //: resize the graph
    void resize(unsigned size)
    {
        if (size!=_cFrags.size()){
            clear();
            _CFTG.clear();
        }
        
        _cFrags.resize(size);
        _pFrags.resize(size);
        
        _CFTG.resize(size);
    }
    
    //: clear the graph
    void clear()
    {
        //delete all the curve fragments
        edgel_chain_list::iterator f_it = _frags.begin();
        for (; f_it != _frags.end(); f_it++)
            delete (*f_it);
        
        _frags.clear();
        _cFrags.clear();
        _pFrags.clear();
        
        _CFTG.clear();
    }
    
    //: add a curve fragment to the graph
    void insert_fragment(edgel_chain* chain)
    {
        edgel* e1 = chain->front();
        edgel* e2 = chain->back();
        
        _cFrags[e1->_id].push_back(chain);
        _pFrags[e2->_id].push_back(chain);
        
        _frags.push_back(chain);
    }
    
    //: remove a curve fragment
    void remove_fragment(edgel_chain* chain)
    {
        edgel* e1 = chain->front();
        edgel* e2 = chain->back();
        
        _pFrags[e2->_id].remove(chain);
        _cFrags[e1->_id].remove(chain);
        
        _frags.remove(chain);
        
        delete chain;
    }
    
    //: just extract a curve fragment from the graph do not delete
    void extract_fragment(edgel_chain* chain)
    {
        edgel* e1 = chain->front();
        edgel* e2 = chain->back();
        
        _pFrags[e2->_id].remove(chain);
        _cFrags[e1->_id].remove(chain);
        
        _frags.remove(chain);
    }
    
    //: return the maximum length of edgel chain in this map
    unsigned edgel_chain_max_length()
    {
        size_t max_len = 0;
        edgel_chain_list::iterator eit;
        for(eit=_frags.begin(); eit!=_frags.end(); eit++)
        {
            max_len = std::max(max_len, (*eit)->size());
        }
        return unsigned(max_len);
    }
    
    //: return the number of curve fragments in this graph
    unsigned num_frags()
    {
        return unsigned(_frags.size());
    }
    
    //: record the curvelet map data in array format
    void get_output_array(arrayi &id_chain, arrayd &/*curvelet_info*/)
    {
        for (edgel_chain_list::iterator lit=_frags.begin(); lit!=_frags.end(); lit++) {
            for (edgel_chain::iterator eit=(*lit)->begin(); eit!=(*lit)->end(); eit++) {
                // anchor edgel, convert to matlab index
                id_chain.set_val(unsigned(std::distance(_frags.begin(), lit)), unsigned(eit-(*lit)->begin()), int((*eit)->_id+1));
            }
        }
        
    }
    
    //: print the graph
    void print(std::ostream& os)
    {
        edgel_chain_list::iterator lit;
        for(lit = _frags.begin(); lit!=_frags.end(); lit++)
        {
            os << "[ ";
            for (edgel_chain::iterator cit=(*lit)->begin(); cit!=(*lit)->end(); cit++) {
                os << (*cit)->_id << ' ';
            }
            os << ']' << std::endl;
        }
    }
    friend class edgemap;
};

#endif // dbdet_curve_fragment_graph_h
