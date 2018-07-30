#ifndef curve_frag_topo_graph_h_
#define curve_frag_topo_graph_h_

/*****************************************************************************
 // file: curveletmap.h
 // brief: Curve Fragment Topology Graph (CFTG) structure
 // Reference: lemsvxl code by Amir Tamrakar
 // author: Xiaoyan Li
 // date: 02/03/2016
 ******************************************************************************/

#include <vector>
#include <list>

#include "edgemap.h"

class CFTG_link
{
public:
    edgel* _eS;
    edgel* _eE;
    edgel_chain_list _cCFs; ///< candidate curve fragments
    double _cost;
    
    CFTG_link(edgel* e1, edgel* e2): _eS(e1), _eE(e2), _cCFs(0), _cost(0.0){}
    ~CFTG_link(){ _eS = 0; _eE = 0; _cCFs.clear(); _cost = 0;}
    
    void add_fragment(edgel_chain* chain)
    {
        _cCFs.push_back(chain);
    }
    
};

typedef std::list<CFTG_link*> CFTG_link_list;

//: This class represents the curve fragment topology graph formed from the contour fragments
//  The links are equivalence classes of curve fragments represented by sets of edgel chains
class curve_frag_topo_graph
{
public:
    std::vector<CFTG_link_list> _cLinks; ///< child links
    std::vector<CFTG_link_list> _pLinks; ///< parent links
    
    CFTG_link_list _Links; ///< redundant single list of all links
    
    //: constructor
    curve_frag_topo_graph(int size=0): _cLinks(size), _pLinks(size){}
    
    //: destructor
    ~curve_frag_topo_graph()
    {
        clear(); //delete everything upon exit
    }
    
    //Access functions
    unsigned size() { return _cLinks.size(); }//should be the same as edgels.size()
    
    //: resize the graph
    void resize(unsigned size)
    {
        if (size!=_cLinks.size())
            clear();
        
        _cLinks.resize(size);
        _pLinks.resize(size);
    }
    
    //: clear the graph
    void clear()
    {
        //delete all the curve fragments
        CFTG_link_list::iterator l_it = _Links.begin();
        for (; l_it != _Links.end(); l_it++)
            delete (*l_it);
        
        _Links.clear();
        _cLinks.clear();
        _pLinks.clear();
    }
    
    //: add a curve fragment to the graph
    void insert_fragment(edgel_chain* chain)
    {
        edgel* e1 = chain->front();
        edgel* e2 = chain->back();
        
        //if there is a link already, add it to the existing link
        CFTG_link* cur_Link = 0;
        if (_cLinks[e1->_id].size()>0){
            //find the link in here
            CFTG_link_list::iterator l_it = _cLinks[e1->_id].begin();
            for (; l_it != _cLinks[e1->_id].end(); l_it++){
                if ((*l_it)->_eE==e2){
                    cur_Link = (*l_it);
                    break;
                }
            }
        }
        
        if (cur_Link){
            cur_Link->add_fragment(chain);
        }
        else //otherwise create a new link and add it to it
        {
            cur_Link = new CFTG_link(e1, e2);
            cur_Link->add_fragment(chain);
            
            //add the link to the graph
            _cLinks[e1->_id].push_back(cur_Link);
            _pLinks[e2->_id].push_back(cur_Link);
            
            _Links.push_back(cur_Link);
        }
    }
    
    
    //: remove a Link
    void remove_link(CFTG_link* link)
    {
        _cLinks[link->_eS->_id].remove(link);
        _pLinks[link->_eE->_id].remove(link);
        
        _Links.remove(link);
        
        delete link;
    }
    
};

#endif // curve_frag_topo_graph_h_
