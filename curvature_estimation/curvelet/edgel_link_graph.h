// This is edgel_link_graph.h
#ifndef edgel_link_graph_h
#define edgel_link_graph_h

/*****************************************************************************
// file: edgel_link_graph.h
// brief: Link graph structure
// author: Xiaoyan Li
// date: 02/04/2016
******************************************************************************/


#include <vector>
#include <list>
#include <iostream>

#include "edgemap.h"
#include "curvelet.h"

//: structure to describe an edgel link in the link graph
class edgel_link {
    
public:
    
    edgel* _pe;                  ///< parent edgel
    edgel* _ce;                  ///< child edgel
    
    double _vote;                      ///< record votes for each link
    
    int _deg_overlap;                  ///< to keep track of the degree of overlap
//    bool _flag;                        ///< a falg for miscellaneous purposes
    
    std::list<curvelet*> _curvelets; ///< list of curvelets containing this link
    
    //: constructor
    edgel_link(edgel* e1, edgel* e2): _pe(e1), _ce(e2), _vote(0.0), _deg_overlap(0), /*_flag(false),*/ _curvelets(0){}
    
    //: count the degree of overlap between pairs of curvelets at a link
    void set_degree_overlap()
    {
        int _deg_overlap = 0;
        
        //go over all pairs of curvelets causing the link
        cvlet_list::iterator c_it = _curvelets.begin();
        for (; c_it != _curvelets.end(); c_it++){
            cvlet_list::iterator c_it2 = c_it;
            c_it2++;
            for (; c_it2 != _curvelets.end(); c_it2++)
            {
                unsigned k1=0, k2=0;
                
                //count the degree of overlap between these two curvelets
                for (unsigned i=0; i<(*c_it)->_edgel_chain.size(); i++){
                    if ((*c_it)->_edgel_chain[i]==_pe){
                        k1 = i;
                        break;
                    }
                }
                
                for (unsigned i=0; i<(*c_it2)->_edgel_chain.size(); i++){
                    if ((*c_it2)->_edgel_chain[i]==_pe){
                        k2 = i;
                        break;
                    }
                }
                
                int deg = count_degree_overlap((*c_it), (*c_it2), k1, k2);
                if (deg>_deg_overlap)
                    _deg_overlap = deg;
            }
        }
    }
    
    //: count the degree of overlap between two curvelets
    int count_degree_overlap(curvelet* cvlet1, curvelet* cvlet2, unsigned k1, unsigned k2)
    {
        //ks are the indices into the edgel chain of the cvlets
        int cnt_overlap=0;
        
        //count the overlaps before the link
        int kk1=k1-1;
        int kk2=k2-1;
        bool continuous = true;
        for (; kk1>=0 && kk2>=0 && continuous; kk1--, kk2--){
            if (cvlet1->_edgel_chain[kk1]==cvlet2->_edgel_chain[kk2])
                cnt_overlap++;
            else
                continuous = false;
        }
        
        //count the overlaps after the link
        kk1=k1+2; kk2=k2+2;
        continuous = true;
        for (; kk1<(int)cvlet1->_edgel_chain.size() && kk2<(int)cvlet2->_edgel_chain.size() && continuous; kk1++, kk2++){
            if (cvlet1->_edgel_chain[kk1]==cvlet2->_edgel_chain[kk2])
                cnt_overlap++;
            else
                continuous = false;
        }
        return cnt_overlap;
    }
    
    
//    //: prune redundant curvelets
//    void prune_redundant_curvelets()
//    {
//        //go over all the curvelets
//        cvlet_list::iterator cv_it = _curvelets.begin();
//        for (; cv_it != _curvelets.end(); cv_it++)
//        {
//            curvelet* cvlet = (*cv_it);
//            
//            cvlet_list::iterator cv_it2 = cv_it;
//            cv_it2++;
//            for (; cv_it2 != _curvelets.end(); cv_it2++)
//            {
//                curvelet* dup_cvlet = (*cv_it2);
//                
//                //determine if this is duplicate
//                if (cvlet->_edgel_chain.size() != dup_cvlet->_edgel_chain.size())
//                    continue;
//                
//                bool duplicate = true; //reset flag
//                for (unsigned k=0; k<cvlet->_edgel_chain.size(); k++)
//                    duplicate = duplicate && (cvlet->_edgel_chain[k]==dup_cvlet->_edgel_chain[k]);
//                
//                if (duplicate){
//                    cv_it2--;
//                    _curvelets.remove(dup_cvlet);
//                }
//            }
//        }
//    }
    
};

typedef std::list<edgel_link*> link_list;

//: This class represents the link graph formed from the edgels
//  The links are described by an adjacency list organized by the id of the
//  edgel
class edgel_link_graph
{
public:
    std::vector<link_list> _cLinks; ///< child links
    std::vector<link_list> _pLinks; ///< parent links
    std::vector<bool> _linked;      ///< flag to signal whether this edgel has been linked(in CFG) already
    
    //temporary link graph created for finding selected edgel chains
    std::vector<link_list> _cLinks2; ///< child links
    std::vector<link_list> _pLinks2; ///< parent links
    
    //: constructor
    edgel_link_graph(int size=0): _cLinks(size), _pLinks(size), _linked(size), _cLinks2(size), _pLinks2(size){}
    
    //: destructor
    ~edgel_link_graph()
    {
        clear();//delete everything upon exit
    }
    
    //: resize the graph
    void resize(unsigned size)
    {
        if (size!=_cLinks.size())
            clear();
        
        _cLinks.resize(size);
        _pLinks.resize(size);
        _linked.resize(size);
        
        _cLinks2.resize(size);
        _pLinks2.resize(size);
    }
    
    //: clear the graph
    void clear()
    {
        //delete all the links
        for (unsigned i=0; i<_cLinks.size(); i++){
            link_list::iterator l_it = _cLinks[i].begin();
            for (;l_it!=_cLinks[i].end(); l_it++)
                delete (*l_it);
        }
        
        _cLinks.clear();
        _pLinks.clear();
        _linked.clear();
        
        _cLinks2.clear();
        _pLinks2.clear();
    }
    
    void clear_linked_flag()
    {
        _linked.clear();
        _linked.resize(_cLinks.size());
    }
    //
    //  //utility functions
    //  bool is_an_end(int id)  {return _cLinks[id].size()==0 || _pLinks[id].size()==0; }
    //  bool is_a_bifurc(int id){return _cLinks[id].size()==2;}
    //  bool is_a_split(int id) {return _cLinks[id].size()>1;}
    //  bool is_a_merge(int id) {return _pLinks[id].size()>1;}
    //  bool is_a_junct(int id) {return _cLinks[id].size()>1 || _pLinks[id].size()>1;} //X-junction

    //: is this edgel legal to start a one-chain
    bool edgel_is_legal_first_edgel(edgel* eA)
    {
        return (!_linked[eA->_id] &&            //not yet part of any edgel chain
                (_pLinks[eA->_id].size()==1 &&   //one link before it and
                 _cLinks[eA->_id].size()==1) &&  //one link after it
                 _pLinks[eA->_id].front()->_pe != _cLinks[eA->_id].front()->_ce //not a loop
                );
    }

    //: is this edgel legal to be on an edgel chain
    bool edgel_is_legal(edgel* eA)
    {
        return (!_linked[eA->_id] &&            //not yet part of any edgel chain
                _pLinks[eA->_id].size()<=1 &&   //one link before it (or endpoint) and
                _cLinks[eA->_id].size()<=1      //one link after it (or endpoint)
                );
    }
    
    //: check if edgels are linked,
    //  if they are linked, return the link
    edgel_link* are_linked(edgel* e1, edgel* e2)
    {
        link_list::iterator cit = _cLinks[e1->_id].begin();
        for(; cit!=_cLinks[e1->_id].end(); cit++){
            if ((*cit)->_ce == e2)
                return (*cit); //return link
        }
        link_list::iterator cit1 = _cLinks[e2->_id].begin();
        for(; cit1!=_cLinks[e2->_id].end(); cit1++){
            if ((*cit1)->_ce == e1)
                return (*cit1); //return link
        }
        return 0; //not _inked
    }
    
    //: link edgels if not already linked
    void link_edgels(edgel* e1, edgel* e2, curvelet* cvlet)
    {
        edgel_link* alink = are_linked(e1, e2);
        if (!alink){
            //create a link
            edgel_link* new_link = new edgel_link(e1, e2);
            
            _cLinks[e1->_id].push_back(new_link);
            _pLinks[e2->_id].push_back(new_link);
            
            //also add a vote
            new_link->_vote += 1.0;
            new_link->_curvelets.push_back(cvlet);
        }
        else {
            //link already exists so simply add to the vote count
            alink->_vote += 1.0;
            
            //and add the cvlet that gave rise to this link
            alink->_curvelets.push_back(cvlet);
        }
    }
    
    //  //: remove the link between these edgels
    //  void remove_link(link* alink)
    //  {
    //    _pLinks[alink->ce->id].remove(link);
    //    _cLinks[alink->pe->id].remove(link);
    //
    //    delete link;
    //  }
    //
    //  //: move a link from the first link graph to the second
    //  void move_link(dbdet_link* link)
    //  {
    //    _pLinks[link->ce->id].remove(link);
    //    _cLinks[link->pe->id].remove(link);
    //
    //    _pLinks2[link->ce->id].push_back(link);
    //    _cLinks2[link->pe->id].push_back(link);
    //  }
    
    //: print the graph
    void print(std::ostream& os)
    {
        for(int i = 0; i<_linked.size(); i++)
        {
            os << '[' << i << "]\t child links: ";
            for (link_list::iterator lit=_cLinks[i].begin(); lit!=_cLinks[i].end(); lit++) {
                os << (*lit)->_ce->_id << ' ';
            }
            os << std::endl;
            
            os << "\t parent links: ";
            for (link_list::iterator lit=_pLinks[i].begin(); lit!=_pLinks[i].end(); lit++) {
                os << (*lit)->_pe->_id << ' ';
            }
            os << std::endl;
        }
    }
};

#endif // dbdet_sel_base_h
