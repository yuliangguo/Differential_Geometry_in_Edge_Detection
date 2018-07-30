// This is edgemap.h
#ifndef edgemap_h
#define edgemap_h

/***************************************************************************************
//file: edgemap.h
//brief: The edge map class
//author: Xiaoyan Li
***************************************************************************************/

#include <map>
#include <vector>
#include <deque>
#include <list>
#include <utility> 
#include <iostream>
#include <math.h>

#include "Array.h"
#include "ref_count.h"

//: class to store edgels
typedef std::pair<double, double> point_2d;

class edgel:public ref_count
{
public:
    point_2d _pt; ///< the location of the edgel
    double _orientation;          ///< the orientation of the edgel
    double _strength;         ///< the strength of the edgel (typically gradient magnitude)
    int _id;                  ///< unique id
    
    edgel(double x, double y, double tan, double edge_strength, unsigned id):
    _pt(x,y), _orientation(tan), _strength(edge_strength), _id(id) {};
    
    // copy constructor
    edgel(edgel & e):_pt(e._pt), _orientation(e._orientation), _strength(e._strength),
    _id(e._id) {};
    
    ~edgel(){};
};

class edgemap
{
private:
    //: size of the edgemap
    unsigned _width, _height;
    
public:
    //: retinotopic map of edgels
    std::map<std::pair<int,int>, std::vector<edgel*> > _map;
    
    //: local list of edgels for easier traversal
    std::map<unsigned, edgel*> _list;
    
    //: constructor
    edgemap(unsigned width, unsigned height, arrayd &edgeinfo): _width(width), _height(height)
    {
        if ( edgeinfo.w() < 4 ){
            std::cerr << "Edge information size error\n";
        }
        for (unsigned i=0; i<edgeinfo.h(); i++)
        {
            insert(edgeinfo.val(i,0), edgeinfo.val(i,1), edgeinfo.val(i,2), edgeinfo.val(i,3), i);
        }
    }
    
    //: constructor
    edgemap(unsigned width, unsigned height): _width(width), _height(height){}
    
    //: destructor
    ~edgemap()
    {
        //go over each cell and delete the edgels
        for (std::map<unsigned,edgel*>::iterator it=_list.begin(); it!=_list.end(); ++it)
            it->second->unref();
        
        _map.clear();
        //also clear the list of edgels
        _list.clear();
    }
    
    //: put an edgel into the edgemap at the prescribed cell
    void insert(double x, double y, double tan, double edge_strength, unsigned id)
    {
        std::pair<int, int> key(round(x),round(y));
        edgel* edge = new edgel(x,y,tan,edge_strength,id);
        if (_map.find(key) == _map.end() ) {
            std::vector<edgel*> cell;
            cell.push_back(edge);
            _map.insert(std::make_pair(key,cell));
        }
        else {
            _map[key].push_back(edge);
        }
        _list.insert(std::make_pair(id,edge));
    }
    
    //Access functions
    unsigned width() const { return _width; }
    unsigned height() const { return _height; }
    unsigned num_edgels() const { return _list.size(); } ///< number of edgels in the edgemap
    
    //: print info to file
    void print(std::ostream& os)
    {
        std::map<unsigned, edgel*>::iterator mit;
        for (mit=_list.begin(); mit!=_list.end(); ++mit){
            os<<mit->second->_pt.first<<","<<mit->second->_pt.second<<","<<mit->second->_orientation
            <<","<<mit->second->_strength<<";"<<std::endl;
        }
    }
     
     /*void print_map()
     {
     for (std::map<std::pair<int,int>, std::vector<edgel*> >::iterator itm=_map.begin(); itm!=_map.end(); ++itm){
     std::cout<<"Cell ("<<itm->first.first<<","<<itm->first.second<<") Edge ID:";
     for(std::vector<edgel*>::iterator itv=itm->second.begin(); itv!=itm->second.end(); ++itv){
     std::cout<<(*itv)->_id<<"  ";
     }
     std::cout<<std::endl;
     }
     }*/
    
};


typedef std::deque<edgel*> edgel_chain; // differnt from lemsvxl
typedef std::list<edgel_chain*> edgel_chain_list;

#endif // edgemap_h
