// This is form_curvelet_process.h
#ifndef form_curvelet_process_h_
#define form_curvelet_process_h_

#include "curveletmap.h"
#include "edgemap.h"
#include "edgel_link_graph.h"
#include "curve_fragment_graph.h"
#include "polyarcmap.h"

class form_curvelet_process
{
public:
    form_curvelet_process(arrayd &edgeinfo, const unsigned &nrows,
                          const unsigned &ncols,const double &rad, const double &gap,
                          const double &dx, const double &dt, const double &token_len,
                          const double max_k,
                          const unsigned max_size_to_goup, bool bCentered_grouping=false,
                          bool bBidirectional_grouping=false);
    // default: ENO style curvelet.
    
    virtual ~form_curvelet_process();
    
    bool execute();
    
    void get_output_size(unsigned &h, unsigned &w, unsigned output_type);
    void get_output_arrary( arrayi &id_chain, arrayd &curvelet_info, unsigned output_type );
    
    // for debug
    edgemap* return_edgemap(){ return &_edgemap; }
    curveletmap* return_curveletmap(){ return &_curveletmap; }
    
private:
    //: form curvelets around each edgel in a greedy fashion
    void build_curvelets_greedy(bool use_flag, bool verbose);
    void build_curvelets_greedy_for_edge(edgel* eA, bool use_flag, bool forward, bool centered, bool leading);
    CC_curve_model_3d* form_a_hypothesis(edgel* ref_e, edgel* e2,  /*double dx, double dt,*/ double max_k,
                                         bool &ref_first, bool forward, bool centered, bool leading, bool verbose);
    void construct_the_link_graph(unsigned min_group_size, bool verbose);
    void form_links_from_a_curvelet(edgel* eA, curvelet* cvlet, unsigned min_group_size);
    void extract_regular_contours_from_the_link_graph(bool verbose);
    void extract_one_chains_from_the_link_graph();
    void clear_all_contours();
    void fit_polyarcs_to_all_edgel_chains( bool verbose );
    void fit_polyarc_to_chain(edgel_chain* chain);
    
protected:
    // input edge map
    edgemap _edgemap;
    // output curvelet map
    curveletmap _curveletmap;
    edgel_link_graph _linkgraph;
    curve_fragment_graph _fragraph;
    polyarcmap _polyarcmap;
    
    //various parameters
    double _rad; // radius of the grouping neighborhood around each edgel
    unsigned _nrad; // the range of edgel cells to search (this should contain all edgels within rad_)
    double _gap; // Distance between two consecutive edges in a link
    double _dx, _dt;
    double _max_k; // maximum curvature of curves in the curve bundle
    // WAITING FOR UPDATE!! should auto choose by program.
    double _token_len;  //Length of the edgel token (puts bounds on legal curvature)
    unsigned _max_size_to_group;
    bool _bCentered_grouping; // curvelets centered on the anchor edgel
    bool _bBidirectional_grouping; // curvelets in both direction
    
    //bool bFormCompleteCvletMap;
    
};

#endif
