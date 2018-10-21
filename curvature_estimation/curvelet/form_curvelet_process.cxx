// This is form_curvelet_process.cxx

#include <iostream>
// #include <cstring>
#include <ctime>
#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>
#include "form_curvelet_process.h"
#include "knot.h"

//extern double To2Pi (double angle);
extern double sq_distance(const point_2d& A, const point_2d& B);
extern double vPointPoint (point_2d& start, point_2d& end);
extern double dot (double v1, double v2);

//: data structure to hold pairwise hypothesis around an edgel
struct sel_hyp {
    bool _ref_first;           ///< the direction of groouping
    double _d;                 ///< distance between them
    bool _flag;                ///< if this hypothesis has been considered before
    edgel* _eN;          ///< neighboring edgel
    CC_curve_model_3d* _cm;    ///< the curve bundle of the appropriate model
};

inline bool comp_dist_hyps_less(const sel_hyp &h1, const sel_hyp &h2)
{
    if (h1._d < h2._d)
        return true;
    else
        return false;
}

//: Constructor
form_curvelet_process::form_curvelet_process(arrayd &edgeinfo, const unsigned &nrows,
                                             const unsigned &ncols,const double &rad, const double &gap,
                                             const double &dx, const double &dt, const double &token_len,
                                             const double max_k,
                                             const unsigned max_size_to_group, bool bCentered_grouping,
                                             bool bBidirectional_grouping):
_edgemap(ncols, nrows, edgeinfo), _curveletmap(&_edgemap), _polyarcmap(&_edgemap), _rad(rad), _gap(gap),
_dx(dx), _dt(dt), _max_k(max_k),  _token_len(token_len),
_max_size_to_group(max_size_to_group), _bCentered_grouping(bCentered_grouping),
_bBidirectional_grouping(bBidirectional_grouping)
{
    _nrad = ceil(rad);
}

//: Destructor
form_curvelet_process::~form_curvelet_process()
{
}

//: Execute the process
bool form_curvelet_process::execute()
{
    //start the timer
    std::clock_t t0;
    t0 = clock();
    
    bool verbose = false;
    
    // build curvelets
    build_curvelets_greedy(false, verbose);
    std::clock_t t = clock();
    double group_time = ((float)(t-t0))/CLOCKS_PER_SEC;
    t0 = t;
    std::cout << "Time taken to form groups and construct curvelet map: " << group_time << " sec" << std::endl;
    
    // construct edgel link graph
    construct_the_link_graph(4, verbose); //default: min_size_to_link = 4
    t = clock();
    group_time = ((float)(t-t0))/CLOCKS_PER_SEC;
    t0 = t;
    std::cout << "Time taken to construct edgel link graph: " << group_time << " sec" << std::endl;

    // construct curve fragment graph
    extract_regular_contours_from_the_link_graph(verbose);
    t = clock();
    group_time = ((float)(t-t0))/CLOCKS_PER_SEC;
    t0 = t;
    std::cout << "Time taken to construct curve fragment graph: " << group_time << " sec" << std::endl;

    // fit polyarc
    fit_polyarcs_to_all_edgel_chains(verbose);
    t = clock();
    group_time = ((float)(t-t0))/CLOCKS_PER_SEC;
    t0 = t;
    std::cout << "Time taken to fit polyarcs: " << group_time << " sec" << std::endl;
    
//    _polyarcmap.resize(2.0);
    return true;
}


/**********************************************************************************
 // functions to form curvelets
**********************************************************************************/

//: form curvelets around each edgel in a greedy fashion
void form_curvelet_process::build_curvelets_greedy(bool use_flag, bool verbose)
{
    if (verbose){
        std::cout << "Building All Possible Curvelets (Greedy) ..." << std::flush;
    }
    
    // Curvelets can be formed in 4 ways now:
    // 1. regular anchor centered
    // 2. anchor centered bidirectional
    // 3. Anchor leading bidirectional
    // 4. ENO style (anchor leading or trailing but in the same direction)
    
    if (_bCentered_grouping) {
        if (_bBidirectional_grouping) {
            for (unsigned i=0; i<_edgemap._list.size(); i++) {
                edgel* eA = _edgemap._list[i];
                // centered_ && bidir_
                build_curvelets_greedy_for_edge(eA, use_flag, true, _bCentered_grouping,false);//first in the forward direction
                build_curvelets_greedy_for_edge(eA, use_flag, false,_bCentered_grouping, false); //then in the other direction
            }
        } else {
            for (unsigned i=0; i<_edgemap._list.size(); i++) {
                edgel* eA = _edgemap._list[i];
                // centered_ && !bidir_
                build_curvelets_greedy_for_edge(eA, use_flag, true, _bCentered_grouping, false); //in the forward direction
            }
        }
    } else {
        if (_bBidirectional_grouping) {
            for (unsigned i=0; i<_edgemap._list.size(); i++) {
                edgel* eA = _edgemap._list[i];
                build_curvelets_greedy_for_edge(eA, use_flag, true, _bCentered_grouping, true); //forward half
                build_curvelets_greedy_for_edge(eA, use_flag, false, _bCentered_grouping, true); //backward half
                // !centered_ && bidir_
            }
        } else {
            for (unsigned i=0; i<_edgemap._list.size(); i++) {
                edgel* eA = _edgemap._list[i];
                build_curvelets_greedy_for_edge(eA,use_flag, true, _bCentered_grouping, true); //forward half
                build_curvelets_greedy_for_edge(eA, use_flag, true, _bCentered_grouping, false); //ENO style forward
                // !centered_ && !bidir_
            }
        }
    }

    if (verbose)
    {
        std::cout << "done!" << std::endl;
        _curveletmap.print(std::cout);
    }
    
}

//: form curvelets around the given edgel in a greedy fashion
void form_curvelet_process::build_curvelets_greedy_for_edge(edgel* eA, bool use_flag, bool forward, bool centered, bool leading)
{
    // 1) construct a structure to temporarily hold the pairwise-hypotheses
    std::vector<sel_hyp> eA_hyps;
    
    //get the grid coordinates of this edgel
    unsigned const ii = round(eA->_pt.first);
    unsigned const jj = round(eA->_pt.second);
    
    unsigned const rad_sqr = _rad*_rad;
    
    bool verbose = false;
    // if (eA->_id==3) verbose = true;
    // 2) iterate over the neighboring cells around this edgel
    for (int xx=(int)ii-(int)_nrad; xx<=(int)(ii+_nrad) ; xx++){
        for (int yy=(int)jj-(int)_nrad; yy<=(int)(jj+_nrad) ; yy++){
            
            if (xx<0 || xx>=_edgemap.width() || yy<0 || yy>=_edgemap.height() )
                continue;
            
            //for all the edgels in its neighborhood
            for (unsigned k=0; k<_edgemap._map[std::make_pair(xx, yy)].size(); k++){
                edgel* eB = _edgemap._map[std::make_pair(xx,yy)][k];
                if (eB->_id==eA->_id) continue;
                
                // 3) do a better check of circular radius because the bucketing neighborhood is very coarse
                if ( sq_distance(eA->_pt, eB->_pt) > rad_sqr) continue;
                
                // 4) form pair-wise hypotheses
                bool ref_first;
                if (verbose) std::cout << '\n' << eB->_id << std::endl;
                CC_curve_model_3d* cm = form_a_hypothesis(eA, eB, /*_dx, _dt,*/ _max_k, ref_first, forward, centered, leading, verbose);
                
                if (cm){
                    if (cm->bundle_is_valid()){ //if legal, record the hypothesis
                        sel_hyp cur_hyp;
                        
                        cur_hyp._eN = eB;
                        cur_hyp._cm = cm;
                        cur_hyp._ref_first = ref_first;
                        cur_hyp._d = sq_distance(eA->_pt, eB->_pt);
                        cur_hyp._flag = false;
                        
                        eA_hyps.push_back(cur_hyp);
                    }
                    else
                        delete cm;
                }
            }
        }
    }
    
    // 5) first sort the pair-wise hyps by distance
    std::sort(eA_hyps.begin(), eA_hyps.end(), comp_dist_hyps_less);
    
    // 6) for each pair-wise hyps formed by this edgel
    for (unsigned h1=0; h1<eA_hyps.size(); h1++)
    {
        //this edgel has been used up already (only if we are using flags)
        if (use_flag && eA_hyps[h1]._flag) continue;
        
        // 7) initialize a new edgel chain that will grow in a greedy depth-first fashion
        edgel_chain cur_edgel_chain; //chain can grow either way
        
        //insert ref edgel first
        cur_edgel_chain.push_back(eA);
        
        // 8)initialize the shrinking curve bundle from the cur hypothesis
        CC_curve_model_3d* cur_cm = eA_hyps[h1]._cm;
        
        if (verbose){
            std::cout << eA_hyps[h1]._eN->_id << std::endl;
            cur_cm->print_bundle(std::cout);
        }
        
        // 9) attempt to integrate the other pair-wise hypotheses (DEPTH-FIRST SEARCH)
        for (unsigned h2=0; h2<eA_hyps.size(); h2++)
        {
            if (h2==h1){
                //add the second edgel of the reference pair (no need for any computation)
                if (eA_hyps[h1]._ref_first) cur_edgel_chain.push_back(eA_hyps[h1]._eN);
                else                        cur_edgel_chain.push_front(eA_hyps[h1]._eN);
                
                // check the size of the grouping
                if (cur_edgel_chain.size() >= _max_size_to_group)
                    break;
                else
                    continue;
            }
            
            // 10) for the others we need to check for consistency with the current hypothesis
            //     (i.e., compute intersection of curve bundles)
            CC_curve_model_3d* new_cm = cur_cm->intersect(eA_hyps[h2]._cm);
            
            // 11) if the intersection is valid, we can add this edgel to the grouping
            if (new_cm->bundle_is_valid())
            {
                //reassign the curve bundle for the growing grouping
                if (verbose){
                    std::cout << eA_hyps[h2]._eN->_id << std::endl;
                    new_cm->print_bundle(std::cout);
                }
                
                if (cur_cm != eA_hyps[h1]._cm)
                    delete cur_cm;
                cur_cm = new_cm;
                
                // 12) add the new edgel to the growing edgel chain
                if (eA_hyps[h2]._ref_first) cur_edgel_chain.push_back(eA_hyps[h2]._eN);
                else                        cur_edgel_chain.push_front(eA_hyps[h2]._eN);
                
                //flag this hypothesis
                eA_hyps[h2]._flag = true;
            }
            else
                delete new_cm; //delete this cb because it is not needed
            
            // 13) check the size of the grouping
            if (cur_edgel_chain.size() >= _max_size_to_group)
                break;
            
        }
        
        // 14) form a new curvelet and assign the curve bundle and the edgel chain to it...
        //     ...if it passes a few tests :
        //                             (a) if it doesn't already exist and
        //                             (b) its a reasonable fit and
        //                             (c) its relatively symmetric
        if (cur_edgel_chain.size()>2 &&
            !_curveletmap.does_curvelet_exist(eA, cur_edgel_chain) &&
            cur_cm->curve_fit_is_reasonable(cur_edgel_chain, eA, _dx))
        {
            curvelet* new_cvlet = new curvelet(eA, cur_cm, cur_edgel_chain, forward);
            
            //compute curvelet quality
            new_cvlet->compute_properties(_rad, _token_len);
            
            //add links from the ref edgel to this curvelet
            _curveletmap.add_curvelet(new_cvlet);
        }
        else {
            //delete cur_cm if no curvelet is formed
            if (cur_cm != eA_hyps[h1]._cm)
                delete cur_cm;
        }
    }
    
    // 15) now delete all the hyps for this edgel, before moving on to a new edgel
    for (unsigned h1=0; h1<eA_hyps.size(); h1++)
        delete eA_hyps[h1]._cm;
    
    eA_hyps.clear();
}

//: form a curve hypothesis of the appropriate model given a pair of edgels
CC_curve_model_3d* form_curvelet_process::form_a_hypothesis(edgel* ref_e, edgel* e2, /*double dx, double dt,*/ double max_k, bool &ref_first,  bool forward, bool centered, bool leading, bool verbose)
{
    // Do the dot product test to determine if this pair of edgels can produce a valid hypothesis
    // The two possible hypotheses are:
    //    a) (ref_e->e2)
    //    b) (e2->ref_e)
    
    double dx = _dx, dt = _dt;
    const double ref_dir = vPointPoint(ref_e->_pt, e2->_pt); //reference dir (ref_e->e2)
    //ref centered grouping
    if (centered) {
        if (dot(ref_dir, ref_e->_orientation)>0) {
            if (forward){
                ref_first = true;
            }
            else {
                ref_first = false;
            }
        }
        else {
            if (forward){
                ref_first = false;
            }
            else {
                ref_first = true;
            }
        }
        // the ref point provide adaptive uncertainty which is not used
        return new CC_curve_model_3d(ref_e, e2, dx, dt, max_k, verbose);
    }
    else { //not centered
        if (dot(ref_dir, ref_e->_orientation)>0 && forward && leading) {
            ref_first = true;
            return new CC_curve_model_3d(ref_e, e2, dx, dt, max_k, verbose);
        }
        if (dot(ref_dir, ref_e->_orientation)<0) {
            if (forward && !leading){
                ref_first = false;
                return new CC_curve_model_3d(ref_e, e2, dx, dt, max_k, verbose);
            }
            if (!forward){
                ref_first = true;
                return new CC_curve_model_3d(ref_e, e2, dx, dt, max_k, verbose);
            }
        }
    }
    return 0;
}


/**********************************************************************************
// functions to construct edgel link graph
**********************************************************************************/

//: form the link graph from the existing edgel groupings
// method = 0 : include all the curvelets
void form_curvelet_process::construct_the_link_graph(unsigned min_group_size, bool verbose)
{
    //1) clear the link graph and form a new one
    _linkgraph.clear();
    _linkgraph.resize(_edgemap.num_edgels());
    
    /*if (!use_anchored_curvelets_){
     // First update the curvelet map to include all links
     _curveletmap.form_full_cvlet_map();
     }*/
    
    // 2) now construct the link graph from the curvelet map
    if (verbose) {
        std::cout << "Constructing the Link Graph using (N >= " << min_group_size << ")..." ;
    }
    
    // 2a) go over all the curvelets and reset the used flags
    for (unsigned i=0; i<_edgemap._list.size(); i++)
    {
        //for all curvelets that are larger than the group size threshold
        cvlet_list::iterator cv_it = _curveletmap.getcurvelets(i).begin();
        for ( ; cv_it!=_curveletmap.getcurvelets(i).end(); cv_it++)
            (*cv_it)->_used = false; //reset this flag
    }
    
    // 2b) go over all the curvelets above the min size and determine which links can be formed
    for (unsigned i=0; i<_edgemap._list.size(); i++)
    {
        edgel* eA = _edgemap._list[i];
        
        //for all curvelets that are larger than the group size threshold
        cvlet_list::iterator cv_it = _curveletmap.getcurvelets(i).begin();
        for ( ; cv_it!=_curveletmap.getcurvelets(i).end(); cv_it++){
            curvelet* cvlet = (*cv_it);
            if (cvlet->order() < min_group_size) continue;
            
            //form all possible links from this curvelet
            form_links_from_a_curvelet(eA, cvlet, min_group_size);
        }
    }
    if (verbose) {
        std::cout << "done!" << std::endl;
        _linkgraph.print(std::cout);
    }
    
    //4) after forming the link graph, determine the degree of overlap of the links in the graph
    for (unsigned i=0; i<_linkgraph._cLinks.size(); i++){
        //all child links of each edgel covers all the links
        link_list::iterator l_it = _linkgraph._cLinks[i].begin();
        for (; l_it != _linkgraph._cLinks[i].end(); l_it++){
            (*l_it)->set_degree_overlap();
        }
    }
    
}

//: form all appropriate links from a curvelet
void form_curvelet_process::form_links_from_a_curvelet(edgel* eA, curvelet* cvlet, unsigned min_group_size)
{
    //Link all edgels in the group
    //
    // Explanation:
    // if xyzAbcd then X-Y, Y-Z, Z-A, A-B, B-C and C-D
    
    for (unsigned k=0; k<cvlet->_edgel_chain.size()-1; k++)
        _linkgraph.link_edgels(cvlet->_edgel_chain[k], cvlet->_edgel_chain[k+1], cvlet);
    
    ///all cvlets are used
    cvlet->_used = true;
}


/**********************************************************************************
// functions to construct curve fragment graph
**********************************************************************************/

//: extract the regular contours from the link graph
void form_curvelet_process::extract_regular_contours_from_the_link_graph(bool verbose)
{
    if (verbose) {
        std::cout << "Constructing the Curve Fragment Graph ..." ;
    }
    
    //first remove any existing contours
    clear_all_contours();
    
    //clear the linked flags
    _linkgraph.clear_linked_flag();
    
    //exract one chains from the main link graph
    extract_one_chains_from_the_link_graph();
    
    if (verbose) {
        std::cout << "done!" << std::endl;
        _fragraph.print(std::cout);
    }
}

//: clear all contours
void form_curvelet_process::clear_all_contours()
{
    //reset linked flags on the link graph
    _linkgraph._linked.assign(_edgemap.num_edgels(), false);
    
    //clear the curve fragment graph
    _fragraph.clear();
    
    //form a new contour fragment graph
    _fragraph.resize(_edgemap.num_edgels());
}

//: extract the one chains from a given link graph (from the primary link grpah of ELG)
void form_curvelet_process::extract_one_chains_from_the_link_graph()
{
    //initialize the curve fragment map
    _fragraph.resize(_edgemap.num_edgels());
    double thres=_gap*_gap;
    //now look for edgel chains
    //Rules:
    //    (a) start from an edgel that is locally legal
    //    (b) trace in both directions until an illegal edgel is reached
    //    (c) prune out short chains
    
    for (unsigned i=0; i<_edgemap.num_edgels(); i++)
    {
        edgel* first_edgel = _edgemap._list[i];
        
        //if it's already linked, ignore
        if (_linkgraph._linked[first_edgel->_id])
            continue;
        
        // Check edgel to see if it is legal to start a chain here
        if (_linkgraph.edgel_is_legal_first_edgel(first_edgel))
        {
            //start a chain from this edgel
            edgel_chain* chain = new edgel_chain();
            
            //add the first edgel to the chain
            chain->push_back(first_edgel);
            _linkgraph._linked[first_edgel->_id] = true; //mark it as linked
            
            //now start tracing FORWARD from its child
            edgel* eA = _linkgraph._cLinks[first_edgel->_id].front()->_ce;
            edgel* eB;
            if(sq_distance(first_edgel->_pt,eA->_pt)<thres)
                chain->push_back(eA);
            
            //trace FORWARD through the link graph until an illegal or terminal edgel is reached
            while (_linkgraph.edgel_is_legal(eA))
            {
                // Mark the last edgel as linked.
                
                // Note:
                //   By doing this here, we can get the edgels at junctions to be added to the contour
                //   without marking them as linked. This means that other contours arriving at
                //   the junction can also claim the junction edgel as being on their chains.
                _linkgraph._linked[eA->_id] = true;
                
                //is this a terminal edgel?
                if (_linkgraph._cLinks[eA->_id].size()==0)
                    break; //terminate chain
                
                //else advance to child node
                eB=eA;
                eA = _linkgraph._cLinks[eA->_id].front()->_ce;
                if(sq_distance(eB->_pt,eA->_pt)<thres)
                    chain->push_back(eA);
                else break;
            }
            
            //Note: Junction edgels will still be marked as unlinked after the tracing is done!
            
            //now start tracing BACKWARD from the first edgel
            
            //first determine if this is a closed contour
            //with closed contours, the chain might already include the first edgel twice
            if (eA != first_edgel){
                //not a closed contour, start tracing
                eA = _linkgraph._pLinks[first_edgel->_id].front()->_pe;
                if(sq_distance(first_edgel->_pt,eA->_pt)<thres)
                    chain->push_front(eA);
            }
            
            while (_linkgraph.edgel_is_legal(eA))
            {
                // Mark the last edgel as linked.
                _linkgraph._linked[eA->_id] = true;
                
                //is this a terminal edge?
                if (_linkgraph._pLinks[eA->_id].size()==0)
                    break; //terminate chain
                
                //else advance to parent node
                eB=eA;
                eA = _linkgraph._pLinks[eA->_id].front()->_pe;
                if(sq_distance(eB->_pt,eA->_pt)<thres)
                    chain->push_front(eA);
                else break;      
            }
            
            //save the current chain on the curve fragment graph
            if (chain->size()>2)
                _fragraph.insert_fragment(chain); //prune out the short ones
            else
                delete chain;
            
        }
    }
}


/**********************************************************************************
// functions to fit polyarcs
**********************************************************************************/

//: fit C^1 polyarcs to all the unambiguous edgel chains in the image
void form_curvelet_process::fit_polyarcs_to_all_edgel_chains( bool verbose )
{
    if (verbose) {
        std::cout << "Fit C^1 polyarcs to all the unambiguous edgel chains ..." ;
    }

    //go over each edgel chain in the link grpah and try to fit polyarcs to them
    edgel_chain_list::iterator f_it = _fragraph._frags.begin();
    for (; f_it != _fragraph._frags.end(); f_it++)
    {
        edgel_chain* chain = (*f_it);
        fit_polyarc_to_chain(chain);
    }
    
    if (verbose)
    {
        std::cout << "done!" << std::endl;
        _polyarcmap.print(std::cout);
    }
}

double compute_BC_saliency(CC_curve_model_3d* cm1, CC_curve_model_3d* cm2)
{
    double vol = 0.0;

    for (unsigned i=0; i<cm1->_Kmax.h(); i++)
        for (unsigned j=0; j<cm1->_Kmax.w(); j++)
            if (cm1->_Kmax.val(i,j)>cm1->_Kmin.val(i,j) && cm2->_Kmax.val(i,j)>cm2->_Kmin.val(i,j))
                vol += 1.0;

    return vol;
}

//: attempt to fit polyarcs to edgel chains (knot-based algorithm)
void form_curvelet_process::fit_polyarc_to_chain(edgel_chain* chain)
{
    unsigned N = chain->size();
    unsigned i = 0, j = 0;
    
    // 1) Grow the edge chain until C2 condition is not satisfied
    edgel *eA, *eB;
    
    while (i<N-1) {
        
        eA = (*chain)[i];
        eB = (*chain)[i+1];
        
        i++;
        bool ref_first; //redundant value, since we already know the order
        //construct the base curvelet then grow the curvelet along the chain
        CC_curve_model_3d* cm = form_a_hypothesis(eA, eB, /*_dx, _dt,*/ _max_k, ref_first, true, false, true, false);
        // holding temporary curve model to be grown
        edgel_chain cur_edgel_chain;
        
        if (cm && cm->bundle_is_valid()) {
            
            j = i;
            //insert base edgels
            cur_edgel_chain.push_back(eA);
            cur_edgel_chain.push_back(eB);
            
            while (j<N-1) {
                eB = (*chain)[j+1];
                
                // attempt to integrate the other edgel along the edgel chain
                CC_curve_model_3d* cm2 = form_a_hypothesis(eA, eB, /*_dx, _dt,*/ _max_k, ref_first, true, false, true, false);

                if (cm2 && cm2->bundle_is_valid()) {
                    
                    // only curve models with the same ref edgel can be intersect
                    CC_curve_model_3d* new_cm = cm->intersect(cm2);
                    
                    if (new_cm->bundle_is_valid()) {
                        //reassign the curve bundle for the growing grouping
                        delete cm;
                        delete cm2;
                        cm = new_cm;
                        
                        // add the new edgel to the growing edgel chain
                        cur_edgel_chain.push_back(eB);
                    }
                    else
                    {
                        delete new_cm;
                        delete cm2;
                        break;
                    }
                }
                else
                {
                    if (cm2) {
                        delete cm2;
                    }
                    break;
                }
                j++;
            }
            
            // form a new curvelet and assign the curve bundle and the edgel chain to it...
            //     ...if it passes a few tests :
            //                             (a) if it doesn't already exist and
            //                             (b) its a reasonable fit and
            //                             (c) its relatively symmetric
            if (cur_edgel_chain.size()>1 &&
                cm->curve_fit_is_reasonable(cur_edgel_chain, eA, _dx))
            {
                curvelet* new_cvlet = new curvelet(eA, cm, cur_edgel_chain, true);
                
                //compute curvelet quality
                new_cvlet->compute_properties(_rad, _token_len);
                
                //add links from the ref edgel to this curvelet
                _polyarcmap.add_curvelet(new_cvlet);
                
                // move to the tail edgel
                i = j;
            }
            else {
                //delete cur_cm if no curvelet is formed
                delete cm;
            }
        }
    }
    
    _polyarcmap.compute_best_fit_according_to_curvelets();
}


void  form_curvelet_process::get_output_size(unsigned &h, unsigned &w, unsigned output_type)
{
    switch (output_type) {
        case 0:
            // out put the curvelet map
            h = _curveletmap.num_curvelets();
            w = _max_size_to_group+1;
            break;
            
        case 1:
            // out put the curve fragment graph
            h = _fragraph.num_frags();
            w = _fragraph.edgel_chain_max_length();
            break;
        case 2:
            // out put the poly arc map
            h = _polyarcmap.num_curvelets();
            w = _polyarcmap.edgel_chain_max_length();
            break;
            
        default:
            std::cerr << "Error: invalid output type." << std::endl;
            break;
    }
}

void form_curvelet_process::get_output_arrary( arrayi &id_chain, arrayd &curvelet_info, unsigned output_type )
{
    unsigned h,w;
    get_output_size(h,w,output_type);
    
    switch (output_type) {
        case 0:
            // out put the curvelet map
            if (h==id_chain.h() && w==id_chain.w() && h==curvelet_info.h() && curvelet_info.w()==10) {
                _curveletmap.get_output_array(id_chain,curvelet_info);
            }
            else
                std::cerr<<"Error: incorrect output array."<<std::endl;
            break;
            
        case 1:
            // out put the curvelet map
            if (h==id_chain.h() && w==id_chain.w()) {
                _fragraph.get_output_array(id_chain,curvelet_info);
            }
            else
                std::cerr<<"Error: incorrect output array."<<std::endl;
            break;
            
        case 2:
            // out put the poly arc map
            if (h==id_chain.h() && w==id_chain.w() && h==curvelet_info.h() && curvelet_info.w()==12) {
                _polyarcmap.get_output_array(id_chain,curvelet_info);
            }
            else
                std::cerr<<"Error: incorrect output array."<<std::endl;
            break;
            
        default:
            break;
    }
   
}

