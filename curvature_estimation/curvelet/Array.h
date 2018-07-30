// This is Array.h
#ifndef Array_h_
#define Array_h_

/*****************************************************************************
// file: curveletmap.h
// brief: trivial array class encapsulating pointer Array
// author: Xiaoyan Li
// date: 01/26/2016
******************************************************************************/

#include <algorithm>
#include <iostream>
#include <assert.h>
#include "ref_count.h"

template <class T> class Array: public ref_count
{
public:
    Array() { _h=_w=0; _data=0; _free=false; }
    Array( const Array &A ){
        // note: don't use init here (_free not initialized)
        _h=A.h(); _w=A.w();
        _data=new T[_h*_w](); _free=true;
        
        for(unsigned r=0; r<_h; r++){
            for(unsigned c=0; c<_w; c++){
                _data[c*_h+r] = A.val(r,c);
            }
        }
    }
    virtual ~Array() { clear(); }
    void init(const unsigned &h, const unsigned &w) { clear(); _h=h; _w=w; _data=new T[h*w](); _free=true; }
    T& val(const unsigned &r, const unsigned &c) const { return _data[c*_h+r]; }
    unsigned h() const {return _h;}
    unsigned w() const {return _w;}
    void set_h(const unsigned &h){ _h=h; }
    void set_w(const unsigned &w){ _w=w; }
    void set_val(const unsigned &r, const unsigned &c, const T& value) { _data[c*_h+r]=value; }
    void fill(const T& value){
        for(unsigned r=0; r<_h; r++){
            for(unsigned c=0; c<_w; c++){
                _data[c*_h+r] = value;
            }
        }
    }
    void array_max(const Array &A1, const Array &A2){
        // must be reference otherwise _data will be delete when deleting the copy
        assert(A1.h()==A2.h());
        assert(A1.w()==A2.w());
        init(A1.h(), A1.w());
        for(unsigned r=0; r<_h; r++){
            for(unsigned c=0; c<_w; c++){
                _data[c*_h+r] = std::max(A1.val(r,c),A2.val(r,c));
            }
        }
    }
    void array_min(const Array &A1, const Array &A2){
        assert(A1.h()==A2.h());
        assert(A1.w()==A2.w());
        init(A1.h(), A1.w());
        for(unsigned r=0; r<_h; r++){
            for(unsigned c=0; c<_w; c++){
                _data[c*_h+r] = std::min(A1.val(r,c),A2.val(r,c));
            }
        }
    }
    
    void print(std::ostream& os) const {
        for(unsigned r=0; r<_h; r++){
            for(unsigned c=0; c<_w; c++){
                os<< val(r, c) << "  ";
            }
            os<<std::endl;
        }
    }
    T *_data;
    
private:
    void clear() { if(_free) delete [] _data; _h=_w=0; _data=0; _free=false; }
    
    unsigned _h, _w;
    bool _free;
};

// convenient typedefs
typedef Array<double> arrayd;
typedef Array<int> arrayi;

#endif
