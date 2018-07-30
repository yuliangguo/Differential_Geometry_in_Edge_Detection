// This is ref_count.h
#ifndef ref_count_h
#define ref_count_h

/************************************************************************
 // file: ref_count.h
 // brief: reference-counting base class template
 // author: Xiaoyan Li
 // date: 01/24/2016
 ************************************************************************/

class ref_count
{
public:
    ref_count() : m_ref_count(0) { }
    
    // Copying an object should not copy the ref count.
    ref_count(ref_count const&) : m_ref_count(0) { }
    
    ref_count& operator=(ref_count const& /*rhs*/)
    { return *this; }
    
    virtual ~ref_count() {}
    
    void ref() { ++m_ref_count; }
    void unref() {
        if (--m_ref_count == 0)
            delete this; }
    long get_references() const { return m_ref_count; }
    bool is_referenced() const { return m_ref_count > 0; }
    
private:
    int m_ref_count;
};

#endif // ref_count_h
