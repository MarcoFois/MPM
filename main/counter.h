#ifndef HAVE_COUNTER_H
#define HAVE_COUNTER_H

#include <iterator>

template <typename N>
struct range {

  struct
  iterator {
    friend class range;
    using difference_type = typename std::make_signed_t<N>;
    using iterator_category = std::input_iterator_tag;
    using value_type = N;
    using reference = N;
    using pointer = const N*;
    iterator &operator ++() { ++i_; return *this; }
    iterator operator ++(int) { iterator copy(*this); ++i_; return copy; }
    reference operator *() const { return i_; }
    bool operator ==(const iterator &other) const { return i_ == *other; }
    bool operator !=(const iterator &other) const { return i_ != *other; }
    /*
    bool operator <(const iterator &other) const { return i_ < *other; }
    difference_type operator -(const iterator &other) const { return static_cast<difference_type>(i_) - static_cast<difference_type>(*other); }
    difference_type operator -(const difference_type other) const { return static_cast<difference_type>(i_) - other; }
    difference_type operator +(const iterator &other) const { return static_cast<difference_type>(i_) + static_cast<difference_type>(*other); }
    difference_type operator +(difference_type other) const { return static_cast<difference_type>(i_) + other; }
    */
    // more operator implementations needed to make
    // this a random_access :
    // [],*,+,-,--,+=,-=,==,>=,<=,>
    protected: explicit iterator(difference_type start) : i_ (start) {};

  private:
    N i_;
  };
  iterator begin () const { return begin_; }
  iterator end () const { return end_; }
  range (N begin, N end) : begin_(begin), end_(end) {}
private:
  iterator begin_, end_;
};

#endif
