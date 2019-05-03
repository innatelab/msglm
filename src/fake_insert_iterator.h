#pragma once

#include <iterator>

template<class Container>
class fake_insert_iterator:
    public std::iterator<std::output_iterator_tag,void,void,void,void>
{
public:
  typedef Container container_type;
  explicit fake_insert_iterator() {}
  fake_insert_iterator<Container>& operator= (const typename Container::value_type& value)
    { return *this; }
  //back_insert_iterator<Container>& operator= (typename Container::value_type&& value)
  //{ return *this; }
  fake_insert_iterator<Container>& operator* ()
    { return *this; }
  fake_insert_iterator<Container>& operator++ ()
    { return *this; }
  fake_insert_iterator<Container> operator++ (int)
    { return *this; }
};