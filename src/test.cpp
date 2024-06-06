#include <iostream>
#include "test.hpp"



template <>
int test<int> (int inp)
{ 
  std::cout << inp << endl;
  return inp;
};

int main() {

  std::cout << "test" << endl;
  int n = test<int>(4);
  cout << n << endl;
  return 0;

};

