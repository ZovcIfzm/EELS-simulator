#include "test_class.h"
#include <iostream>
using namespace std;

test_class::test_class(){}
test_class::test_class(double input){
	TCVariable = input;
	print();
}

void test_class::print(){
	cout << TCVariable << endl;
}