#pragma once
#include "main_sequence.h"
//TODO
//Test limits of roundoff error for variable checking
//
//


//Process code moved to main_sequence.cpp rather than keeping it in source, 
//since other files (data_processing.cpp) must include it. (Can't include source, it copies the main function)
int main(){
	mainSequence();
	return 0;
}