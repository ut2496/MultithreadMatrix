#include "matrixnn.hpp"
 
#include <cmath>
#include <thread>
#include <vector>
#include <cstdio>
#include <cstring>
#include <random>
#include <chrono>
 
using namespace std::chrono;

int main(int argc, char *argv[]){
	std::vector<int> ranges;
	printf("%d\n",argc)	;
	if(argc == 3){
		ranges.push_back(atoi(argv[2]));
	} else {
		//ranges = {1,2,5,10,15,25,35};
		ranges = {3};
	}
	for(auto i: ranges){	
		
		MatrixNN a(i),b(i),c1(1),c2(1);
		a.FillRandom();
		b.FillRandom();
		
		a.print_mat();
		printf("\n");
		b.print_mat();
		printf("\n");
			
		double det = a.determ();
		//c2 = a.QuickTrans();
		//c2 = a.inverse();
		//c2 = a.QuickMult(b);
		//c2.print_mat();	
		printf("%f\n",det);		
	}
  return 1;
}
