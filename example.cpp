#include "GA.hpp"

int main()
{
	GA ga("bays29.txt", 100, 2000, 2.0, 0.5, 0.8, 0.05);
	ga.solve();

	return 0;
}