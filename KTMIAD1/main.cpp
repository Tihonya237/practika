#include "grid.h"
#include "fem.h"

int main()
{
	setlocale(LC_ALL, "Rus");
	FEM fm;
	fm.Compute();
	return 0;
}