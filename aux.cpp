#include <iostream>
#include <TNL/Containers/Array.h>
#include <TNL/Containers/Vector.h>
#include <TNL/Meshes/Mesh.h>

using namespace TNL;
using namespace TNL::Containers;

const int size = 10;

void f1()
{
	Array< float, Devices::Host > a( size );
	for(int i = 0; i < size; i++)
	{
		a[i] = size - i;
	}
	std::cout << a << std::endl;
}

void f2()
{
	Array< int /*Devices::Host by default*/> b( size );
	b = 69;
	std::cout << b << std::endl;
}

void f3()
{
	Vector< double, Devices::Host> u = { 1, 2, 3 };
	std::cout << u << std::endl;
}

void f4()
{
	Vector< int, Devices::Host > v( size );
	for( int i = 0; i < size; i++ )
	{
		v[i] = i;
	}
	std::cout << v << std::endl;
}

int main()
{
	f1();
	f2();

	f3();
	f4();

	return 0;
}
