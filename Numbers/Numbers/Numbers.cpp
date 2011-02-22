#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
void main()
{
	ofstream myfile;
	myfile.open("numbers.txt");
	
	stringstream ss;
	//int count=0;

	for (int a=0; a <10; a++) {
		for (int b=0; b<10; b++) {
			for (int c=0; c<10; c++) {
				for (int d=0; d<10; d++) {
					for (int e=0; e<10; e++) {
						for (int f=0; f<10; f++) {
							for (int g=0; g<10; g++) {
								
								ss << "301" << a << b << c << d << e << f << g;
								myfile << ss.str() << "\n";
								ss.str("");
								//count++;
								//cout << 100*count/10000000 << endl;

							}
						}
					}
				}
			}
		}
	}


	myfile.close();
}