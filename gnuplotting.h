#ifndef GNUPLOTTING.H
#define GNUPLOTTING.H

#include <string>
#include <cstdio>
using namespace std;

/*
Gnuplotting & operator << (Gnuplotting & outstream, char * command){
	outstream.cmd(command);
	return outstream;
}
*/

class Gnuplotting {

	public:
		Gnuplotting();
		Gnuplotting(string filename);
		Gnuplotting(char *filename);
		~Gnuplotting();
		void refreshPlot();
		void clearCanvas();
		void clear();
		void xrange(double &xmin, double &xmax);
		void yrange(double &xmin, double &xmax);
		void cmd(string & command);
		void cmd(char *  command);
		Gnuplotting & operator<<(char * command);




	private:
		string filename;
		FILE * pipe;




};










#endif

