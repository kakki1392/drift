#include "gnuplotting.h"
#include <cstdio>
#include <string>

using namespace std;

Gnuplotting::Gnuplotting(){
	pipe = popen("gnuplot -persist", "w");
}

Gnuplotting::~Gnuplotting(){
	pclose(pipe);
}

void Gnuplotting::cmd(char * command){
	fprintf(pipe,"%s\n",command);
	fflush(pipe);
}

void Gnuplotting::cmd(string & command){
	fprintf(pipe,"%s\n",command.c_str());
	fflush(pipe);
}

Gnuplotting & Gnuplotting::operator<<(char * command){
	cmd(command);
	return *this;
}
