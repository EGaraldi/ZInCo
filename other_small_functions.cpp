/*============================================================

* File Name : other_small_functions.cpp

* Purpose : Here some small functions to perform various tasks:
			print_info, display_info, get_particle_index,
			build_file_name, my_fread, my_fwrite, my_exit,
			error_flag_check, spheres_intersect, DegreesToRadiant

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include<cstdarg>
#include"global_variables.h"
#include"functions_declaration.h"

//print_info is basically a fprintf + fflush to info.txt which
//is performed only by proc0
void print_info(const char * format, ...){
	//print info to the info.txt file
	va_list args;
	va_start(args, format);
	vfprintf(info, format, args);
	va_end(args);
	fflush(info);
}


//display_info is basically a printf + fflush to info.txt which
//is performed only by proc0
void display_info(const char * format, ...){
	//display info to screen
	va_list args;
	va_start(args, format);
	vprintf(format, args);
	va_end(args);
}


//get_particle_index given a position in the particles array, return the species (0-5)
//comparing to Ngas, Nhalo, ...
int get_particle_index(int position){
	if(position < Ngas) return 0;
	else if(position < Ngas+Nhalo) return 1;
	else if(position < Ngas+Nhalo+Ndisk) return 2;
	else if(position < Ngas+Nhalo+Ndisk+Nbulge) return 3;
	else if(position < Ngas+Nhalo+Ndisk+Nbulge+Nstars) return 4;
	else /*if(position < Ngas+Nhalo+Ndisk+Nbulge+Nstars+Nbndry)*/ return 5;
}


//build_file_name, incredibly, build the file name for ICs/SNAP
void build_file_name(char directory[], char name[], int file_number, int current_file, char output[]){
	//build the ICs/SNAP filename given some info
	char f[200];
	sprintf(f, "%s/%s", directory, name);
	if(file_number > 1) sprintf(output,"%s.%i",f,current_file); else sprintf(output,"%s",f);
}


//my_fread is fread + check everything is fine
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream){
	//fread+error check
	size_t nread;
	if((nread = fread(ptr, size, nmemb, stream)) != nmemb){
		display_info("ERROR: I/O error (fread) has occured: end of file\n");
		my_exit(1);
	}

	return nread;
}


//my_fread is fwrite + check everything is fine
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream){
	//fwrite+error check
	size_t nwrite;
	if((nwrite = fwrite(ptr, size, nmemb, stream)) != nmemb){
		display_info("ERROR: I/O error (fwrite) has occured\n");
		my_exit(1);
	}

	return nwrite;
}


void my_exit(int code){
	exit(code);
}


//check everything is fine
void error_flag_check(){
	if(error_flag != 0) my_exit(error_flag);
}


//spheres_intersect check for two given spheres if they intersect
bool spheres_intersect(Sphere one, Sphere two){
	//check for spheres intersection
	return ( pow(one.center[0]-two.center[0],2) + pow(one.center[1]-two.center[1],2) + pow(one.center[2]-two.center[2],2) <= pow(one.radius+two.radius, 2) );
}


//DegreesToRadiant, also incredibly, convert from degrees to radiant
float DegreesToRadiant(int degrees){
	const float kPi = 3.14159265359f;
	return degrees*kPi/180;
}
