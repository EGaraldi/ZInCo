/*============================================================

* File Name : global_variables.h

* Purpose : contains definitions for the global variables

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#ifndef GLOBAL_VARIABLES_H
 #define GLOBAL_VARIABLES_H

 #ifndef LONGIDS_IN_HIGH
  #define LOIinHigh int
 #else
  #define LOIinHigh long
 #endif


 #ifndef LONGIDS_IN_LOW
  #define LOIinLow int
 #else
  #define LOIinLow long
 #endif


 #ifndef LONGIDS_OUT
  #define LOIout int
 #else
  #define LOIout long
 #endif

#endif



#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<cstring>

#include<vector>
using namespace std;
#include <algorithm>
using namespace std;

//#include<sys/time.h>
#include<ctime>



//path variables

extern char snap_dir[200];
extern char snap_name[200];
extern char Lic_dir[200];
extern char Lic_name[200];
extern char Hic_dir[200];
extern char Hic_name[200];
extern char output_dir[200];
extern int Lin_fnr;
extern int Hin_fnr;
extern int Sin_fnr;
extern int out_fnr;
extern int Lin_ftype;
extern int Hin_ftype;
extern int Sin_ftype;
extern int out_ftype;

//geometric variables
extern int cubes_per_side;
extern float c[3];	//coordinates of the center of the halo
extern const float PI;

//struct used to read the header
struct io_header{
	int npart[6];
	double massarr[6];
	double time;
	double redshift;
	int flag_sfr;
	int flag_feedback;
	int npartTotal[6];
	int flag_cooling;
	int num_files;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	int flag_age;
	int flag_metals;
	int npartTotalHW[6];
	int flag_entropy_ics;
	char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4 - 4 /*=60*/];	/* fills to 256 Bytes */

	io_header(){}

};
extern io_header header1;

//struct that will contains the data of a single particle
struct particle_data{
	float pos[3];
	float vel[3];
	float mass;
	float internal_energy;

	particle_data(){
		pos[0]=0.0f;
		pos[1]=0.0f;
		pos[2]=0.0f;
		vel[0]=0.0f;
		vel[1]=0.0f;
		vel[2]=0.0f;
		mass = 0.0f;
		internal_energy = -1.0f; //default to -1 in order to recognize when read
	}

	particle_data(float p0, float p1, float p2, float v0, float v1, float v2, float m, float ie){
		pos[0]=p0;
		pos[1]=p1;
		pos[2]=p2;
		vel[0]=v0;
		vel[1]=v1;
		vel[2]=v2;
		mass = m;
		internal_energy = ie;
	}

};


//struct to desibe a Sphere
struct Sphere{
	Sphere(){
		center[0]=0.0f;
		center[1]=0.0f;
		center[2]=0.0f;
		radius=0.0f;
	}
	Sphere(float (&c)[3], float r){
		center[0]=c[0];
		center[1]=c[1];
		center[2]=c[2];
		radius=r;
	}
	Sphere(float c0, float c1, float c2, float r){
		center[0]=c0;
		center[1]=c1;
		center[2]=c2;
		radius=r;
	}

	bool contains(float (&point)[3]){
		return ( pow(center[0]-point[0],2) + pow(center[1]-point[1],2) + pow(center[2]-point[2],2) <= radius*radius );
	}

	bool contains(float dist){
		return ( dist <= radius );
	}

	void setCenter(float (&c)[3]){
		center[0]=c[0];
		center[1]=c[1];
		center[2]=c[2];
	}

	float center[3];
	float radius;
};


//struct to store resolution info
struct resolution_info{
	int cubeID;
	int level;
	bool processed;
};

//arrays for the particles in the zoomed ICs 
extern vector< vector< particle_data > > new_particles;
extern int Ngas, Nhalo, Ndisk, Nbulge, Nstars, Nbndry;
extern LOIout lastID;

//geometric variables
extern float BoxSize;		//side of the simulation box (read from snapshot/ICs files)
extern float lambda;		//side of the cubic cell
extern float r_high_IC, r_medium_IC;	//radius of the new high- and mid-res region in the ICs
extern float c_IC[3];		//center of the new high- and mid-res region in the ICs
extern int mtot,ntot,ltot;	//number of cubes in each side of the BoxSize
extern resolution_info ***resolution_matrix;		//resolution matrix


//variables to read particles from a file
extern vector<particle_data> particles_in;	//array of particles properties
extern vector<LOIinHigh> idH;			//array of the ids
extern vector<LOIinLow> idL;                   //array of the ids


	
//other stuff
extern int particles_per_file;	//number of particles to write in each file
extern int error_flag;
extern float shift[3];		//displacement to put c at the center of the Box


// MPI variables
extern int world_size, my_rank;


//resolution structure variables
extern int levels_number, species_number;	//number of resolution levels and species to process
extern int level_size[6];                      //number of bubbles each level should contain
extern int level_cubic_region_side[6];         //side of the cubic region to merge in cells unit
extern float level_bubbles_radii[6];           //radius of the bubble of each level
extern char run_type[50];
extern bool cascade, zoom, dilution;            //variable to determine the type of the run
extern int cascade_random_seed;                //random seed for the cascade bubbles
extern int max_cascade_iter;                   //max number of iteration for each bubble in the loop which places the bubbles
extern vector< vector< Sphere > > resolution_bubbles;	//vector which contains the bubbles

extern double init_time;
extern int rcounter;	//counter of mid- and low-res regions produced
extern int hcounter;	//counter of high-res cubes produced
extern FILE *info;	//file in which I'll write the program information

extern char fname[200];

#ifdef PRINT_IDS_CORRESPONDENCE
extern vector< vector< vector<LOIinHigh> > > old_ids;
extern FILE * id_corr;
#endif
