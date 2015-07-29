/*============================================================

* File Name : global variables.cpp

* Purpose : initialization of global variables

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include"global_variables.h"

//path variables

char snap_dir[200];
char snap_name[200];
char Lic_dir[200];
char Lic_name[200];
char Hic_dir[200];
char Hic_name[200];
char output_dir[200];
int Lin_fnr;
int Hin_fnr;
int Sin_fnr;
int out_fnr;

//geometric variables
int cubes_per_side;
float c[3];

//struct used to read the header
struct io_header header1;
MPI_Datatype MPI_HEADER_t;

//struct that will contains the data of a single particle
struct particle_data;
MPI_Datatype MPI_PARTICLE_t;

//struct to describe a sphere
struct Sphere;
MPI_Datatype MPI_SPHERE_t;

//struct for resolution
struct resolution_info;

//arrays for the particles in the zoomed ICs 
vector< vector< particle_data > > new_particles;
int Ngas, Nhalo, Ndisk, Nbulge, Nstars, Nbndry;
LOIout lastID;

//geometric variables
float BoxSize;		//side of the simulation box (read from snapshot/ICs files)
float lambda;		//side of the cubic cell
float r_high_IC, r_medium_IC;	//radius of the new high- and mid-res region in the ICs
float c_IC[3];		//center of the new high- and mid-res region in the ICs
int mtot,ntot,ltot;	//number of cubes in each side of the BoxSize
resolution_info ***resolution_matrix;		//resolution matrix


//variables to read particles from a file
vector<particle_data> particles_in;	//array of particles properties
vector<LOIinHigh> idH;			//array of the ids
vector<LOIinLow> idL;                   //array of the ids


	
//other stuff
int particles_per_file;	//number of particles to write in each file
int error_flag;
float shift[3];		//displacement to put c at the center of the Box

// MPI variables
int world_size, my_rank;
int mpi_tag;


//resolution structure variables
int levels_number, species_number;	//number of resolution levels and species to process
int level_size[6];			//number of bubbles each level should contain
int level_cubic_region_side[6];		//side of the cubic region to merge in cells unit
float level_bubbles_radii[6];		//radius of the bubble of each level
char run_type[50];
bool cascade, zoom, dilution;		//variable to determine the type of the run
int cascade_random_seed;		//random seed for the cascade bubbles
int max_cascade_iter;			//max number of iteration for each bubble in the loop which places the bubbles
vector< vector< Sphere > > resolution_bubbles;	//vector which contains the bubbles

double init_time;
int rcounter;	//counter of mid- and low-res regions produced
int hcounter;	//counter of high-res cubes produced
FILE *info;	//file in which I'll write the program information

char fname[200];

#ifdef PRINT_IDS_CORRESPONDENCE
vector< vector< vector<LOIinHigh> > > old_ids;
FILE * id_corr;
#endif
