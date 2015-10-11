/*============================================================

* File Name : initialize_check_finalize_variables.cpp

* Purpose : The following functions deal with variables during the run,
			mainly initializing (MPI-related stuff, timers and arguments of the call),
			checking the user inputs are consistent, checking for ignored particles (i.e.
			particles outside ALL the resolution regions defined) and finalizing MPI.

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include"global_variables.h"
#include"functions_declaration.h"

//The following function is called at the beginning of the run and deals (mainly) with:
// -initializing MPI
// -building MPI custom datatypes
// -checking the call sequence
void initialize(int argc){

	//MPI-related variables
	mpi_tag=0;
	error_flag = 0;
	lastID=0;

	//Initialize MPI
	MPI_Init(NULL, NULL);

	// Get the number of processes and the current process' rank
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	//starting the clock
	init_time = MPI_Wtime();

	//check the call sequence
	if(argc<2){ display_info("Missing parameters file, please run using: <exec> <parameters file name>\n"); error_flag=1; }
	error_flag_check();

	//create a MPI_type for struct particle_data 
	{const int nitems=4;
	int blocklengths[4] = {3,3,1,1};
	MPI_Datatype types[4] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
	MPI_Aint offsets[4] = {
				(MPI_Aint)offsetof(particle_data, pos), 
				(MPI_Aint)offsetof(particle_data, vel), 
				(MPI_Aint)offsetof(particle_data, mass),
				(MPI_Aint)offsetof(particle_data, internal_energy)
				};

	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_PARTICLE_t);
	MPI_Type_commit(&MPI_PARTICLE_t);
	}

	// create a MPI_type for struct io_header
	{const int nitems=18;
	int blocklengths[18] = {6, 6, 1, 1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 60};
	MPI_Datatype types[18] = { 
				  MPI_INT, 
				  MPI_DOUBLE,  
				  MPI_DOUBLE,  
				  MPI_DOUBLE,  
				  MPI_INT,  
				  MPI_INT,  
				  MPI_INT,  
				  MPI_INT,  
				  MPI_INT,  
				  MPI_DOUBLE,  
				  MPI_DOUBLE,  
				  MPI_DOUBLE,  
				  MPI_DOUBLE,  
				  MPI_INT,  
				  MPI_INT,  
				  MPI_INT,  
				  MPI_INT,  
				  MPI_CHAR 
				  };
	MPI_Aint offsets[18] = {
				(MPI_Aint)offsetof(io_header, npart), 
				(MPI_Aint)offsetof(io_header, massarr), 
				(MPI_Aint)offsetof(io_header, time), 
				(MPI_Aint)offsetof(io_header, redshift), 
				(MPI_Aint)offsetof(io_header, flag_sfr), 
				(MPI_Aint)offsetof(io_header, flag_feedback), 
				(MPI_Aint)offsetof(io_header, npartTotal), 
				(MPI_Aint)offsetof(io_header, flag_cooling), 
				(MPI_Aint)offsetof(io_header, num_files), 
				(MPI_Aint)offsetof(io_header, BoxSize), 
				(MPI_Aint)offsetof(io_header, Omega0), 
				(MPI_Aint)offsetof(io_header, OmegaLambda), 
				(MPI_Aint)offsetof(io_header, HubbleParam), 
				(MPI_Aint)offsetof(io_header, flag_age), 
				(MPI_Aint)offsetof(io_header, flag_metals), 
				(MPI_Aint)offsetof(io_header, npartTotalHW), 
				(MPI_Aint)offsetof(io_header, flag_entropy_ics), 
				(MPI_Aint)offsetof(io_header, fill)
				};

	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_HEADER_t);
	MPI_Type_commit(&MPI_HEADER_t);
	}

	// create a MPI_type for struct Sphere
	{const int nitems=2;
	int blocklengths[2] = {3,1};
	MPI_Datatype types[2] = {MPI_FLOAT, MPI_FLOAT};
	MPI_Aint offsets[2] = {(MPI_Aint)offsetof(Sphere, center), (MPI_Aint)offsetof(Sphere, radius)};
	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_SPHERE_t);
	MPI_Type_commit(&MPI_SPHERE_t);
	}

	//the vector contaiing the new particles created is a 6*? vector, since GADGET allows (max) 6 particles types
	new_particles.resize(6);

	#ifdef PRINT_IDS_CORRESPONDENCE
	old_ids.resize(6);
	#endif
}


//The following function check the internal consistency of the user inputs
void check_variables(){
	int i;

	//check the new ICs require maximum 6 species
	if(levels_number * species_number > 6){
		display_info("LevelsNumber and SpeciesNumber NOT consistent (they must satisfy: LevelsNumber * SpeciesNumber <= 6). Exiting...\n");
		error_flag = 10;
	}
	error_flag_check();

	//check only one run is selected
/*	if(dilution+zoom+cascade != 1){
		display_info("Undetermined action to be done! (ONLY ONE among dilution, zoom and cascade must be set to 1, the others must be 0). Exiting...\n");
		error_flag = 10;
	}
*/
	dilution=false; zoom=false; cascade=false;
	if(strcmp(run_type,"dilution") == 0) dilution=true;
	else if(strcmp(run_type,"zoom") == 0) zoom=true;
	else if(strcmp(run_type,"cascade") == 0) cascade=true;
	else{
		display_info("Undetermined action to be done! (allowed values for RunType are: dilution, zoom and cascade). Exiting...\n");
                error_flag = 10;
	}
	error_flag_check();

	//check positiveness of radii
	for(i=0; i<levels_number; i++) if(level_bubbles_radii[i] <= 0){
		display_info("Bubbles radii MUST be positive!");
		error_flag = 10;
	}
	error_flag_check();

	//check higher levels have lower radii
	for(i=0; i<levels_number-1; i++) if(level_bubbles_radii[i] < level_bubbles_radii[i+1]){
		display_info("Bubbles radii NOT decreasing!");
		error_flag = 10;
	}
	error_flag_check();

	//check random seed
	if(cascade_random_seed == -1) cascade_random_seed = time(NULL);
	if(cascade_random_seed < -1){
		display_info("cascade random seed MUST be >= -1");
		error_flag = 10;
	}
	error_flag_check();

	//check region side is >= 0
	if(level_cubic_region_side[i] < 0){
		display_info("CubesPerSide for each level MUST be >= 0 ");
		error_flag = 10;
	}
	error_flag_check();

	//check file type is supported
	if (Lin_ftype < 0 || Lin_ftype > 2){
		display_info("LowResICsFileType MUST be 1 or 2");
		error_flag = 10;
	}
	if (Hin_ftype < 0 || Hin_ftype > 2){
		display_info("HighResICsFileType MUST be 1 or 2");
		error_flag = 10;
	}
	if (Sin_ftype < 0 || Sin_ftype > 2){
		display_info("SnapFileType MUST be 1 or 2");
		error_flag = 10;
	}
	if (out_ftype < 0 || out_ftype > 2){
		display_info("OutputFileType MUST be 1 or 2");
		error_flag = 10;
	}
	error_flag_check();
}


//The following function finalize the run
void finalize(){

	if(my_rank == 0) fclose(info);
  	MPI_Finalize();	// Finalize the MPI environment.
}


//The following function is used to check if there are particle outside ALL the regions defined
//this is simply done by checking level0 since all the other levels are contained in it
void check_ignored_particles(){
	Sphere biggest_sphere(c[0], c[1], c[2], level_bubbles_radii[0]);
	float vertex[8][3] = {
		{ BoxSize, BoxSize, BoxSize },
		{ BoxSize, BoxSize, 0.0f },
		{ BoxSize, 0.0f, BoxSize },
		{ 0.0f, BoxSize, BoxSize },
		{ BoxSize, 0.0f, 0.0f },
		{ 0.0f, BoxSize, 0.0f },
		{ 0.0f, 0.0f, BoxSize },
		{ 0.0f, 0.0f, 0.0f }
	};

	//Print a warning if there are particle ignored
	int i;
	bool warning = false;
	for(i=0; i<8; i++) if(!biggest_sphere.contains(vertex[i])) warning = true;
	if(warning){
		display_info("WARNING: your spheres does NOT contain all the simulation box, some particles may be ignored!\n");
		print_info("WARNING: your spheres does NOT contain all the simulation box, some particles may be ignored!\n");
	}
}
