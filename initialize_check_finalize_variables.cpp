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

#include<omp.h>
#include"global_variables.h"
#include"functions_declaration.h"

//The following function is called at the beginning of the run and deals (mainly) with:
// -initializing MPI
// -building MPI custom datatypes
// -checking the call sequence
void initialize(int argc){

	error_flag = 0;
	lastID=0;

	// Get the number of processes and the current process' rank
	world_size = omp_get_num_threads();
	my_rank = omp_get_thread_num();

	//starting the clock
	init_time = omp_get_wtime();

	//check the call sequence
	if(argc<2){ display_info("Missing parameters file, please run using: <exec> <parameters file name>\n"); error_flag=1; }
	error_flag_check();

	//the vector containing the new particles created is a 6*? vector, since GADGET allows (max) 6 particles types
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
