/*============================================================

* File Name : main.cpp

* Purpose : main function for the ZInCo code

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include"global_variables.h"
#include"functions_declaration.h"

int main(int argc, char** argv){

	#ifdef VDEBUG
	printf("proc ?   debug 0\n");
	#endif

	//initialize everything (mainly MPI)
	initialize(argc);

	//read parameters from file
	display_info("Reading the parameters file...\n");
	char ParameterFile[200];
	strcpy(ParameterFile, argv[1]);
	Read_parameter_file(ParameterFile);

	//check everything is fine with the read parameters
	check_variables();

	#ifdef VDEBUG
	printf("proc %i   debug 1\n",my_rank);
	#endif

	//open the file where printing information
	if(my_rank == 0){
		char filename_info[200];
		sprintf(filename_info,"%s/info.txt",output_dir);
		info=fopen(filename_info,"w");
		if (info == NULL){ display_info("ERROR while opening file %s\n", filename_info); error_flag=2; }
	}
	error_flag_check();	

	//print header
	print_info("// ================================================================================================= //\n");
	print_info("//                                                ZInCo                                              //\n");
	print_info("// ================================================================================================= //\n\n");
	if(dilution) print_info("  This is a DILUTION  \n\n");
	else if(zoom) print_info("  This is a ZOOM  \n\n");
	else if(cascade) print_info("  This is a CASCADE  \n\n");
	display_info("Starting...\n");

	#ifdef VDEBUG
	printf("proc %i   debug 2\n",my_rank);
	#endif

	#if defined(VDEBUG) || defined(DEBUG)
	MPI_Barrier(MPI_COMM_WORLD);
	#endif

	//get info from ICs header
	print_info("Reading the header to get BoxSize and compute mtot, ntot, ltot...\n");
	display_info("Reading the header to get BoxSize and compute mtot, ntot, ltot...\n");

	build_file_name(Hic_dir, Hic_name, Hin_fnr, 0, fname);

	bool things_to_read[6] = { true, false, false, false, false, false };
	Read_ICs_File(fname, Hin_ftype, header1, particles_in, idH, false, things_to_read);

	BoxSize=(float)header1.BoxSize;

	#ifdef VDEBUG
	printf("proc %i   debug 3\n",my_rank);
	#endif

	//send BoxSize to all the procs
	MPI_Bcast(&BoxSize,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	//print a warning for ignored particles (i.e. particles outside *all* the resolution regions specified)
	check_ignored_particles();

	#if defined(VDEBUG) || defined(DEBUG)
	MPI_Barrier(MPI_COMM_WORLD);
	#endif

	//initialize the desired run
	if(dilution) create_dilution();
	if(zoom) create_zoom();
	if(cascade) create_cascade();

	//for the zoom run: from now on, whenever a particle is read, its position is shifted in order to have the center (c) exactly at the center of the box
	//(this is useful becaus in GADGET2 the load balancing is far easier if the high-res region is connected)
	if(zoom){
		shift[0]=c[0]-BoxSize/2;
		shift[1]=c[1]-BoxSize/2;
		shift[2]=c[2]-BoxSize/2;
	} else{
		shift[0]=0.0f;
		shift[1]=0.0f;
		shift[2]=0.0f;
	}

	#ifdef DEBUG
	display_info("shift= %f %f %f\n",shift[0],shift[1],shift[2]);
	#endif

	//compute the number of cells per side
	//(in GADGET the simlation box is cubic, so mtot=ntot=ltot. They are different only to ease possible future developments)
	lambda=(float) BoxSize/cubes_per_side;
	mtot=cubes_per_side;
	ntot=cubes_per_side;
	ltot=cubes_per_side;

	#ifdef DEBUG
	printf("(process %i) BoxSize:%f  cubes_per_side:%i  lambda:%f\n",my_rank,BoxSize,cubes_per_side,lambda);
	#endif

	print_info("Time elapsed: %f\n\n",MPI_Wtime()-init_time);

	//if we are zooming, we have to adjust the resolution regions in the ICs. This is done by tracking the particles from the region defined 
	//in the snapshot to the corresponding one in the corresponding ICs
	#ifndef ZOOM_FROM_ICS
	if(zoom){
		print_info("Calling New_regions_finder...\n");
		display_info("Calling New_regions_finder...\n");

		New_regions_finder();
		print_info("Time elapsed: %f\n\n",MPI_Wtime()-init_time);
		
		#ifdef VDEBUG
		printf("proc %i   debug 4\n",my_rank);
		#endif
	}
	#endif

	#ifdef DEBUG
	printf("(process %i) r_high_IC:%f  r_medium_IC:%f  c_IC:%f  %f  %f\n",my_rank,r_high_IC,r_medium_IC,c_IC[0],c_IC[1],c_IC[2]);
	#endif

	#if defined(VDEBUG) || defined(DEBUG)
	MPI_Barrier(MPI_COMM_WORLD);
	#endif

	//Now we slice the box and assign each cell to the corresponding resolution region
	print_info("Calling Slicer...\n");
	display_info("Calling Slicer...\n");
	Slicer();
	print_info("Time elapsed: %f\n\n",MPI_Wtime()-init_time);

	#ifdef DEBUG
	//print out the resolution matrix
	display_info("(process %i) %i %i %i res:\n",my_rank,mtot,ntot,ltot);
	for(int l=0;l<ltot;l++){
		for(int n=0;n<ntot;n++){
			for(int m=0;m<mtot;m++) display_info("%i ",resolution_matrix[m][n][l].cubeID);
			display_info("\n");
		}
		display_info("\n");
	}
	display_info("(process %i) %i %i %i res:\n",my_rank,mtot,ntot,ltot);
	for(int l=0;l<ltot;l++){
		for(int n=0;n<ntot;n++){
			for(int m=0;m<mtot;m++) display_info("%i ",resolution_matrix[m][n][l].level);
			display_info("\n");
		}
		display_info("\n");
	}
	#endif	
	
	#ifdef VDEBUG
	printf("proc %i   debug 5\n",my_rank);
	#endif

	#if defined(VDEBUG) || defined(DEBUG)
	MPI_Barrier(MPI_COMM_WORLD);
	#endif

	//process the cells and produce the particles for the ICs
	print_info("Calling New_particles_maker...\n");
	display_info("Calling New_particles_maker...\n");
	New_particles_maker();
	print_info("Time elapsed: %f\n\n",MPI_Wtime()-init_time);

	#ifdef VDEBUG
	printf("proc %i   debug 6\n",my_rank);
	#endif

	#if defined(VDEBUG) || defined(DEBUG)
	MPI_Barrier(MPI_COMM_WORLD);
	#endif

	//write the new ICs file(s)
	print_info("Calling New_files_writer...\n");
	display_info("Calling New_files_writer...\n");
	New_files_writer();
	print_info("Time elapsed: %f\n\n",MPI_Wtime()-init_time);

	MPI_Barrier(MPI_COMM_WORLD);
	
	display_info("Done!\n");

	#ifdef VDEBUG
	printf("proc %i   debug 7\n",my_rank);
	#endif

	//finalize everything and exit
	finalize();
	return 0;
}
