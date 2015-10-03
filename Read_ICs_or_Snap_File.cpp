/*============================================================

* File Name : Read_ICs_or_Snap_File.cpp

* Purpose : Two functions to read the data from an ICs or SNAPSHOT file.
			The information to read are defined in an array (things_to_read).
			The data can be BCast to other procs.

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include"global_variables.h"
#include"functions_declaration.h"


bool ReadBlockLabel(){
	int blockheader1, blockheader2, blocksize;
	char label[4];
	my_fread(&blockheader1, sizeof(int), 1, file);
	my_fread(&label, sizeof(char), 4, file);
	my_fread(&blocksize, sizeof(int), 1, file);
	my_fread(&blockheader1, sizeof(int), 1, file);
	return blockheader1 == blockheader2
}

void Read_ICs_File(char fname[], int ftype,  io_header& header, vector<particle_data>& particles, vector<LOIinHigh>& id, bool share_data, bool things_to_read[]){

/*  share_data = true if data must be sent to other process
    things_to_read[i] = true if the field must be read
        fields: header, pos, vel, id, masses                    */

	//we MUST read the header to have information on the content of the file.
	if(!things_to_read[0]) return;

	int blockheader1, blockheader2;
	int N, Nwithmass;
	int k;

	//proc 0 read and then send data to the others if share_data == true
	//this is faster than parallel reading
	if(my_rank == 0){

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.1\n",my_rank);
		#endif

		//open the file to read
		FILE *file;
		file=fopen(fname,"r");
		if (file == NULL){ display_info("ERROR while opening file %s\n", fname); error_flag=3; }

		error_flag_check();

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.2\n",my_rank);
		#endif

		//read the header
		if (ftype == 2)
			if (!ReadBlockLabel()){
				printf("fatal error: blockheaders not equal! Stopped at: <header> label.\n");
				error_flag = 4;
				error_flag_check();
			}
		my_fread(&blockheader1, sizeof(int), 1, file);
		my_fread(&header, sizeof(io_header), 1, file);
		my_fread(&blockheader2, sizeof(int), 1, file);
		if(blockheader1 != blockheader2){
			printf("fatal error: blockheaders not equal! Stopped at: <header> block.\n");
			error_flag=4;
		}

		error_flag_check();

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.3\n",my_rank);
		#endif

		N=0;	//total number of particle in this file
		for(k=0; k<6; k++) N += header.npart[k];

		//allocate memory for the particles to read
		if (things_to_read[1] || things_to_read[2] || things_to_read[3] || things_to_read[4])
			particles.resize(N);
		else return;

		//positions
		
		#ifdef DEBUG
		printf("N: %i   particles.size: %zu\n",N,particles.size());
		#endif

		if (ftype == 2)
			if (!ReadBlockLabel()){
				printf("fatal error: blockheaders not equal! Stopped at: <position> label.\n");
				error_flag = 4;
				error_flag_check();
			}
		my_fread(&blockheader1, sizeof(int), 1, file);

		//in the zoom run, we shift the positions in order to have the center of the zoom regions at the centre of the box
		if (things_to_read[1]){
			for(k=0;k<N;k++){
				my_fread(&particles[k].pos[0], sizeof(float), 3, file);
				particles[k].pos[0] -= shift[0];
				if(particles[k].pos[0]<0) particles[k].pos[0]+=BoxSize; if(particles[k].pos[0]>BoxSize) particles[k].pos[0] -= BoxSize;
				particles[k].pos[1] -= shift[1];
				if(particles[k].pos[1]<0) particles[k].pos[1]+=BoxSize; if(particles[k].pos[1]>BoxSize) particles[k].pos[1] -= BoxSize;
				particles[k].pos[2] -= shift[2];
				if(particles[k].pos[2]<0) particles[k].pos[2]+=BoxSize; if(particles[k].pos[2]>BoxSize) particles[k].pos[2] -= BoxSize;
			}
		} else{
			float trash;
			my_fread(&trash, sizeof(float), 3*N, file);
		}
		my_fread(&blockheader2, sizeof(int), 1, file);
		if(blockheader1 != blockheader2){
			display_info("fatal error: blockheaders not equal! Stopped at: <position> block.\n");
			error_flag=4;
		}

		error_flag_check();
		

		//velocities

		if (!things_to_read[2] && !things_to_read[3] && !things_to_read[4]) return;

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.4\n",my_rank);
		#endif

		if (ftype == 2)
			if (!ReadBlockLabel()){
				printf("fatal error: blockheaders not equal! Stopped at: <velocities> label.\n");
				error_flag = 4;
				error_flag_check();
			}
		my_fread(&blockheader1, sizeof(int), 1, file);

		if(things_to_read[2])
			for(k=0;k<N;k++) my_fread(&particles[k].vel[0], sizeof(float), 3, file);
		else{
			float trash;
			my_fread(&trash, sizeof(float), 3 * N, file);
		}

		my_fread(&blockheader2, sizeof(int), 1, file);

		if(blockheader1 != blockheader2){
			display_info("fatal error: blockheaders not equal! Stopped at: <velocities> block.\n");
			error_flag=4;
		}
		error_flag_check();


		//IDs

		if (!things_to_read[3] && !things_to_read[4]) return;

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.5\n",my_rank);
		#endif

		id.resize(N);

		if (ftype == 2)
			if (!ReadBlockLabel()){
				printf("fatal error: blockheaders not equal! Stopped at: <ids> label.\n");
				error_flag = 4;
				error_flag_check();
			}
		my_fread(&blockheader1, sizeof(int), 1, file);

		if(things_to_read[3])
			for(k=0;k<N;k++) {/*id.push_back(-1);*/ my_fread(&id[k], sizeof(LOIinHigh), 1, file);}
		else{
			LOIinHigh trash;
			my_fread(&trash, sizeof(float), N, file);
		}
		my_fread(&blockheader2, sizeof(int), 1, file);
		if(blockheader1 != blockheader2){
			display_info("fatal error: blockheaders not equal! Stopped at: <ids> block.\n");
			error_flag=4;
		}

		error_flag_check();

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.6\n",my_rank);
		#endif

		//From now on process particles type-by-type
		Ngas=header.npart[0];
		Nhalo=header.npart[1];
		Ndisk=header.npart[2];
		Nbulge=header.npart[3];
		Nstars=header.npart[4];
		Nbndry=header.npart[5];


		//masses

		if (!things_to_read[4]) return;

		//Determine the number of masses to be read in order to determine if the block (=> blockeheader) is present
		Nwithmass= 0;
	    for(k=0;k<6;k++){
        	if((header.npart[k] > 0) && (header.massarr[k] == 0.0f)) Nwithmass += header.npart[k];
        }

		if (Nwithmass > 0){
			if (ftype == 2)
				if (!ReadBlockLabel()){
					printf("fatal error: blockheaders not equal! Stopped at: <masses> label.\n");
					error_flag = 4;
					error_flag_check();
				}
			my_fread(&blockheader1, sizeof(int), 1, file);
		}

		if(things_to_read[4]){
			//For each species present, if the corresponding massarr is zero the corresponding masses are read, otherwise their mass is set to the one in massarr
			if(Ngas>0){
				if(header.massarr[0] == 0.0f) for(k=0;k<Ngas;k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for(k=0;k<Ngas;k++) particles[k].mass=header.massarr[0];
			}
			if(Nhalo>0){
				if(header.massarr[1] == 0.0f) for(k=Ngas;k<Ngas+Nhalo;k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for(k=Ngas;k<Ngas+Nhalo;k++) particles[k].mass=header.massarr[1];
			}
			if(Ndisk>0){
				if(header.massarr[2] == 0.0f) for(k=Ngas+Nhalo;k<Ngas+Nhalo+Ndisk;k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for(k=Ngas+Nhalo;k<Ngas+Nhalo+Ndisk;k++) particles[k].mass=header.massarr[2];
			}
			if(Nbulge>0){
				if(header.massarr[3] == 0.0f) for(k=Ngas+Nhalo+Ndisk;k<Ngas+Nhalo+Ndisk+Nbulge;k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for(k=Ngas+Nhalo+Ndisk;k<Ngas+Nhalo+Ndisk+Nbulge;k++) particles[k].mass=header.massarr[3];
			}
			if(Nstars>0){
				if(header.massarr[4] == 0.0f) for(k=Ngas+Nhalo+Ndisk+Nbulge;k<Ngas+Nhalo+Ndisk+Nbulge+Nstars;k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for(k=Ngas+Nhalo+Ndisk+Nbulge;k<Ngas+Nhalo+Ndisk+Nbulge+Nstars;k++) particles[k].mass=header.massarr[4];
			}
			if(Nbndry>0){
				if(header.massarr[5] == 0.0f) for(k=Ngas+Nhalo+Ndisk+Nbulge+Nstars;k<Ngas+Nhalo+Ndisk+Nbulge+Nstars+Nbndry;k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for(k=Ngas+Nhalo+Ndisk+Nbulge+Nstars;k<Ngas+Nhalo+Ndisk+Nbulge+Nstars+Nbndry;k++) particles[k].mass=header.massarr[5];
			}
		} else{
			float trash;
			my_fread(&trash, sizeof(float), Nwithmass, file);
		}

		if (Nwithmass > 0) my_fread(&blockheader2, sizeof(int), 1, file);

		if(blockheader1 != blockheader2){
			display_info("fatal error: blockheaders not equal! Stopped at: <masses> block.\n");
			error_flag=4;
		}

		error_flag_check();

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.7\n",my_rank);
		#endif

		//end_function:

		fclose(file);
	} else{
		if(!things_to_read[0]) return;
		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.1\n",my_rank);
		#endif
		error_flag_check();

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.2\n",my_rank);
		#endif
		if (ftype == 2)	error_flag_check();
		error_flag_check();

		if(things_to_read[1]){
			#ifdef VDEBUG
			printf("proc %i   debug RIoSF.3\n",my_rank);
			#endif
			if (ftype == 2)	error_flag_check();
			error_flag_check();
		}

		if(things_to_read[2]){
			#ifdef VDEBUG
			printf("proc %i   debug RIoSF.4\n",my_rank);
			#endif
			if (ftype == 2)	error_flag_check();
			error_flag_check();
		}

		if(things_to_read[3]){
			#ifdef VDEBUG
			printf("proc %i   debug RIoSF.5\n",my_rank);
			#endif
			if (ftype == 2)	error_flag_check();
			error_flag_check();
		}

		if(things_to_read[4]){
			#ifdef VDEBUG
			printf("proc %i   debug RIoSF.6\n",my_rank);
			#endif
			if (ftype == 2)	error_flag_check();
			error_flag_check();
		}

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.7\n",my_rank);
		#endif

		//end_function:
	}

	if(share_data){ //BCast information to other procs

		#if defined(DEBUG) || defined(VDEBUG) 
		printf("proc %i at MPI_Barrier\n",my_rank);
		MPI_Barrier(MPI_COMM_WORLD);
		#endif


		//header always read => always BCast
		MPI_Bcast(&blockheader1, 1, MPI_INT, 0,MPI_COMM_WORLD);
		MPI_Bcast(&header, 1, MPI_HEADER_t, 0,MPI_COMM_WORLD);
		MPI_Bcast(&blockheader2, 1, MPI_INT, 0,MPI_COMM_WORLD);			
		MPI_Bcast(&N, 1, MPI_INT, 0,MPI_COMM_WORLD);
		//update variables
		Ngas=header.npart[0];
		Nhalo=header.npart[1];
		Ndisk=header.npart[2];
		Nbulge=header.npart[3];
		Nstars=header.npart[4];
		Nbndry=header.npart[5];

		//positions
		if(things_to_read[1]){
			#if defined(DEBUG) || defined(VDEBUG)  
			printf("proc %i at MPI_Barrier\n",my_rank);
			MPI_Barrier(MPI_COMM_WORLD);
			#endif

			//resize array to contain the informations
			if(my_rank != 0){
				particles.resize(N);
				if(things_to_read[3]) id.resize(N);
			}

			#ifdef DEBUG
			printf("proc %i  N= %i particles.size= %zu  id.size= %zu\n",my_rank,N,particles.size(),id.size());
			#endif

			for(k=0;k<N;k++){
				MPI_Bcast(&particles[k], 1, MPI_PARTICLE_t, 0, MPI_COMM_WORLD);
				if(things_to_read[3]) MPI_Bcast(&id[k], 1, MPI_LOIinHigh, 0,MPI_COMM_WORLD);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		#ifdef DEBUG
		display_info("End broadcasting\n");
		#endif
	}
}

//same as above but for snapshot file
void Read_Snap_File(char fname[], int ftype, io_header& header, vector<particle_data>& particles, vector<LOIinLow>& id, bool share_data, bool things_to_read[]){

/*  share_data = true if data must be sent to other process
    things_to_read[i] = true if the field must be read
        fields: header, pos, vel, id, masses                    */

	if(!things_to_read[0]) return;

	int blockheader1, blockheader2;
	int N, Nwithmass;
	int k;

	if(my_rank == 0){

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.1\n",my_rank);
		#endif

		FILE *file;
		file=fopen(fname,"r");
		if (file == NULL){ display_info("ERROR while opening file %s\n", fname); error_flag=3; }

		error_flag_check();

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.2\n",my_rank);
		#endif

		if (ftype == 2)
			if (!ReadBlockLabel()){
				printf("fatal error: blockheaders not equal! Stopped at: <header> label.\n");
				error_flag = 4;
				error_flag_check();
			}
		my_fread(&blockheader1, sizeof(int), 1, file);
		my_fread(&header, sizeof(io_header), 1, file);
		my_fread(&blockheader2, sizeof(int), 1, file);
		if(blockheader1 != blockheader2){
			printf("fatal error: blockheaders not equal! Stopped at: <header> block.\n");
			error_flag=4;
		}

		error_flag_check();

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.3\n",my_rank);
		#endif

		N=0;	//total number of particle in this file
		for(k=0; k<6; k++) N += header.npart[k];

		if (things_to_read[1] || things_to_read[2] || things_to_read[3] || things_to_read[4])
			particles.resize(N);
		else return;

		#ifdef DEBUG
		printf("N: %i   particles.size: %zu\n",N,particles.size());
		#endif

		if (ftype == 2)
			if (!ReadBlockLabel()){
				printf("fatal error: blockheaders not equal! Stopped at: <positions> label.\n");
				error_flag = 4;
				error_flag_check();
			}
		my_fread(&blockheader1, sizeof(int), 1, file);

		if (things_to_read[1]){
			for (k = 0; k<N; k++){
				my_fread(&particles[k].pos[0], sizeof(float), 3, file);
				particles[k].pos[0] -= shift[0];
				if(particles[k].pos[0]<0) particles[k].pos[0]+=BoxSize; if(particles[k].pos[0]>BoxSize) particles[k].pos[0] -= BoxSize;
				particles[k].pos[1] -= shift[1];
				if(particles[k].pos[1]<0) particles[k].pos[1]+=BoxSize; if(particles[k].pos[1]>BoxSize) particles[k].pos[1] -= BoxSize;
				particles[k].pos[2] -= shift[2];
				if(particles[k].pos[2]<0) particles[k].pos[2]+=BoxSize; if(particles[k].pos[2]>BoxSize) particles[k].pos[2] -= BoxSize;
			} 
		}else{
			float trash;
			my_fread(&trash, sizeof(float), 3 * N, file);
		}

		my_fread(&blockheader2, sizeof(int), 1, file);
			
		if(blockheader1 != blockheader2){
			display_info("fatal error: blockheaders not equal! Stopped at: <position> block.\n");
			error_flag=4;
		}

		error_flag_check();
		


		if (!things_to_read[2] && !things_to_read[3] && !things_to_read[4]) return;

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.4\n",my_rank);
		#endif

		my_fread(&blockheader1, sizeof(int), 1, file);
		if (things_to_read[2])
			for(k=0;k<N;k++) my_fread(&particles[k].vel[0], sizeof(float), 3, file);
		else{
			float trash;
			my_fread(&trash, sizeof(float), 3 * N, file);
		}

		if (ftype == 2)
			if (!ReadBlockLabel()){
				printf("fatal error: blockheaders not equal! Stopped at: <velocities> label.\n");
				error_flag = 4;
				error_flag_check();
			}
		my_fread(&blockheader2, sizeof(int), 1, file);
		if(blockheader1 != blockheader2){
			display_info("fatal error: blockheaders not equal! Stopped at: <velocities> block.\n");
			error_flag=4;
		}
		error_flag_check();

		if (!things_to_read[3] && !things_to_read[4]) return;

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.5\n",my_rank);
		#endif

		id.resize(N);

		if (ftype == 2)
			if (!ReadBlockLabel()){
				printf("fatal error: blockheaders not equal! Stopped at: <ids> label.\n");
				error_flag = 4;
				error_flag_check();
			}
		my_fread(&blockheader1, sizeof(int), 1, file);

		if (things_to_read[3])
			for(k=0;k<N;k++) {/*id.push_back(-1);*/ my_fread(&id[k], sizeof(LOIinLow), 1, file);}
		else{
			LOIinLow trash;
			my_fread(&trash, sizeof(float), N, file);
		}
			
		my_fread(&blockheader2, sizeof(int), 1, file);
			
		if(blockheader1 != blockheader2){
			display_info("fatal error: blockheaders not equal! Stopped at: <ids> block.\n");
			error_flag=4;
		}

		error_flag_check();

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.6\n",my_rank);
		#endif
		//process particles type-by-type
		Ngas=header.npart[0];
		Nhalo=header.npart[1];
		Ndisk=header.npart[2];
		Nbulge=header.npart[3];
		Nstars=header.npart[4];
		Nbndry=header.npart[5];

		if (!things_to_read[4]) return;

		//Determine the number of masses to be read in order to determine if the block (=> blockeheader) is present
		Nwithmass = 0;
		for (k = 0; k<6; k++){
			if ((header.npart[k] > 0) && (header.massarr[k] == 0.0f)) Nwithmass += header.npart[k];
		}

		if (Nwithmass > 0) my_fread(&blockheader1, sizeof(int), 1, file);

		if (things_to_read[4]){
			//For each species present, if the corresponding massarr is zero the corresponding masses are read, otherwise their mass is set to the one in massarr
			if (Ngas>0){
				if (header.massarr[0] == 0.0f) for (k = 0; k<Ngas; k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for (k = 0; k<Ngas; k++) particles[k].mass = header.massarr[0];
			}
			if (Nhalo>0){
				if (header.massarr[1] == 0.0f) for (k = Ngas; k<Ngas + Nhalo; k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for (k = Ngas; k<Ngas + Nhalo; k++) particles[k].mass = header.massarr[1];
			}
			if (Ndisk>0){
				if (header.massarr[2] == 0.0f) for (k = Ngas + Nhalo; k<Ngas + Nhalo + Ndisk; k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for (k = Ngas + Nhalo; k<Ngas + Nhalo + Ndisk; k++) particles[k].mass = header.massarr[2];
			}
			if (Nbulge>0){
				if (header.massarr[3] == 0.0f) for (k = Ngas + Nhalo + Ndisk; k<Ngas + Nhalo + Ndisk + Nbulge; k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for (k = Ngas + Nhalo + Ndisk; k<Ngas + Nhalo + Ndisk + Nbulge; k++) particles[k].mass = header.massarr[3];
			}
			if (Nstars>0){
				if (header.massarr[4] == 0.0f) for (k = Ngas + Nhalo + Ndisk + Nbulge; k<Ngas + Nhalo + Ndisk + Nbulge + Nstars; k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for (k = Ngas + Nhalo + Ndisk + Nbulge; k<Ngas + Nhalo + Ndisk + Nbulge + Nstars; k++) particles[k].mass = header.massarr[4];
			}
			if (Nbndry>0){
				if (header.massarr[5] == 0.0f) for (k = Ngas + Nhalo + Ndisk + Nbulge + Nstars; k<Ngas + Nhalo + Ndisk + Nbulge + Nstars + Nbndry; k++) my_fread(&particles[k].mass, sizeof(float), 1, file);
				else for (k = Ngas + Nhalo + Ndisk + Nbulge + Nstars; k<Ngas + Nhalo + Ndisk + Nbulge + Nstars + Nbndry; k++) particles[k].mass = header.massarr[5];
			}
		} else{
			float trash;
			my_fread(&trash, sizeof(float), Nwithmass, file);
		}


		if (Nwithmass > 0){
			if (ftype == 2)
				if (!ReadBlockLabel()){
					printf("fatal error: blockheaders not equal! Stopped at: <masses> label.\n");
					error_flag = 4;
					error_flag_check();
				}
			my_fread(&blockheader1, sizeof(int), 1, file);
		}

		if (blockheader1 != blockheader2){
			display_info("fatal error: blockheaders not equal! Stopped at: <masses> block.\n");
			error_flag = 4;
		}

		error_flag_check();

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.7\n",my_rank);
		#endif

		//end_function:

		fclose(file);
	} else{
		if(!things_to_read[0]) return;
		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.1\n",my_rank);
		#endif
		error_flag_check();

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.2\n",my_rank);
		#endif
		if (ftype == 2)	error_flag_check();
		error_flag_check();

		if(things_to_read[1]){
			#ifdef VDEBUG
			printf("proc %i   debug RIoSF.3\n",my_rank);
			#endif
			if (ftype == 2)	error_flag_check();
			error_flag_check();
		}

		if(things_to_read[2]){
			#ifdef VDEBUG
			printf("proc %i   debug RIoSF.4\n",my_rank);
			#endif
			if (ftype == 2)	error_flag_check();
			error_flag_check();
		}

		if(things_to_read[3]){
			#ifdef VDEBUG
			printf("proc %i   debug RIoSF.5\n",my_rank);
			#endif
			if (ftype == 2)	error_flag_check();
			error_flag_check();
		}

		if(things_to_read[4]){
			#ifdef VDEBUG
			printf("proc %i   debug RIoSF.6\n",my_rank);
			#endif
			if (ftype == 2)	error_flag_check();
			error_flag_check();
		}

		#ifdef VDEBUG
		printf("proc %i   debug RIoSF.7\n",my_rank);
		#endif

		//end_function:
	}

	if(share_data){

		#if defined(DEBUG) || defined(VDEBUG) 
		printf("proc %i at MPI_Barrier\n",my_rank);
		MPI_Barrier(MPI_COMM_WORLD);
		#endif

		MPI_Bcast(&blockheader1, 1, MPI_INT, 0,MPI_COMM_WORLD);
		MPI_Bcast(&header, 1, MPI_HEADER_t, 0,MPI_COMM_WORLD);
		MPI_Bcast(&blockheader2, 1, MPI_INT, 0,MPI_COMM_WORLD);			
		MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

		Ngas = header.npart[0];
		Nhalo = header.npart[1];
		Ndisk = header.npart[2];
		Nbulge = header.npart[3];
		Nstars = header.npart[4];
		Nbndry = header.npart[5];


		if(things_to_read[1]){
			#if defined(DEBUG) || defined(VDEBUG)  
			printf("proc %i at MPI_Barrier\n",my_rank);
			MPI_Barrier(MPI_COMM_WORLD);
			#endif

			if(my_rank != 0){
				particles.resize(N);
				if(things_to_read[3]) id.resize(N);
			}

			#ifdef DEBUG
			printf("proc %i  N= %i particles.size= %zu  id.size= %zu\n",my_rank,N,particles.size(),id.size());
			#endif

			for(k=0;k<N;k++){
				MPI_Bcast(&particles[k], 1, MPI_PARTICLE_t, 0, MPI_COMM_WORLD);
				if(things_to_read[3]) MPI_Bcast(&id[k], 1, MPI_LOIinLow, 0,MPI_COMM_WORLD);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		#ifdef DEBUG
		display_info("End broadcasting\n");
		#endif
	}
}
