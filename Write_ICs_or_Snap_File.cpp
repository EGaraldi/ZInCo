/*============================================================

* File Name : Write_ICs_or_Snap_File.cpp

* Purpose : This function write a new ICs or SNAPSHOT file, given the file
			name, a header, a particles array and a list of things to write.
			The IDs are produced in a squential way.

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include"global_variables.h"
#include"functions_declaration.h"

void WriteBlockLabel(FILE* file, char label[4], int blocksize){
	int blockheader = 4*sizeof(char) + sizeof(int);
	my_fwrite(&blockheader, sizeof(int), 1, file);
	my_fwrite(&label, sizeof(char), 4, file);
	my_fwrite(&blocksize, sizeof(int), 1, file);
	my_fwrite(&blockheader, sizeof(int), 1, file);
}

void Write_ICs_or_Snap_File(char fname[], int ftype, io_header& header, vector< vector<particle_data> >& particles, bool things_to_write[]){
	//things_to_write is an array of bbol corresponding to: header, pos, vel, id, mass. The only optional is the last, all the other 
	//must be true to have a working ICs file. Nevertheless, it is used this way to allow an easy extension to SPH variables and
	//to allow to produce non-working ICs file, whatever the reason.

	//first some variables used everywhere
	int tot_npart=0;
	int i, j, dim;

	FILE * file;
	file=fopen(fname,"w");
	if (file == NULL){ display_info("ERROR while opening file %s\n", fname); error_flag=7; }
	error_flag_check();

	//The structure of the function is simple: for each true item in thigs_to_write the corresponding
	//blockheader+block+blockheader is printed

	//header
	if (things_to_write[0]){
		dim = 256;
		if (ftype == 2) WriteBlockLabel(file, "HEAD", dim + 2*sizeof(int));
		my_fwrite(&dim, sizeof(int), 1, file);
		my_fwrite(&header, sizeof(io_header), 1, file);
		my_fwrite(&dim, sizeof(int), 1, file);
	}

	#ifdef VDEBUG
	printf("proc %i   debug 6.5\n",my_rank);
		#endif

	//the total number of particles in this file is needed everywhere, so:
	for(i=0;i<6;i++) tot_npart += header.npart[i];

	//positions
	if(things_to_write[1]){
		#ifdef NO_CENTERING_OUTPUT
		//particles back to their original position
		for(i=0; i<particles.size(); i++) for(j=0;j<particles[i].size();j++){
			particles[i][j].pos[0] += shift[0];
			if(particles[i][j].pos[0]<0) particles[i][j].pos[0]+=BoxSize; if(particles[i][j].pos[0]>BoxSize) particles[i][j].pos[0] -= BoxSize;
			particles[i][j].pos[1] += shift[1];
			if(particles[i][j].pos[1]<0) particles[i][j].pos[1]+=BoxSize; if(particles[i][j].pos[1]>BoxSize) particles[i][j].pos[1] -= BoxSize;
			particles[i][j].pos[2] += shift[2];
			if(particles[i][j].pos[2]<0) particles[i][j].pos[2]+=BoxSize; if(particles[i][j].pos[2]>BoxSize) particles[i][j].pos[2] -= BoxSize;
		}
		#endif

		dim = tot_npart*sizeof(float)*3;
		if (ftype == 2) WriteBlockLabel(file, "POS ", dim + 2*sizeof(int));
		my_fwrite(&dim, sizeof(int), 1, file);
		//write the last particles in the array in order to make the erase of the written element faster (no need to shift the following elements
		//to the positions of the cancelled ones)
		for(i=0;i<6;i++) for(j=1;j<=header.npart[i];j++) my_fwrite(&particles[i][particles[i].size()-j].pos[0], sizeof(float), 3, file);
		my_fwrite(&dim, sizeof(int), 1, file);

		#ifdef VDEBUG
		printf("proc %i   debug 6.6\n",my_rank);
		#endif
	}

	//velocities
	if(things_to_write[2]){
		dim = tot_npart*sizeof(float)*3;
		if (ftype == 2) WriteBlockLabel(file, "VEL ", dim + 2*sizeof(int));
		my_fwrite(&dim, sizeof(int), 1, file);
		for(i=0;i<6;i++) for(j=1;j<=header.npart[i];j++) my_fwrite(&particles[i][particles[i].size()-j].vel[0], sizeof(float), 3, file);
		my_fwrite(&dim, sizeof(int), 1, file);

		#ifdef VDEBUG
		printf("proc %i   debug 6.7\n",my_rank);
		#endif
		}

	//IDs
	if(things_to_write[3]){
		dim = tot_npart*sizeof(LOIout);
		if (ftype == 2) WriteBlockLabel(file, "IDS ", dim + 2*sizeof(int));
		my_fwrite(&dim, sizeof(int), 1, file);
		//the IDs are sequentially produced
		for(i=0;i<6;i++) for(j=1;j<=header.npart[i];j++){ my_fwrite(&lastID, sizeof(LOIout), 1, file); lastID++; }
		my_fwrite(&dim, sizeof(int), 1, file);

		#ifdef VDEBUG
		printf("proc %i   debug 6.8\n",my_rank);
		#endif
	}

	//masses
	if(things_to_write[4]){
		dim = tot_npart*sizeof(float);
		if (ftype == 2) WriteBlockLabel(file, "MASS", dim + 2*sizeof(int));
		my_fwrite(&dim, sizeof(int), 1, file);
		for(i=0;i<6;i++) for(j=1;j<=header.npart[i];j++) my_fwrite(&particles[i][particles[i].size()-j].mass, sizeof(float), 1, file);
		my_fwrite(&dim, sizeof(int), 1, file);
	}

	//fake internal energies
	#ifdef ADD_INTERNAL_ENERGIES
	if(header.npart[0] > 0){
		float ie=0.0f;
		dim = header.npart[0]*sizeof(float);
		if (ftype == 2) WriteBlockLabel(file, "U   ", dim + 2*sizeof(int));
        my_fwrite(&dim, sizeof(int), 1, file);
        for(j=0;j<header.npart[0];j++) my_fwrite(&ie, sizeof(float), 1, file);
        my_fwrite(&dim, sizeof(int), 1, file);
	}
	#endif

	fclose(file);

	#ifdef PRINT_IDS_CORRESPONDENCE
	int k;
	for(i=0;i<6;i++) for(j=1;j<=header.npart[i];j++){
		fprintf(id_corr, "%zu", lastID);
		for(k=0; k<old_ids[i][particles[i].size()-j].size(); k++) fprintf(id_corr, "\t\t%zu\n", old_ids[i][particles[i].size()-j][k]);
	}
	#endif

	#ifdef VDEBUG
	printf("proc %i   debug 6.9\n",my_rank);
	#endif
}
