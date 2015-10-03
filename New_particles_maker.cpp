	/*============================================================

* File Name : New_particles_maker.cpp

* Purpose : The following function reads the particles in the input ICs and
			produces the particles for the new ICs by either merging, copying
			or ignoring particles, depending on the region they belong to.

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include"global_variables.h"
#include"functions_declaration.h"

void New_particles_maker(){

	#ifdef VDEBUG
	printf("proc %i   debug 5.0\n",my_rank);
	#endif

	int i, j, k, l, m, n; //loop variables
	int ps_min = -1, ps_max; //variables containing the minimum and maximum domain id for the current proc
	float chunk = 0.0f;	//how many domains each proc analyses

	#ifdef DEBUG
	printf("proc %i  world_size=%i\n",my_rank,world_size);
	#endif

	//now compute the domains range this proc will analyses
	chunk = (rcounter-hcounter-1)/world_size;
	ps_min = int(my_rank*chunk)+hcounter+1;
	//make sure the last proc analyses all the remaining domains
	ps_max = my_rank == world_size-1 ? rcounter+1 : int((my_rank+1)*chunk)+hcounter+1;

	#ifdef VDEBUG
	printf("proc %i   debug 5.1\n",my_rank);
	#endif
	
	#ifdef DEBUG
	printf("proc %i  ps_min %i  ps_max %i\n",my_rank,ps_min,ps_max);
	printf("trying to allocate: %i x %i particle_data values\n", 6-species_number, ps_max-ps_min+1);
	#endif

	//arrays that will contain the file-averaged particles for this proc
	//we use a new variables instead of new_particles because we first sum up all the relevant quantities
	//and at the end (when the total mass of the merged particles is known) perform the average. Then using
	//a new variable avoid the problem of keeping track of the partially-processed particles in new_particles.
	particle_data **CMs;
	CMs = new particle_data*[6];
	for(i=0; i<6; i++){
		CMs[i] = new particle_data[ps_max-ps_min+1];
		for(j=0; j<ps_max-ps_min+1; j++) CMs[i][j] = particle_data();
	}

	#ifdef PRINT_IDS_CORRESPONDENCE
	vector< vector< vector< LOIinHigh > > > IDs;
	IDs.resize(6-species_number);
	for(i=0; i<6-species_number; i++){
	//	IDs.push_back(vector< vector< LOIinHigh > >());
		for(j=0; j<ps_max-ps_min+1; j++) IDs[i].resize(ps_max-ps_min);//IDs[i].push_back(vector< LOIinHigh >());
	}
	#endif

	#ifdef VDEBUG
	printf("proc %i   debug 5.2\n",my_rank);
	#endif

	#ifdef DEBUG
	long int tot_part_check = 0;		//check all the particles are processed
	long double tot_mass_check1 = 0.0;	//check for the mass of the copied particles
	long double tot_mass_check2 = 0.0;	//check for the mass of the merged particles
	#endif

	//here I read the ICs files to average the properties of particles inside the cubes
	for(i=0;i<Hin_fnr;i++){

		build_file_name(Hic_dir, Hic_name, Hin_fnr, i, fname);

		print_info("New_particles_maker: processing ICs file %s\n",fname);

		#ifdef VDEBUG
		printf("proc %i   debug 5.3\n",my_rank);
		#endif

		//read "old" particles
		bool things_to_read[5] = { true, true, true, true, true }; /* header, pos, vel, id, mass */
		Read_ICs_File(fname, Hin_ftype, header1, particles_in, idH, true, things_to_read);
		#ifdef DEBUG
		printf("Ngas %i Nhalo %i Hdisk %i Nbulge %i Nstars %i Nbndry %i\n",Ngas, Nhalo, Ndisk, Nbulge, Nstars, Nbndry);
		#endif

		#ifdef VDEBUG
		printf("proc %i   debug 5.9\n",my_rank);
		#endif

		//loop over all the particles
		for(k=0; k<(int)particles_in.size(); k++){

			//compute in which cube the particle is
			//(the "if" is necessary to handle the case pos[i] == BoxSide, which results in m == mtot)
                        //(any other value of m/n/l not in [0, m/n/ltot-1] is an error)
                        m=(int) particles_in[k].pos[0]/lambda; if(m == mtot) m--;
                        n=(int) particles_in[k].pos[1]/lambda; if(n == ntot) n--;
                        l=(int) particles_in[k].pos[2]/lambda; if(l == ltot) l--;
			
			//if the cell should be ignored, skip to the next particle
			if (resolution_matrix[m][n][l].level == -1) continue;

			//if it does not belong to this proc, skip to the next particle
			if (resolution_matrix[m][n][l].cubeID < ps_min || resolution_matrix[m][n][l].cubeID >= ps_max) continue;

			//get the index (0-5) of the particle
			int particle_index = get_particle_index(k);

			#ifdef DEBUG
			if(m>=mtot || n>=ntot || l>=ltot || m<0 || n<0 || l<0) 
				printf("proc %i problem at k=%i pos=%f  %f  %f  l=%f  mtot=%i  ntot=%i  ltot=%i\n",my_rank,k,particles_in[k].pos[0],particles_in[k].pos[1],particles_in[k].pos[2],lambda,mtot,ntot,ltot);
			tot_part_check++;
			tot_mass_check1 += particles_in[k].mass;
			#endif

			//process the particle
			if(resolution_matrix[m][n][l].cubeID < 0){ //copy particles
				new_particles[particle_index].push_back(particles_in[k]);	//add a particle
				#ifdef DEBUG
				tot_mass_check2 += particles_in[k].mass;
				#endif
				#ifdef PRINT_IDS_CORRESPONDENCE
				old_ids[particle_index].push_back(vector< LOIinHigh >());
				old_ids[particle_index][old_ids[particle_index].size()-1].push_back(idH[k]);
				#endif
			} else if(resolution_matrix[m][n][l].cubeID >= 0){ //average over particles
				//compute the first (0 <-> 5) and second (0 <-> ps_max-ps_min+1) index for the CMs array,
				//first_index: put the highest-level particles in the first NSpecies slots, the
				//second-highest-levels ones in the second NSpecies slots and so on...
				int first_index = particle_index + (6 - species_number * (resolution_matrix[m][n][l].level + 1));
				int second_index = resolution_matrix[m][n][l].cubeID-ps_min;

				CMs[first_index][second_index].mass += particles_in[k].mass;
				//I put in pos (vel) the sum of pos(vel)*mass, and then divide by the total mass at the end
				CMs[first_index][second_index].pos[0] += particles_in[k].pos[0]*particles_in[k].mass;
				CMs[first_index][second_index].pos[1] += particles_in[k].pos[1]*particles_in[k].mass;
				CMs[first_index][second_index].pos[2] += particles_in[k].pos[2]*particles_in[k].mass;
				CMs[first_index][second_index].vel[0] += particles_in[k].vel[0]*particles_in[k].mass;
				CMs[first_index][second_index].vel[1] += particles_in[k].vel[1]*particles_in[k].mass;
				CMs[first_index][second_index].vel[2] += particles_in[k].vel[2]*particles_in[k].mass;
				#ifdef PRINT_IDS_CORRESPONDENCE
				IDs[first_index][second_index].push_back(idH[k]);
				#endif
			}
		}

		#ifdef VDEBUG
		printf("proc %i   debug 5.10\n",my_rank);
		#endif
 
		//free used memory
		particles_in.clear();
		idH.clear();
	}

	//free vectors and resolution_matrix (now useless)
	vector<particle_data>().swap(particles_in);
	vector<LOIinHigh>().swap(idH);

	for(m=0;m<mtot;m++){
		for(n=0;n<ntot;n++) delete[] resolution_matrix[m][n];
		delete[] resolution_matrix[m];
	}
	delete[] resolution_matrix;

	#ifdef VDEBUG
	printf("proc %i   debug 5.11\n",my_rank);
	#endif

	//now the last average for the merged particles
	for(i=0; i<6; i++){
		for(j=0;j<ps_max-ps_min;j++){
			if(CMs[i][j].mass > 0.0f){
				#ifdef DEBUG
				tot_mass_check2 += CMs[i][j].mass;
				#endif
				new_particles[i].push_back(particle_data(CMs[i][j].pos[0]/CMs[i][j].mass, CMs[i][j].pos[1]/CMs[i][j].mass, CMs[i][j].pos[2]/CMs[i][j].mass,
								CMs[i][j].vel[0]/CMs[i][j].mass, CMs[i][j].vel[1]/CMs[i][j].mass, CMs[i][j].vel[2]/CMs[i][j].mass,
								CMs[i][j].mass));
				#ifdef PRINT_IDS_CORRESPONDENCE
				old_ids[i].push_back(IDs[i][j]);
				#endif
			}
		}
	}

	#ifdef DEBUG
    printf("tot_part_check= %zu\n", tot_part_check);
	printf("tot_mass_check1 = %Lf\n", tot_mass_check1);
    printf("tot_mass_check2 = %Lf\n", tot_mass_check2);	
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank==0){
		printf("proc 0\n");
		for(i=0;i<new_particles.size();i++) for(j=0; j<new_particles[i].size(); j++) 
			printf("%f  %f  %f     %f\n",new_particles[i][j].pos[0],new_particles[i][j].pos[1],new_particles[i][j].pos[2],new_particles[i][j].mass);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank==1){	
		printf("proc 1\n");
		for(i=0;i<new_particles.size();i++) for(j=0; j<new_particles[i].size(); j++) 
			printf("%f  %f  %f     %f\n",new_particles[i][j].pos[0],new_particles[i][j].pos[1],new_particles[i][j].pos[2],new_particles[i][j].mass);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	#endif

	#ifdef VDEBUG
	printf("proc %i   debug 5.12\n",my_rank);
	#endif

	//free memory
	for(i=0; i<6-species_number; i++) delete[] CMs[i];
	delete[] CMs;
	#ifdef PRINT_IDS_CORRESPONDENCE
	IDs.clear();
	vector< vector< vector< LOIinHigh > > >.swap(IDs);
	#endif

	#ifdef VDEBUG
	printf("proc %i   debug 5.15\n",my_rank);
	#endif
}
