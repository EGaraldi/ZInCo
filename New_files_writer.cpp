/*============================================================

* File Name : New_files_writer.cpp

* Purpose : This function handles the writing of new files. It is used to compute the
			particles to write in each file, to exchange the particles data between
			the procs and update the header with the new values. Then it calls
			Write_ICs_or_Snap_File to actually write them.

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include"global_variables.h"
#include"functions_declaration.h"

//This function handles the writing of new files.
void New_files_writer(){

	#ifdef VDEBUG
	printf("proc %i   debug 6.0\n",my_rank);
	#endif

	//some variables
	int i,j,k,next_process_to_contact=1,gimme_data=0;
	LOIout particles_available=0;
	long maxint=1; maxint = maxint << 32;	
	
	//struct to contain the header
	io_header wh;

	build_file_name(Hic_dir, Hic_name, Hin_fnr, 0, fname);

	#ifdef VDEBUG
	printf("proc %i   debug 6.1\n",my_rank);
	#endif

	print_info("New_files_writer: reading the first ICs file (%s) to get the header variables\n",fname);

	//read the existing header and then modify the variables changed by ZInCo
	bool things_to_read[5] = { true, false, false, false, false }; /* header, pos, vel, id, mass */
	Read_ICs_File(fname, wh, particles_in, idH, false, things_to_read);

	#ifdef DEBUG
	display_info("npart:  %i  %i  %i  %i  %i  %i\n",wh.npart[0],wh.npart[1],wh.npart[2], wh.npart[3],wh.npart[4],wh.npart[5]);
	display_info("massarr:  %f  %f  %f  %f  %f  %f\n",wh.massarr[0],wh.massarr[1],wh.massarr[2], wh.massarr[3],wh.massarr[4],wh.massarr[5]);
	display_info("time:  %f\n",wh.time);
	display_info("redshift:  %f\n",wh.redshift);
	display_info("flag_sfr:  %i\n",wh.flag_sfr);
	display_info("flag_feedback:  %i\n",wh.flag_feedback);
	display_info("npartTotal:  %i  %i  %i  %i  %i  %i\n",wh.npartTotal[0],wh.npartTotal[1],wh.npartTotal[2], wh.npartTotal[3],wh.npartTotal[4],wh.npartTotal[5]);

	printf("proc %i  new_particles dimensions=%zu  %zu  %zu  %zu  %zu  %zu\n",my_rank,new_particles[0].size(),new_particles[1].size(),new_particles[2].size(),
											new_particles[3].size(),new_particles[4].size(),new_particles[5].size());
	#endif
		
	#ifdef VDEBUG
	printf("proc %i   debug 6.2\n",my_rank);
	#endif

	//The file(s) will be written by proc 0, so first of all I need to collect the data from the other procs
	if(my_rank==0){
		//receiving data from other nodes
		LOIout **other_process_particles_number;	//number of particles of each process
		LOIout total_other_particles = 0;			//total number of particles received from other procs
		#ifdef DEBUG
		printf("trying to allocate: %i x 6 LOIout values\n",world_size-1);
		#endif
		other_process_particles_number=new LOIout*[world_size-1];
		for(i=0;i<world_size-1;i++) other_process_particles_number[i]=new LOIout[6];

		for(i=0;i<world_size-1;i++) for(j=0;j<6;j++){
			MPI_Recv(&other_process_particles_number[i][j],1,MPI_LOIout,i+1,mpi_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			total_other_particles += other_process_particles_number[i][j];
		}

		//Now compute how many particles per file should be written
		particles_per_file = total_other_particles;
		for(i=0;i<6;i++) particles_per_file += new_particles[i].size();
		if(particles_per_file == 0){
			printf("ERROR: no particles to be produced! Exiting...\n");
			error_flag = 11;
		}
		else particles_per_file /= out_fnr;

		error_flag_check();
		#ifdef DEBUG
                printf("particles_per_file %i\n",particles_per_file);
		printf("time:  %f\n",wh.time);
		printf("redshift:  %f\n",wh.redshift);
		#endif

		//update header variables 
		long int long_total_npart[6];
		for(i=0;i<6;i++){
			long_total_npart[i] = new_particles[i].size();
			for(j=0;j<world_size-1;j++) {
				long_total_npart[i] += other_process_particles_number[j][i];
			}
		}
		for(i=0;i<6;i++){
			//massarr is 0 everywhere because the mass is always written (because it is different for each particles)
			//DONE?: easy improvement: if the particles COPIED have the same mass => write the mass in massarr and not in the mass block
			int species = i%species_number;
			int level = (species + 6 - i) / species_number + 1;
			if (level_cubic_region_side[level] != 0 || wh.massarr[species] == 0)
				wh.massarr[i] = 0.0f;
			wh.npartTotalHW[i] = int(long_total_npart[i] / maxint);
			wh.npartTotal[i] = int(long_total_npart[i] % maxint);
		}

		//now an option to write in a separate file the correspondence between old and new IDs
		//NOTE: Not really tested yet
		#ifdef PRINT_IDS_CORRESPONDENCE
		sprintf(fname,"%s/IDs_correspondence.txt", output_dir);
		id_corr = fopen(fname,"w");
		if (id_corr == NULL){ display_info("ERROR while opening file %s\n", fname); error_flag=7; }
		error_flag_check();
		#endif

		//Now we actually write the file(s)
		int filenum = 0; //variable to keep track of the number of file written
		for(filenum=0; filenum<out_fnr; filenum++){

			//check if proc 0 has sufficient particles for the current file
			particles_available = 0;
			for(i=0; i<6; i++) particles_available += new_particles[i].size();
			#ifdef DEBUG
            printf("particles_available %i\n",particles_available);
			#endif

			//if it's not the case, ask for data and collect them
			while((particles_available < particles_per_file) && (next_process_to_contact < world_size)){
				#ifdef DEBUG
				printf("requesting particles to proc %i\n",next_process_to_contact);
				#endif
				//first tell the next proc that data are needed
				MPI_Send(&gimme_data,1,MPI_INT,next_process_to_contact,mpi_tag,MPI_COMM_WORLD);
				print_info("New_files_writer: Requesting particles to proc %i...",next_process_to_contact); fflush(info);

				//then receive the data sent
				for(j=0; j<6; j++) for(k=0;k<other_process_particles_number[next_process_to_contact-1][j];k++){
					particle_data next_particle;
					MPI_Recv(&next_particle,1,MPI_PARTICLE_t,next_process_to_contact,mpi_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
					new_particles[j].push_back(next_particle);
				}

				#ifdef PRINT_IDS_CORRESPONDENCE
				LOIinHigh received_id;
				for(j=0; j<6; j++) for(k=0;k<other_process_particles_number[next_process_to_contact-1][j];k++){
					old_ids[j].push_back(vector<LOIinHigh>());
					int old_ids_size;
					MPI_Recv(&old_ids_size, 1, MPI_INT, next_process_to_contact, mpi_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for(i=0; i<old_ids_size; i++){
						MPI_Recv(&received_id, 1, MPI_LOIinHigh, next_process_to_contact, mpi_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						old_ids[j][k].push_back(received_id);
					}
				}
				#endif

				//switch to the next procs
				next_process_to_contact++;
			
				print_info("recieved!\n"); fflush(info);

				#ifdef DEBUG
				printf("proc %i   newp.size()=%zu   newp:\n",my_rank,newp.size());
				for(i=0;i<newp.size();i++) printf("%f  %f  %f  %f\n",newp[i].pos[0],newp[i].pos[1],newp[i].pos[2],newp[i].mass);
				#endif
			}
	
			//Once the main proc has sufficient particles, it starts writing the particles...

			#ifdef VDEBUG
			printf("proc %i   debug 6.3\n",my_rank);
			#endif

			#ifdef DEBUG
			printf("newp sizes: %i  %i  %i  %i  %i  %i\n",new_particles[0].size(),new_particles[1].size(),new_particles[2].size(),
									new_particles[3].size(),new_particles[4].size(),new_particles[5].size());
			#endif
	
			//determine the particles to write in this file and the correspondent npart
			int left = particles_per_file;	//particles left to write in the current file

			if(filenum != out_fnr-1){
				//here we parse all the particles types and add the available number of particle
				//to the 'list' of particle to write = header.npart and consequentely decrease left
				for (j = 0; j < 6; j++){
					if ((int)new_particles[j].size() < left){
						wh.npart[j] = new_particles[j].size();
						left -= new_particles[j].size();
					} else{
						wh.npart[j] = left;
						break; //because now left = 0
					}
				}
			} else{
				//if the file is the last one, I write ALL the particles
				for (j = 0; j<6; j++) wh.npart[j] = new_particles[j].size();
			}

			#ifdef DEBUG
			printf("wh_npart: %i  %i  %i  %i  %i  %i\n",wh.npart[0],wh.npart[1],wh.npart[2],wh.npart[3],wh.npart[4],wh.npart[5]);
			#endif

			//update the header with the number of files in output
			wh.num_files = out_fnr;

			#ifdef VDEBUG
			printf("proc %i   debug 6.4\n",my_rank);
			#endif

			//Build file name
			char f[200];
			if(dilution) sprintf(f,"%s/%s.dilution%i",output_dir,Hic_name,cubes_per_side);
			else if(zoom) sprintf(f,"%s/%s.zoom",output_dir,Hic_name);
			else if(cascade) sprintf(f,"%s/%s.cascade",output_dir,Hic_name);
			if(out_fnr > 1) sprintf(fname,"%s.%i",f,filenum); else sprintf(fname,"%s",f);

			//Call the function to actually write data
			print_info("New_files_writer: writing on file %s\n",fname);
			bool things_to_write[5] = { true, true, true, true, true }; /* header, pos, vel, id, mass */
			Write_ICs_or_Snap_File(fname, wh, new_particles, things_to_write);

			//Remove the particles written
			if(out_fnr == 1) new_particles.clear();
			else for(j=0; j<6; j++) if(wh.npart[j] > 0) new_particles[j].erase(new_particles[j].end()-wh.npart[j], new_particles[j].end());

			#ifdef PRINT_IDS_CORRESPONDENCE
			if(out_fnr == 1) old_ids.clear();
			else for(j=0; j<6; j++) if(wh.npart[j] > 0) old_ids[j].erase(old_ids[j].end()-wh.npart[j], old_ids[j].end());
			#endif

			#ifdef VDEBUG
			printf("proc %i   debug 6.10\n",my_rank);
			#endif
		}

		#ifdef PRINT_IDS_CORRESPONDENCE
		fclose(id_corr);
		#endif

		//delete array now useless
		for(i=0;i<world_size-1;i++) delete[] other_process_particles_number[i];
		delete[] other_process_particles_number;

	} else{
		//send particles available
		for(i=0;i<6;i++){
			LOIout new_particle_size = new_particles[i].size();
			MPI_Send(&new_particle_size,1,MPI_LOIout,0,mpi_tag,MPI_COMM_WORLD);
		}

		error_flag_check();

		#ifdef PRINT_IDS_CORRESPONDENCE
		error_flag_check();
		#endif

		//receive the gimme_data variables: this always happen, so the proc waits until it receives the gimme_data variable and then send its particles
		MPI_Recv(&gimme_data,1,MPI_INT,0,mpi_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		//send data when requested
		for(j=0; j<6; j++) for(k=0;k<(int)new_particles[j].size();k++) MPI_Send(&new_particles[j][k],1,MPI_PARTICLE_t,0,mpi_tag,MPI_COMM_WORLD);

		#ifdef PRINT_IDS_CORRESPONDENCE
		for(j=0; j<6; j++) for(k=0;k<(int)new_particles[j].size();k++){
			int size_to_send = old_ids[j][k].size();
			MPI_Send(&size_to_send, 1, MPI_INT, 0, mpi_tag, MPI_COMM_WORLD);
			for(i=0; i<(int)old_ids[j][k].size(); i++) MPI_Send(&old_ids[j][k][i], 1, MPI_LOIinHigh, 0, mpi_tag, MPI_COMM_WORLD);
		}
		#endif

		error_flag_check(); //match the ones in Write_ICs_or_Snap_file called by my_rank==0
	}
}
