/*============================================================

* File Name : New_regions_finder.cpp

* Purpose : Used for zoom only. This function trace the particles in the
			highest-level region back to their position in the ICs, compute
			the new highest-level region as the smallest regions which contains
			all the highest-level particles and update the other-level regions
			mantaining the ratio between their radius and the highest-level one.

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include"global_variables.h"
#include"functions_declaration.h"

void New_regions_finder(){

	#ifdef VDEBUG
	printf("proc %i   debug 3.0\n",my_rank);
	#endif

	int i,k,h,N;

	//everything done by proc 0 and then broadcast results
	if(my_rank == 0){
		//boundaries of the high-res region in the ICs initially set to the whole box
		float x_minIC =  BoxSize+1;
		float x_maxIC = -BoxSize-1;
		float y_minIC =  BoxSize+1;
		float y_maxIC = -BoxSize-1;
		float z_minIC =  BoxSize+1;
		float z_maxIC = -BoxSize-1;
	

		//variables to store the particles:
		vector<LOIinHigh> high_res_ID;	// in the snapshot's high-res region
		vector<particle_data> rr;	// in the ICs' high-res region

		//here I read the SNAPSHOT to detect the particles inside the high-res region
		for(i=0;i<Sin_fnr;i++){
	
			//build up the file name
			build_file_name(snap_dir, snap_name, Sin_fnr, i, fname);

			#ifdef VDEBUG
			printf("proc %i   debug 3.1\n",my_rank);
			#endif
		
			print_info("New_regions_finder: processing snapshot file %s\n",fname);fflush(info);

			//actually read the file
			bool things_to_read[5] = { true, true, true, true, false }; /* header, pos, vel, id, mass */
			Read_Snap_File(fname, header1, particles_in, idL, false, things_to_read);
			N = Ngas + Nhalo + Ndisk + Nbulge + Nstars + Nbndry; //total number of particles

			#ifdef VDEBUG
			printf("proc %i   debug 3.7\n",my_rank);
			#endif

			//because of the shift of the particles, the center is in (BoxSize/2, BoxSize/2, BoxSize/2)
			Sphere high_res_region(BoxSize/2, BoxSize/2, BoxSize/2, level_bubbles_radii[levels_number-1]);

			//saving the ID of the particles inside the high-res region
			for(k=0; k<N; k++){
				if(high_res_region.contains(particles_in[k].pos)){
					high_res_ID.push_back(idL[k]);
				}
			}

			#ifdef VDEBUG
			printf("proc %i   debug 3.8\n",my_rank);
			#endif

			//free used memory
			particles_in.clear();
			idL.clear();
		}

		#ifdef DEBUG
		printf("high_res_ID.size()= %zu\n",high_res_ID.size());
		#endif

		#ifdef VDEBUG
		printf("proc %i   debug 3.9\n",my_rank);
		#endif

		//reading the ICs's files to find the high-res region in the ICs
		for(i=0; i<Lin_fnr;i++){

			build_file_name(Lic_dir, Lic_name, Lin_fnr, i, fname);

			#ifdef VDEBUG
			printf("proc %i   debug 3.10\n",my_rank);
			#endif

			print_info("New_regions_finder: processing ICs file %s\n",fname); fflush(info);
			
			//actually read file
			bool things_to_read[5] = { true, true, true, true, false }; /* header, pos, vel, id, mass */
			Read_Snap_File(fname, header1, particles_in, idL, false, things_to_read);
			N = Ngas + Nhalo + Ndisk + Nbulge + Nstars + Nbndry;

			#ifdef VDEBUG
			printf("proc %i   debug 3.16\n",my_rank);
			#endif

			//find the particles of this file in the high-res region
			for (k = 0; k<N; k++){	//loop over all ids
				for (h = 0; h<(int)high_res_ID.size(); h++){
					if (idL[k] == high_res_ID[h]){
						#ifdef DEBUG
						printf("SNAP particle located in ICs\n");
						#endif
						x_maxIC = max(x_maxIC, particles_in[k].pos[0]);
						x_minIC = min(x_minIC, particles_in[k].pos[0]);
						y_maxIC = max(y_maxIC, particles_in[k].pos[1]);
						y_minIC = min(y_minIC, particles_in[k].pos[1]);
						z_maxIC = max(z_maxIC, particles_in[k].pos[2]);
						z_minIC = min(z_minIC, particles_in[k].pos[2]);
						rr.push_back(particles_in[k]);
					}
				}
			}

	/*		//faster matching algorithm, to be tested (works only with sorted arrays)
	
			sort(high_res_ID.begin(),high_res_ID.end());
			sort(begin(id),end(id));

			for(k=0, h=0; k<N && h<high_res_ID.size(); ){
				if(id[k]==high_res_ID[h]){
					#ifdef DEBUG
					printf("SNAP particle located in ICs\n");
					#endif
					x_maxIC = max(x_maxIC, particles_in[k].pos[0]);
					x_minIC = min(x_minIC, particles_in[k].pos[0]);
					y_maxIC = max(y_maxIC, particles_in[k].pos[1]);
					y_minIC = min(y_minIC, particles_in[k].pos[1]);
					z_maxIC = max(z_maxIC, particles_in[k].pos[2]);
					z_minIC = min(z_minIC, particles_in[k].pos[2]);
					rr.push_back(particles_in[k]);
					k++;
					h++;
				}
				else if(id[k]<high_res_ID[h]){
					k++;
	    			}
				else{
					h++;
				}
			}
	*/		

			#ifdef VDEBUG
			printf("proc %i   debug 3.17\n",my_rank);
			#endif

			particles_in.clear();
			idL.clear();
		}

		high_res_ID.clear();
		vector<LOIinLow>().swap(idL);
		vector<LOIinHigh>().swap(high_res_ID);


		#ifdef VDEBUG
		printf("proc %i   debug 3.18\n",my_rank);
		#endif

		//center of the high-res region in the ICs
		c_IC[0]=0.5*(x_maxIC+x_minIC);
		c_IC[1]=0.5*(y_maxIC+y_minIC);
		c_IC[2]=0.5*(z_maxIC+z_minIC);
		//radius of the high-res region in the ICs
		r_high_IC=0;
		float nr;

		#ifdef VDEBUG
		printf("proc %i   debug 3.19\n",my_rank);
		#endif

		//compute new highest-level radius
		for(k=0;k<(int)rr.size();k++){
			nr=pow(rr[k].pos[0]-c_IC[0],2)+pow(rr[k].pos[1]-c_IC[1],2)+pow(rr[k].pos[2]-c_IC[2],2);
			r_high_IC = max(r_high_IC, nr);
		}

		r_high_IC = pow(r_high_IC, 0.5);

		//adjust other radii to maintain the ratio set in parameters file
		for(i=0; i<levels_number; i++){
			resolution_bubbles[i][0].setCenter(c_IC);
			resolution_bubbles[i][0].radius = r_high_IC * resolution_bubbles[i][0].radius / resolution_bubbles[levels_number-1][0].radius ;
		}

		rr.clear();
		vector<particle_data>().swap(rr);
	}

	//sending the computed quantities to the other process
    for(i=0; i<levels_number; i++) MPI_Bcast(&resolution_bubbles[i],1,MPI_SPHERE_t,0,MPI_COMM_WORLD);

	#ifdef VDEBUG
	printf("proc %i   debug 3.20\n",my_rank);
	#endif
}
