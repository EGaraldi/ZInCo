/*============================================================

* File Name : Slicer.cpp

* Purpose : Function to divide the simulation box in cubic cells using a resolution matrix
			and assign values depending on the resolution region they belong to.

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include"global_variables.h"
#include"functions_declaration.h"

void Slicer(){

	#ifdef VDEBUG
	printf("proc %i   debug 4.1\n",my_rank);
	#endif

	#ifdef DEBUG
	printf("trying to allocate: %i x %i x %i LOIout values\n",mtot, ntot, ltot);
	#endif

	//The resolution_matrix represents the cells in the box. It is used to determine the particle to merge. 
	//Each cubic cell is represented by two different fields in the matrix:
	// -level is the resolution level it belongs to
	// -cubeID is an unique
	resolution_matrix=new resolution_info**[mtot];
	for(int m=0;m<mtot;m++){
		resolution_matrix[m]=new resolution_info*[ntot];
		for(int n=0;n<ntot;n++) resolution_matrix[m][n]=new resolution_info[ltot];
	}

	#ifdef VDEBUG
	printf("proc %i   debug 4.2\n",my_rank);
	#endif
	
	//k=half the cube's diagonal
	float half_cube_diagonal = 0.5*1.732050808*lambda; // sqrt(3)
	float R;

	print_info("creating the resolution matrix\n");

	//now I loop over the cubic cells in the box, compute the distance of its center from the center of each resolution 
	//bubble. Then half of the cube's diagonal is subtracted in order to make sure that also cells with 
	//dist(centre_cell-centre_bubble) > R_bubble but dist(vertex_cell-centre_bubble) < R_bubble are assigned to the bubble.
	//This could lead to assign cell completely outside the bubble to the bubble itself, but this is a fair trade to avoid
	//a significantly higher complexity in checking for intersection between cell and bubble. Moreover, this can not reduce 
	//the high-res region but only increase it.
    #pragma omp parallel for collapse(3)
	for(int m=0;m<mtot;m++) for(int n=0;n<ntot;n++) for(int l=0;l<ltot;l++){
		//first set level to -1 to identify ignored cells
		resolution_matrix[m][n][l].level = -1;
		resolution_matrix[m][n][l].processed = false;
		//loop over levels
		for(int i=0; i<levels_number; i++){
			//loop over level's bubbles
			for(int j=0; j<level_size[i]; j++){
				//computing the distance
				R=pow(pow(resolution_bubbles[i][j].center[0]-(2*m+1)*0.5*lambda,2) + 
					pow(resolution_bubbles[i][j].center[1]-(2*n+1)*0.5*lambda,2) + 
					pow(resolution_bubbles[i][j].center[2]-(2*l+1)*0.5*lambda,2),0.5) - 
					half_cube_diagonal;
				//assign level if intersect. Since the level loop goes from 0 to LevelNumber, each cell 
				//is assigned to the higher level it intersect with
				if(resolution_bubbles[i][j].contains(R)) resolution_matrix[m][n][l].level = i;
			}
		}
	}

	#ifdef VDEBUG
	printf("proc %i   debug 4.3\n",my_rank);
	#endif

	//now set two counters to count for copied particles and domains to merge (rcounter counts with positive numbers, 
	//while hcounter uses negative numbers in order to allow to to use the counter to assign the domains ID and 
	//determine the action to take (merge or copy) just checking the sign of the domain ID.
	rcounter = 0;
	hcounter = -1;


	//Now assign an unique ID to each subdomain and count the particles to be (potentially) produced. 
	//To do so, loop again and try to merge adjacent domains.
	//{not easily OMP parallelizable since for a given (m,n,l) also neighbors are modified!}
	for (int m = 0; m<mtot; m++) for (int n = 0; n<ntot; n++) for (int l = 0; l<ltot; l++){
		if (level_cubic_region_side[resolution_matrix[m][n][l].level] > 0){ //merge particles
			//here we loop over the number of cubic cells to merge (the -1 is due to the fact that k is the number of
			//cells merged to the current one). We loop to make sure that if there are not enough adjacent cells to be 
			//merged, we try again with k reduced by one until we reach k=0 <-> the cell alone is merged.
			for (int k = level_cubic_region_side[resolution_matrix[m][n][l].level] - 1; k >= 0; k--){
				//now we loop on the k adjacent cells in all directions and check if all the adjacent cells share with 
				//the current one the same level (and have not yet been processed = merged). If not, go to the next cell
				for (int m1 = 0; m1 <= k; m1++) for (int n1 = 0; n1 <= k; n1++) for (int l1 = 0; l1 <= k; l1++){
					if ((m + m1 >= mtot) and(n + n1 >= ntot) and(l + l1 >= ltot)) goto check_failed; //check NOT out of bound
					if (resolution_matrix[m + m1][n + n1][l + l1].processed) goto check_failed; //check the neighbors are not yet processed
				}

				#ifdef VINFO
				if (my_rank == 0) print_info("Correspondence found! mid-res %i x %i x %i cubes in (%i,%i,%i). Marked!\n", k + 1, k + 1, k + 1, m, n, l);
				#endif

				//set cubeID and mark the cubes as already processed
				for (int m1 = 0; m1 <= k; m1++) for (int n1 = 0; n1 <= k; n1++) for (int l1 = 0; l1 <= k; l1++){
					resolution_matrix[m + m1][n + n1][l + l1].cubeID = rcounter;
					resolution_matrix[m + m1][n + n1][l + l1].processed = true;
				}

				//update rcounter
				rcounter++;

				check_failed:;
			}
		}
		else if (level_cubic_region_side[resolution_matrix[m][n][l].level] == 0){ //copy particles
			resolution_matrix[m][n][l].cubeID = hcounter;
			resolution_matrix[m][n][l].processed = true;

			//update hcounter
			hcounter--;
		} /* else if(level_cubic_region_side[i] < 0) particles to ignore => resolution_matrix[m][n][l].level = -1 but it is already set */
	}

	#ifdef DEBUG
	printf("rcounter=%i, hcounter=%i\n",rcounter,hcounter);
	#endif

	#ifdef VDEBUG
	printf("proc %i   debug 4.5\n",my_rank);
	#endif	
}
