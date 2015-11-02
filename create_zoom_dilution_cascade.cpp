/*============================================================

* File Name : create_zoom_dilution_cascade.cpp

* Purpose : These functions are used to start up a zoom, dilution or cascade ICs production,
			setting or randomizing the resolutions regions following the user's input.

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include"global_variables.h"
#include"functions_declaration.h"

//CASCADE: here after the functions used to build the resolution bubbles in
// a cascade run

//the first level is created separately from the others because we don't have a
//super-bubble but we have to deal with the entire simulation box
Sphere create_level_0_bubble(float radius){

	float x,y,z;

	//place the center randomly in the box
	x = BoxSize*float(rand())/float(RAND_MAX);
	y = BoxSize*float(rand())/float(RAND_MAX);
	z = BoxSize*float(rand())/float(RAND_MAX);

	return Sphere(x,y,z,radius);
}


//Now a function to create a new bubble inside an already-existing bubble.
Sphere create_sub_bubble(float old_center[], float old_radius, float new_radius){
	//random position in spherical coordinates
	float random1 = float(rand())/float(RAND_MAX);		/* 0 <-> 1 */
	float random2 = 180*float(rand())/float(RAND_MAX) - 90;	/* -90 <-> 90 */
	float random3 = 360*float(rand())/float(RAND_MAX);	/* 0 <-> 360 */

	//The new bubble MUST lie completely inside the existing one, so:
	float Dr = (random1)*(old_radius-new_radius);
	float theta = DegreesToRadiant(random2);
	float phi = DegreesToRadiant(random3);

	return Sphere(old_center[0]+Dr*sin(theta)*sin(phi), old_center[1]+Dr*sin(theta)*cos(phi), old_center[2]+Dr*cos(theta), new_radius);
}


//The following function will create the resolution bubbles of the cascade calling the above functions
//and checking for bubbles intersection (which must be avoided to allow smooth resolution transition)
void create_cascade(){

	//eveything will be done by the proc ranked 0, which will then send everything
	//to the other procs

	if (my_rank == 0){
		//set up random seed
		srand(cascade_random_seed);

		int i, j, k, h;

#ifdef DEBUG
		printf("creating level 0...\n");
#endif

		//write info in a seprate file
		FILE * cascade_info;
		sprintf(fname, "%s/Cascade_info.txt", output_dir);
		cascade_info = fopen(fname, "w");
		if (cascade_info == NULL){ display_info("ERROR while opening file %s\n", fname); /*error_flag=7;*/ }

		//first create bubbles for level0
		fprintf(cascade_info, "random seed used for the cascade: %i\n", cascade_random_seed);
		fprintf(cascade_info, "level 0 :\n");
		fprintf(cascade_info, "  (no super-bubble since this is the first level)\n");
		//I declare a vector which will contain all the bubble for a given level
		//Once all the bubbles are produced, the vector will be appended to the bubble 2D-vector
		vector<Sphere> sub;
		//set up a counter to check for the number of tries
		int bubble_tries = 0;
		for (j = 0; j < level_size[0];){
			Sphere new_sphere = create_level_0_bubble(level_bubbles_radii[0]);
			//check if it intersect the other spheres
			bool intersect = false;
			for (k = 0; k < j; k++) if (spheres_intersect(new_sphere, sub[k])) intersect = true;
			if (intersect){
				//try again or print advice for the maximum number of tries
				bubble_tries++;
				if (bubble_tries <= max_cascade_iter) continue; // = retry
				else{
					print_info("CRITICAL WARNING: Reached the maximum number of tries for placing bubbles on level 0\n");
					break;
				}
			}
			else{
				//add the sphere and continue
				sub.push_back(new_sphere);
				fprintf(cascade_info, "    bubble n. %i :\n", j);
				fprintf(cascade_info, "      center: %f %f %f :\n", new_sphere.center[0], new_sphere.center[1], new_sphere.center[2]);
				fprintf(cascade_info, "      radius: %f :\n", new_sphere.radius);
				j++;
			}
		}
		//once everything is done, we put all the bubbles of this level in the bubbles vector
		resolution_bubbles.push_back(sub);
		sub.clear();

#ifdef DEBUG
		printf("creating other levels...\n");
#endif
		//Now build the other levels: for each existing bubble, LevelSize bubble are created. They must lies completely inside the parent one
		for (i = 1; i < levels_number; i++){
#ifdef DEBUG
			printf("  level %i...\n", i);
#endif
			fprintf(cascade_info, "\nlevel %i :\n", i);
			//loop over the existing bubble of the previous level
			for (j = 0; j < (int)resolution_bubbles[i - 1].size(); j++){
				fprintf(cascade_info, "  super-bubble n. %i :\n", j);
				//loop over the current leve size
				for (k = 0; k < level_size[i];){
					Sphere new_sphere = create_sub_bubble(resolution_bubbles[i - 1][j].center, level_bubbles_radii[i - 1], level_bubbles_radii[i]);
					// check if it intersect the others
					bool intersect = false;
					for (h = 0; h < k; h++) if (spheres_intersect(new_sphere, sub[h])) intersect = true;
					if (intersect){
						//try again or print advice for the maximum number of tries
						bubble_tries++;
						if (bubble_tries <= max_cascade_iter) continue; // = retry
						else{
							print_info("CRITICAL WARNING: Reached the maximum number of tries for placing bubbles on level %i\n", i);
							break;
						}
					}
					else{
						//add the sphere and continue
						sub.push_back(new_sphere);
						fprintf(cascade_info, "    bubble n. %i :\n", k);
						fprintf(cascade_info, "      center: %f %f %f :\n", new_sphere.center[0], new_sphere.center[1], new_sphere.center[2]);
						fprintf(cascade_info, "      radius: %f :\n", new_sphere.radius);
						k++;
					}
				}
			}
			//once everything is done, we put all the bubbles of this level in the bubbles vector
			resolution_bubbles.push_back(sub);
			sub.clear();
		}
		fclose(cascade_info);
	}
}

//ZOOM: here after the function to initialize a zoom run
void create_zoom(){
	int i;
	for(i=0; i<levels_number; i++){ 
		//zoom must have 1 bubble per level, all with the same center
		level_size[i] = 1;
		//the resolution_bubbles vector still a 2D vector, so:
		vector<Sphere> sub;
		sub.push_back(Sphere(c, level_bubbles_radii[i]));
		resolution_bubbles.push_back(sub);
		sub.clear();
	}
}

//DILUTION: here after the function to initialize a dilution run
void create_dilution(){
	levels_number = 1;  //only one level of resolution
	level_size[0] = 1;  //with only one bubble
	level_cubic_region_side[0] = 1; //with merging parameter equal to 1
	//the resolution_bubbles vector still a 2D vector, so:
	vector<Sphere> sub;
	//the bubble must cover all the box, so:
	sub.push_back(Sphere(BoxSize/2.0f, BoxSize/2.0f, BoxSize/2.0f, BoxSize));
	resolution_bubbles.push_back(sub);
	sub.clear();
}
