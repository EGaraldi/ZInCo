/*============================================================

* File Name : Read_parameter_file.cpp

* Purpose : Function to read the parameter file and assign the read 
			values to the corresponding variables 
			(adjusted from the one used in GADGET-2 by 
			[V. Springel, MNRAS, 364 (2005) 1105–1134])

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/

#include"global_variables.h"
#include"functions_declaration.h"

void Read_parameter_file(char *ParameterFile){
      
	#define REAL 1
	#define STRING 2
	#define INT 3
	#define MAXTAGS 300

	FILE *fd, *fdout;
	char buf[200], buf1[200], buf2[200], buf3[400];
	int i, j, nt;
	int id[MAXTAGS];
	void *addr[MAXTAGS];
	char tag[MAXTAGS][50];
	int  errorFlag = 0;


	if(sizeof(long long) != 8){
		display_info("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
		my_exit(5);
	}

	if(sizeof(int) != 4){
		display_info("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
		my_exit(5);
	}

	if(sizeof(float) != 4){
		display_info("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
		my_exit(5);
	}

	if(sizeof(double) != 8){
		display_info("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
		my_exit(5);
	}

	if(my_rank == 0){		/* read parameter file on process 0 */
		nt = 0;

/*		strcpy(tag[nt], "Dilution");
		addr[nt] = &dilution;
		id[nt++] = INT;

		strcpy(tag[nt], "Zoom");
		addr[nt] = &zoom;
		id[nt++] = INT;

		strcpy(tag[nt], "Cascade");
		addr[nt] = &cascade;
		id[nt++] = INT;
*/	
		strcpy(tag[nt], "RunType");
		addr[nt] = &run_type;
		id[nt++] = STRING;

		strcpy(tag[nt], "HighResICsDir");
		addr[nt] = &Hic_dir;
		id[nt++] = STRING;

		strcpy(tag[nt], "HighResICsName");
		addr[nt] = &Hic_name;
		id[nt++] = STRING;

		strcpy(tag[nt], "HighResICsFilesNumber");
		addr[nt] = &Hin_fnr;
		id[nt++] = INT;

		strcpy(tag[nt], "OutputDir");
		addr[nt] = &output_dir;
		id[nt++] = STRING;

		strcpy(tag[nt], "OutputFilesNumber");
		addr[nt] = &out_fnr;
		id[nt++] = INT;

		strcpy(tag[nt], "CubesPerSide");
		addr[nt] = &cubes_per_side;
		id[nt++] = INT;

		strcpy(tag[nt], "LevelsNumber");
		addr[nt] = &levels_number;
		id[nt++] = INT;

		strcpy(tag[nt], "SpeciesNumber");
		addr[nt] = &species_number;
		id[nt++] = INT;

		strcpy(tag[nt], "Level0Size");
		addr[nt] = &level_size[0];
		id[nt++] = INT;

		strcpy(tag[nt], "Level1Size");
		addr[nt] = &level_size[1];
		id[nt++] = INT;

		strcpy(tag[nt], "Level2Size");
		addr[nt] = &level_size[2];
		id[nt++] = INT;

		strcpy(tag[nt], "Level3Size");
		addr[nt] = &level_size[3];
		id[nt++] = INT;

		strcpy(tag[nt], "Level4Size");
		addr[nt] = &level_size[4];
		id[nt++] = INT;

		strcpy(tag[nt], "Level5Size");
		addr[nt] = &level_size[5];
		id[nt++] = INT;

		strcpy(tag[nt], "Level0BubbleRadius");
		addr[nt] = &level_bubbles_radii[0];
		id[nt++] = REAL;

		strcpy(tag[nt], "Level1BubbleRadius");
		addr[nt] = &level_bubbles_radii[1];
		id[nt++] = REAL;

		strcpy(tag[nt], "Level2BubbleRadius");
		addr[nt] = &level_bubbles_radii[2];
		id[nt++] = REAL;

		strcpy(tag[nt], "Level3BubbleRadius");
		addr[nt] = &level_bubbles_radii[3];
		id[nt++] = REAL;

		strcpy(tag[nt], "Level4BubbleRadius");
		addr[nt] = &level_bubbles_radii[4];
		id[nt++] = REAL;

		strcpy(tag[nt], "Level5BubbleRadius");
		addr[nt] = &level_bubbles_radii[5];
		id[nt++] = REAL;

		strcpy(tag[nt], "Level0CubesPerSide");
		addr[nt] = &level_cubic_region_side[0];
		id[nt++] = INT;

		strcpy(tag[nt], "Level1CubesPerSide");
		addr[nt] = &level_cubic_region_side[1];
		id[nt++] = INT;

		strcpy(tag[nt], "Level2CubesPerSide");
		addr[nt] = &level_cubic_region_side[2];
		id[nt++] = INT;

		strcpy(tag[nt], "Level3CubesPerSide");
		addr[nt] = &level_cubic_region_side[3];
		id[nt++] = INT;

		strcpy(tag[nt], "Level4CubesPerSide");
		addr[nt] = &level_cubic_region_side[4];
		id[nt++] = INT;

		strcpy(tag[nt], "Level5CubesPerSide");
		addr[nt] = &level_cubic_region_side[5];
		id[nt++] = INT;

		strcpy(tag[nt], "LowResICsDir");
		addr[nt] = &Lic_dir;
		id[nt++] = STRING;

		strcpy(tag[nt], "LowResICsName");
		addr[nt] = &Lic_name;
		id[nt++] = STRING;

		strcpy(tag[nt], "LowResICsFilesNumber");
		addr[nt] = &Lin_fnr;
		id[nt++] = INT;

		strcpy(tag[nt], "SnapshotDir");
		addr[nt] = &snap_dir;
		id[nt++] = STRING;

		strcpy(tag[nt], "SnapshotName");
		addr[nt] = &snap_name;
		id[nt++] = STRING;

		strcpy(tag[nt], "SnapshotFilesNumber");
		addr[nt] = &Sin_fnr;
		id[nt++] = INT;

		strcpy(tag[nt], "x_c");
		addr[nt] = &c[0];
		id[nt++] = REAL;

		strcpy(tag[nt], "y_c");
		addr[nt] = &c[1];
		id[nt++] = REAL;

		strcpy(tag[nt], "z_c");
		addr[nt] = &c[2];
		id[nt++] = REAL;

		strcpy(tag[nt], "CascadeRandomSeed");
		addr[nt] = &cascade_random_seed;
		id[nt++] = INT;

		strcpy(tag[nt], "CascadeMaxIter");
		addr[nt] = &max_cascade_iter;
		id[nt++] = INT;

		if((fd = fopen(ParameterFile, "r"))){
			sprintf(buf, "%s%s", ParameterFile, "-usedvalues");
			if(!(fdout = fopen(buf, "w"))){
				printf("error opening file '%s' \n", buf);
				errorFlag = 1;
			} else{
				while(!feof(fd)){
					//char *ret;
					*buf = 0;
					/*ret = */fgets(buf, 200, fd);
					if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2) continue;
					if(buf1[0] == '%') continue;

					for(i = 0, j = -1; i < nt; i++) if(strcmp(buf1, tag[i]) == 0){
						j = i;
						tag[i][0] = 0;
						break;
					}

					if(j >= 0){
						switch (id[j]){
							case REAL:
								*((float *) addr[j]) = atof(buf2);
								fprintf(fdout, "%-35s%g\n", buf1, *((float *) addr[j]));
								break;
							case STRING:
								strcpy((char *) addr[j], buf2);
								fprintf(fdout, "%-35s%s\n", buf1, buf2);
								break;
							case INT:
								*((int *) addr[j]) = atoi(buf2);
								fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
								break;
						}
					} else{
						fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",ParameterFile, buf1);
						errorFlag = 1;
					}
				}
				fclose(fd);
				fclose(fdout);
			}
		} else{
			display_info("Parameter file %s not found.\n", ParameterFile);
			errorFlag = 1;
		}

		for(i = 0; i < nt; i++){
			if(*tag[i]){
				display_info("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], ParameterFile);
				errorFlag = 1;
			}
		}
	}

	MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(errorFlag) my_exit(9);

/*	MPI_Bcast(&dilution,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&zoom,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&cascade,1,MPI_INT,0,MPI_COMM_WORLD);
*/
	MPI_Bcast(&run_type,sizeof(run_type)/sizeof(char),MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&Hic_dir,sizeof(Hic_dir)/sizeof(char),MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&Hic_name,sizeof(Hic_name)/sizeof(char),MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&Hin_fnr,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&output_dir,sizeof(output_dir)/sizeof(char),MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&out_fnr,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&cubes_per_side,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&levels_number,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&species_number,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&level_size,6,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&level_bubbles_radii,6,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&level_cubic_region_side,6,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lic_dir,sizeof(Lic_dir)/sizeof(char),MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lic_name,sizeof(Lic_name)/sizeof(char),MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&Lin_fnr,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&snap_dir,sizeof(snap_dir)/sizeof(char),MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&snap_name,sizeof(snap_name)/sizeof(char),MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&Sin_fnr,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&c,3,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Bcast(&cascade_random_seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&max_cascade_iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	#undef REAL
	#undef STRING
	#undef INT
	#undef MAXTAGS
}
