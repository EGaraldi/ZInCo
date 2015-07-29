/*============================================================

* File Name : functions_declaration.h

* Purpose : functions prototypes

* Created By :  Enrico Garaldi
(PhD student,
Argelander Institut fur Astronomie,
Bonn University,
Bonn, Germany)

=============================================================*/


//create_zoom_dilution_cascade.cpp
Sphere create_level_0_bubble(float);
Sphere create_sub_bubble(float[], float, float);
void create_cascade();
void create_zoom();
void create_dilution();

//initialize_check_finalize_variables.cpp
void check_ignored_particles();
void initialize(int);
void finalize();
void check_variables();

//New_files_writer.cpp
void New_files_writer();

//New_particles_maker.cpp
void New_particles_maker();

//New_regions_finder.cpp
void New_regions_finder();

//other_small_functions.cpp
size_t my_fread(void *, size_t, size_t, FILE *);
size_t my_fwrite(void *, size_t, size_t, FILE *);
void my_exit(int);
void error_flag_check();
void print_info(const char*, ...);
void display_info(const char*, ...);
int get_particle_index(int);
void build_file_name(char[], char[], int, int, char[]);
bool spheres_intersect(Sphere, Sphere);
float DegreesToRadiant(int);

//Read_ICs_or_Snap_File.cpp
void Read_ICs_File(char[], io_header&, vector<particle_data>&, vector<LOIinHigh>&, bool, bool[]);
void Read_Snap_File(char[], io_header&, vector<particle_data>&, vector<LOIinLow>&, bool, bool[]);

//Read_parameter_file.cpp
void Read_parameter_file(char *);

//Slicer.cpp
void Slicer();

//Write_ICs_or_Snap_File.cpp
void Write_ICs_or_Snap_File(char[], io_header&, vector< vector<particle_data> >&, bool[]);