/* ************************************************************************
PatchMapper.c
A tool to identify and map contiguous patches in a raster image, and perform
	certain other analyses on the patches in that image
TO DO:
1. get rid of parameter m for missing. Missing is hardwired as zero in some places, and it doesn't
matter because re-coding is possible, so there's no reason to have a parameter for missing values,
only a requirement that missing MUST be zero, by recode if necessary.
2. check into parallelizing patch definer.
3. add file writes for basic pixel and patch stats that are now spewed to the terminal.


Compilation: 
gcc -std=c99 -m64 -Wall -fopenmp PatchMapper.c -o patchmapper  -O2 -lm

Usage notes
Missing values:
	should be zero on the input dataset, or recoded to zero via parameters.recode=1
	after any re-coding, will be set to 'internal' code 255 in {main}.
	
Required files:
	size.txt
	inputmap
	parfile.txt
Optional files:
	recode.txt	If parameter Z = 1; a lookup table for re-coding input pixel values. One line per change. Format <old code><space><new code>
Output files:
	If parameter X = 1; patchmap - an image containing sequential patch numbers from 1...to N total patches
	If parameter P = 1; patchstatistics.csv - a text file with info about each patch:
		If parameter N = 4 (4-neighbor)
			Patch_number,Patch_type,Num_pixels,Out_edges,In_edges,Out_pixels,Row_origin,Col_origin,Touch_border,Row_max,Col_min,Colmax
		If parameter N = 8 (8-neighbor)
			Patch_number,Patch_type,Num_pixels,Row_origin,Col_origin,Row_max,Col_min,Colmax
		
Parameter file format, each line has a letter (case insensitive) and a paramter value separated by spaces
m	0	pixel value representing missing value (after re-coding if parameter Z = 1; default=0)
z	0	re-code input pixels? (0=no; 1=yes; default = 0)
n	4	neighbor rule (either 4 or 8; default = 4)
x	0	output patch map? (0=no; 1=yes; default = 0)
p   0	calculate and export patch-level statistics? (0=no, 1=yes)
	
Misc notes.
	Here is the magic to access a 1-D array of size rows*columns using [row][col] syntax.
		Example. for unsigned char *mat_in of size rxc, and want to access as Matrix_1[r][c], declare
		unsigned char (*Matrix_1)[ncols] = (unsigned char (*)[ncols])mat_in;
		This casting requires the same typedef on both sides
		Then Matrix_1[r][c] is the same as *(mat_in + (r*ncols) + c)
		Reference: https://stackoverflow.com/questions/32317605/treating-one-dimensional-array-as-two-dimensional-at-run-time

Changelog
Version 0.91  February 2024
Added the accumulation and text export of individual patch statistics (size, class, location of patch origin)
Version 0.9 October 2021
Began re-construction of key capabilities from LandStat.c (ca 1991-1998), a program
	that was used to calculate the landscape pattern metrics described in the appendix of
	Riitters et al 1995 (https://doi.org/10.1007/BF00158551). As implemented here, the metrics
	may differ from the published versions. Important code changes include:
	1. implementation of 8-neighbor connectivity in addition to 4-neighbor.
	2. certains aspects parallelized
	This version sets up a few basic routines (I/O, parameter handling, patch identification, re-coding input pixels,
		some basic whole-image statistics). As before, patch identification is achieved using a linked-list, flood-
		filling algorithm.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include  <omp.h>
#define max(a,b) (((a) > (b)) ? (a) : (b)) 
#define min(a,b) (((a) < (b)) ? (a) : (b))

// prototypes 
void Patch_Definer(unsigned int, unsigned int); 
void Find_Neighbors(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int); 
void Add_Queue(unsigned int, unsigned int);
int Read_Parameter_File(FILE *);
int Pixel_Analyzer(unsigned int, unsigned int) ;
void Patch_Analyzer_4neighbor(unsigned int, unsigned int); 
void Patch_Analyzer_8neighbor(unsigned int, unsigned int); 

// make some things visible everywhere 
unsigned char *mat_in;
unsigned int *mat_patnum;
unsigned char *mat_util;
unsigned int *edge_matrix;
// run or control parameters
struct options  {              	// Holds run-time requests or default values  
	 int Missing_input;        	// 'Missing' data code value used in the input data matrix 
	 int Rows;                  // number of rows and columns as found in the data set 
	 int Cols;
}; 
struct options control =
   {0, 0, 0};    /* Defaults */ 

struct queue_1 {					// holds a linked-list queue 
    unsigned int Row;			// the row and column of an entry define a pixel in the queue 
	unsigned int Col;
	struct queue_1 *next_ptr;		//pointer to next element in the queue  
}; 
struct queue_1 *first_ptr;			// pointer to first element in the queue, Add_Queue expects this be declared in common  

struct run_parameters {
	int missing_value_code;   	//after recode, if any 
	int recode;       		  	// 0 = no recode, 1 = recode, supply recode values in recode file
	int neighbor_rule;			// either 4 or 8 
	int write_patch_map;		// 0=no, 1=yes 
	int patch_analyzer;			// 0=no, 1=yes
};                              // Default parameter values 
struct run_parameters parameters = {0,0,4,0,0};

// overall image statistics
 struct statistics_1 {
    unsigned int Num_pixels; 				// Number of pixels with non-missing data  
	int Num_types;             					// Number of non-missing types in matrix  
	float Pix_percents [256];          			// Percent of pixels in matrix by type  
	unsigned int Pix_counts [256]; 		// Number of pixels in matrix by type  
	unsigned int Num_patches;   			// Number of patches of all sizes in the matrix  
	unsigned int Num_type_patches [256]; 	// Number of patches of all sizes by type.  
	float Average_patch_size [256];
	float Average_patch_size_overall;
};     
struct statistics_1 matrix_stats =  
     { 0,0,{0.},{0},0,{0},{0.},0.}; 

// patch statistics, The patches are numbered from 1...N total patches. A reduced set is available for 8-neighbor
struct statistics_2 {             
    unsigned char Patch_type; 	// The internal type of this patch, 0 - 254 
	unsigned int Num_pixels;    // Number of pixels in patch       
	unsigned int Out_edges;     // Number of outside edges, ie. perimeter 
	unsigned int In_edges;      // Number of inside edges                   
	unsigned int Out_pixels;    // Number of outside pixels (# to enclose)  
	unsigned int Row_org;       // Starting row of this patch   Note that cell Row_org, Col_org is guaranteed to be in a given patch   
	unsigned int Col_org;       // Starting column of this patch   
	unsigned char Touch_border;	// 1 if patch touches border OR a missing value. This is called "touch_flag" on csv output
	unsigned int Row_max;       // Max row where this patch occurs. (Row_min is == Row_org) 
	unsigned int Col_min;       // Min column where this patch occurs. 
	unsigned int Col_max;       // Max column where this patch occurs.
}; 
struct statistics_2 *patch_stats;   

int main(int argc, char **argv)
{
	static int row, col, nrows, ncols, index, ret_val, vold, vnew;
	static int dumint, ind1, counter, recode_table[256], temp_int;
	static int Out_to_in[256];      // Array(x) = internal code for external code x 
	static int In_to_out[256];      // Array(x) = external code for internal code x 
	char filename_in[100], filename_siz[100], filename_par[100], 
		 filename_rec[100], filename_pat[100],header_line[30], filename_patstats[100];
	FILE *infile, *sizfile, *parfile, *recfile, *patfile, *patstatsfile;
//			/* as needed ... windows and linux are different
			printf("\nMake sure an int is 4 bytes");
			index = sizeof(char);
			printf("\n char is %d bytes\n", index);
			index = sizeof(short int);
			printf("\n short int is %d bytes\n", index);
			index = sizeof(int);
			printf("\n int is %d bytes\n", index);
			index = sizeof(long int);
			printf("\n long int is %d bytes\n", index);
			index = sizeof(long long int);
			printf("\n long long int is %d bytes\n", index);
			index = sizeof(float);
			printf("\n float is %d bytes\n", index);
//			*/
	printf("\nLoaded and running");
	strcpy(filename_siz, "size.txt");
	strcpy(filename_in, "inputmap");
	strcpy(filename_rec, "recode.txt");
	strcpy(filename_par, "parfile.txt");
	strcpy(filename_pat, "patchmap");
	strcpy(filename_patstats, "patchstatistics.csv");
	setbuf(stdout, NULL);
	// Open some files 
	if( (infile = fopen(filename_in, "rb") ) == NULL){printf("\nError opening %s\n", filename_in);exit(15);}
	if( (sizfile= fopen(filename_siz, "r") ) == NULL){printf("\nError opening %s\n", filename_siz);fclose(infile);exit(17);}
	if( (parfile= fopen(filename_par, "r") ) == NULL){printf("\nError opening %s\n", filename_par);fclose(infile);fclose(sizfile);exit(7);}
	if( (ret_val = Read_Parameter_File(parfile)) != 0){
		printf("\nError reading parameter file");
		fclose(infile);fclose(sizfile);exit(19);
	}
	fclose(parfile);
	printf("\nParameters: m %d  z %d  n %d  x %d p %d", parameters.missing_value_code, parameters.recode, parameters.neighbor_rule,
		parameters.write_patch_map,  parameters.patch_analyzer);
	// Read the siz file and save nrows and ncols in control structure
	if(fscanf(sizfile,"%s %d", header_line, &nrows) != 2) {printf("\nError reading size file"); exit(18);}
	if(fscanf(sizfile,"%s %d", header_line, &ncols) != 2) {printf("\nError reading size file"); exit(18);}
	fclose(sizfile);
	control.Rows = nrows;
	control.Cols = ncols;
	printf("\nAllocating memory for input data, patchmap, and utility matrix");
	// mat_in is the input map 
	if( (mat_in = (unsigned char *)calloc( (nrows * ncols), sizeof(unsigned char) ) ) == NULL ) {printf("\nMalloc failed, mat_in.\n");exit(20);}
	// mat_patnum is the patch number map
	if( (mat_patnum = (unsigned int *)calloc( (nrows * ncols), sizeof(unsigned int) ) ) == NULL ){printf("\nMalloc failed, mat_patnum.\n");exit(21);}
	// mat_util is a utility matrix, note it's buffered by one pixel all around 
	if( (mat_util = (unsigned char *)calloc( ((nrows+2) * (ncols+2)), sizeof(unsigned char) ) ) == NULL ) {printf("\nMalloc failed, mat_util.\n");exit(22);}
	// read the input data
	printf("\nReading %d columns and %d rows from file %s.", ncols, nrows, filename_in);
	if(fread(mat_in, 1, (nrows * ncols), infile) != (nrows * ncols) ) {printf("\nError reading input map.\n"); free(mat_in); exit(23);}
	fclose(infile);
	// re-code if requested 
	// If recoding, open and process the recode file 
	if(parameters.recode == 1) {
		if( (recfile = fopen(filename_rec, "r") ) == NULL) {
			printf("\nError opening %s\n", filename_rec); exit(5);
		}	
		for(index = 0; index < 256; index ++) {
			recode_table[index] = index; 
		}
		printf("\nRe-coding input map");
		printf("\n   Old code ---> New code\n");
		while(fscanf(recfile, "%d %d", &vold, &vnew) == 2) {
			recode_table[vold] = vnew;
			printf("   %8d ---> %3d\n", vold, vnew);
		}
		fclose(recfile);
// OMP the recode		
#pragma omp parallel  for  	 private ( row, index)		
		for (row=0; row < nrows; row++) {
			index = row * ncols;
#pragma omp parallel  for  	 private ( col, dumint)
			for (col=0; col < ncols; col++) {
				dumint = *(mat_in + index + col);
				*(mat_in + index + col) = recode_table[dumint];
			}
		}
	}
	// Recode input so that type codes start at zero and go 
	//  to t-1 types, with missing values set to type code 255 
	// Initialize the lookup tables  
	for(ind1 = 0; ind1 < 256; ind1++) { 
		Out_to_in [ind1] = -9;
		In_to_out [ind1] = -9;
	} 
	Out_to_in [0] = 255;
	In_to_out [255] = 0;
	counter = -1;
	for(row = 0; row < nrows; row++) {
		index = row * ncols;
		for(col=0; col < ncols; col++) {
			dumint = *(mat_in + index + col);
			if(Out_to_in[dumint] == -9) {
					counter++;  // assign a new non-missing code 
					Out_to_in[dumint] = counter; 
					In_to_out[counter] = dumint;
				
			}
			// store converted type code as a character 
			*(mat_in + index + col) = Out_to_in[dumint];
		} 
	} 
	// Set the number of non-missing types found in the matrix 
	matrix_stats.Num_types = counter + 1;
	printf("\nNumber of non-missing types: %d", matrix_stats.Num_types);
	printf("\nInternal and External codes:");
	for(ind1 = 0; ind1 < matrix_stats.Num_types; ind1++) {
			printf("\n%4d %4d", ind1, In_to_out[ind1]);
	}
	printf("\n 255 %4d", In_to_out[255]);
	//basic pixel statistics 
	printf("\nGetting basic pixel statistics");
	if( (ret_val = Pixel_Analyzer(nrows, ncols)) != 0){
		printf("\nError in Pixel Analyzer"); exit(31);
	}
	printf("\nNumber of non-missing pixels: %d", matrix_stats.Num_pixels);
	printf("\npixel counts and percents by type:");
	for(ind1 = 0; ind1 < matrix_stats.Num_types; ind1++) {
			printf("\n%4d %10d  %7.5f", In_to_out[ind1], matrix_stats.Pix_counts [ind1], matrix_stats.Pix_percents [ind1]);
	}
	// Identify the patches. 
	printf("\nIdentifying patches with %d neighbor rule", parameters.neighbor_rule);
	Patch_Definer(nrows, ncols);
	// print n patches 
	printf("\nNumber of patches by type:");
	for (index = 0; index < 256; index++){
		temp_int = matrix_stats.Num_type_patches[index];
		if(temp_int > 0) {
			printf("\n%4d  %10d", In_to_out[index], temp_int);
		}
	}
	printf("\nAll types: %d", matrix_stats.Num_patches);
	//average patch sizes 
	for (index = 0; index < 256; index++){
		temp_int = matrix_stats.Num_type_patches[index];
		if(temp_int > 0) {
			matrix_stats.Average_patch_size [index] = (1.0*matrix_stats.Pix_counts [index])/(1.0*matrix_stats.Num_type_patches[index]);
		}
	}
	matrix_stats.Average_patch_size_overall = (1.0* matrix_stats.Num_pixels)/(1.0*matrix_stats.Num_patches);
	printf("\nAverage patch size by type:");
	for (index = 0; index < 256; index++){
		temp_int = matrix_stats.Num_type_patches[index];
		if(temp_int > 0) {
			printf("\n%4d  %10.2f", In_to_out[index],matrix_stats.Average_patch_size [index] );
		}
	}
	printf("\nAll types: %10.2f", matrix_stats.Average_patch_size_overall);
		// write patch map if requested
	if(parameters.write_patch_map == 1) {
		if( (patfile = fopen(filename_pat, "wb") ) == NULL) {printf("\nError opening %s\n", filename_pat);exit(16);	
		}
		printf("\nWriting patch number map");
		if(fwrite(mat_patnum, 4, (nrows * ncols), patfile) != (ncols * nrows) ) {
			printf("\nError writing patch number map file.\n");	exit(24);
		}
	}
	// call the routine to analyze each patch
	if(parameters.patch_analyzer == 1) {
		// Allocate the memory for the patch statistics to come.  Need to malloc one 
		//   extra structure because the patches are indexed at base 1               
		patch_stats = (struct statistics_2 *)calloc((matrix_stats.Num_patches + 1), sizeof(struct statistics_2)); 
		if(patch_stats == NULL) {printf("\nError: Malloc failed in Patch_Analyzer\n");exit(8);} 
		if( (patstatsfile = fopen(filename_patstats, "w") ) == NULL) {
			printf("\nError opening %s\n", filename_patstats);exit(16);	
		}
		printf("\nAnalyzing patches");
		if(parameters.neighbor_rule == 4) {
			Patch_Analyzer_4neighbor(nrows, ncols);
			// write the patch stats file
			printf("\nWriting patch statistics file");
			fprintf(patstatsfile, "Patch_number,Patch_type,Num_pixels,Out_edges,In_edges,Out_pixels,Row_org,Col_org,Touch_flag,Row_max,Col_min,Col_max");
			for(index = 1; index <= matrix_stats.Num_patches; index++) {
				fprintf(patstatsfile, "\n%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d", index, In_to_out[(*(patch_stats+index)).Patch_type],
					(*(patch_stats+index)).Num_pixels, (*(patch_stats+index)).Out_edges, (*(patch_stats+index)).In_edges,
					(*(patch_stats+index)).Out_pixels, (*(patch_stats+index)).Row_org, (*(patch_stats+index)).Col_org,
					(*(patch_stats+index)).Touch_border, (*(patch_stats+index)).Row_max, (*(patch_stats+index)).Col_min,
					(*(patch_stats+index)).Col_max);
				
			}
		}
		if(parameters.neighbor_rule == 8) {
			Patch_Analyzer_8neighbor(nrows, ncols);
			// write the patch stats file
			printf("\nWriting patch statistics file");
			fprintf(patstatsfile, "Patch_number,Patch_type,Num_pixels,Row_org,Col_org,Row_max,Col_min,Col_max");
			for(index = 1; index <= matrix_stats.Num_patches; index++) {
				fprintf(patstatsfile, "\n%d,%d,%d,%d,%d,%d,%d,%d", index, In_to_out[(*(patch_stats+index)).Patch_type],
					(*(patch_stats+index)).Num_pixels, (*(patch_stats+index)).Row_org, (*(patch_stats+index)).Col_org,
					(*(patch_stats+index)).Row_max, (*(patch_stats+index)).Col_min, (*(patch_stats+index)).Col_max);
				
			}
		}
		fclose(patstatsfile);
		free(patch_stats);	
	}
	// clean up
	free(mat_in);
	free(mat_patnum);
	free(mat_util);
	printf("\nNormal Finish.\n");
	exit(0);
}
/* Patch_Definer: FIND AND CODE PATCHES IN THE MATRIX
    Matrix_1 comes in with up to 255 type codes, with 255 indicating missing
	Matrix_2 goes out with patch numbers, or zero if that pixel was missing
    Matrix_3 is used to code missing, border, and completed pixels--The
		offset of (+1,+1) relative to Matrix_1 is needed to handle testing
		borders during the Neighbor tests. 
	This routine uses a linked-list queue to keep track of the to-do list.
	This cannot be parallelized because of the linked list queue, 
		and because two threads could be writing to same memory location.
*/
void Patch_Definer(unsigned int nrows, unsigned int ncols) 
{ 
	unsigned int patch_num, r1, c1, t1, r, c;
	struct queue_1 *current_ptr; // queue_1 is defined in common
	
	// the magic, see main 
	unsigned char (*Matrix_1)[ncols] = (unsigned char (*)[ncols])mat_in;
	unsigned int (*Matrix_2)[ncols] = (unsigned int (*)[ncols])mat_patnum;
	unsigned char (*Matrix_3)[ncols] = (unsigned char (*)[ncols])mat_util;

// Initialize Matrix_2 to patch number 0 
// Initialize interior of Matrix_3 to NULL if it's not a missing pixel, or to X if it is a missing pixel 
	for (r = 0; r < nrows; r++) { 
		for (c = 0; c < ncols; c++) { 
			Matrix_2 [r] [c] = 0; 
			if(Matrix_1 [r] [c] != 255) {  // If it's not missing 
				Matrix_3 [r+1] [c+1] = '\0'; 
			} 
			else {   // If it is missing  
				Matrix_3 [r+1] [c+1] = 'X'; 
			} 
		} 
	} 
// Initialize one-pixel borders of Matrix_3   
	r = 0; 
	for(c = 0; c <= ncols + 1; c++) { 
		Matrix_3 [r] [c] = 'X';} 
	r = nrows + 1; 
	for(c = 0; c <= ncols + 1; c++) { 
		Matrix_3 [r] [c] = 'X';} 
	c = 0; 
	for(r = 0; r <= nrows + 1; r++) { 
		Matrix_3 [r] [c] = 'X';} 
	c = ncols + 1; 
	for(r = 0; r <= nrows + 1; r++) { 
		Matrix_3 [r] [c] = 'X';} 
// Find a member of a new patch   
	patch_num = 0; // patch numbers will start at 1
	for(r = 0; r < nrows; r++) { 
		for(c = 0; c < ncols; c++) { 
			// skip if part of a patch already 
			if (Matrix_2 [r] [c] != 0) { continue; } 
			// skip if it's a missing value  
			if (Matrix_3 [r+1] [c+1] == 'X') { continue; } 
			// This pixel is part of a new patch  
			patch_num++; 
			// Check for limit to number of patches.
			if(patch_num > 4294967295) { 
				printf("\nMore than 4,294,967,295 patches found; exceeds program capacity.\n"); 
				free(mat_in); free(mat_patnum); free(mat_util);
				exit(8); 
			}
			first_ptr = NULL;   // first_ptr is declared in the common area  
			// Seed the patch with the patch number.  Other pixels in this
			//	patch will get their patch number in Find_Neighbors.  
			Matrix_2 [r] [c] = patch_num; 
			t1 = Matrix_1 [r] [c];  // convert the type to int 
			matrix_stats.Num_type_patches [t1] ++; 
			// Find neighbors to this pixel, and clear this pixel  
			Find_Neighbors(ncols, r, c, t1, patch_num); 
			// Now iterate the neighbor function to clear this patch  
			// Use the queue of neighbors to guide which pixels to clear 
			// Find_Neighbor will:
			//	code patch members it finds as 'x' in Matrix_3
			//	add these pixels into the queue
			//	code patch members it clears as 'X' in Matrix_3
			//	assign patch numbers to Matrix_2  
			current_ptr = first_ptr; 
			while(current_ptr != NULL) { 
				// find next pixel to check  
			    r1 = (*current_ptr).Row; 
			    c1 = (*current_ptr).Col; 
			    // remove that pixel from the queue and free the memory 
			    first_ptr = (*current_ptr).next_ptr; // pull out first_ptr 
			    free(current_ptr); 
			    Find_Neighbors(ncols, r1, c1, t1, patch_num); // this may update first_ptr 
			    current_ptr = first_ptr; 
			} 
		} 
	} 
	matrix_stats.Num_patches = patch_num; 
return; 
} 
/* Find_Neighbors(ncols, row, col, type, patch): 
		Find and mark neighbors of pixel at (row, col) that are of
		type "type", and mark them as patch number "patch"
		The marks are made in Matrix_3 and then in a queue  
		Then clear the pixel at (row, col)
*/ 
void Find_Neighbors(unsigned int ncols, unsigned int row, unsigned int col, unsigned int type, unsigned int patch) 
{
	// the magic, see main 
	unsigned char (*Matrix_1)[ncols] = (unsigned char (*)[ncols])mat_in;
	unsigned int (*Matrix_2)[ncols] = (unsigned int (*)[ncols])mat_patnum;
	unsigned char (*Matrix_3)[ncols] = (unsigned char (*)[ncols])mat_util;

	unsigned int o_row, o_col;  // local indices for Matrix_3 are offset from those of Matrices_1 and _2 
	o_row = row + 1; 
	o_col = col + 1; 
	// 4 neighbor first
	// Check neighbor above  
	if(Matrix_3 [o_row - 1] [o_col] != 'X') {   //skip if it's cleared, missing, or border 
		if(Matrix_3 [o_row - 1] [o_col] != 'x') {   //skip if it's in the queue already 
		if(Matrix_1 [row - 1] [col] == type) {   // check its type  
			Matrix_3 [o_row - 1] [o_col] = 'x';  // mark for later pass  
			Matrix_2 [row - 1] [col] = patch;    // assign patch number 
			Add_Queue( (row-1), col ); 
		} 
	   } 
	} 
	// Check neighbor below  
	if(Matrix_3 [o_row + 1] [o_col] != 'X') {  //  see comments above  
	   if(Matrix_3 [o_row + 1] [o_col] != 'x') {    
		if(Matrix_1 [row + 1] [col] == type) { 
			Matrix_3 [o_row + 1] [o_col] = 'x'; 
			Matrix_2 [row + 1] [col] = patch; 
			Add_Queue( (row+1), col); 
		} 
	   } 
	} 
	// Check neighbor to the left  
	if(Matrix_3 [o_row] [o_col - 1] != 'X')  { 
	   if(Matrix_3 [o_row] [o_col - 1] != 'x') { 
		if(Matrix_1 [row] [col - 1] == type) { 
			Matrix_3 [o_row] [o_col - 1] = 'x';  
			Matrix_2 [row] [col - 1] = patch; 
			Add_Queue(row, (col-1) ); 
		} 
	   } 
	} 
	// Check neighbor to the right  
	if(Matrix_3 [o_row] [o_col + 1] != 'X')  { 
	   if(Matrix_3 [o_row] [o_col + 1] != 'x') { 
		if(Matrix_1 [row] [col + 1] == type) { 
			Matrix_3 [o_row] [o_col + 1] = 'x';  
			Matrix_2 [row] [col + 1] = patch; 
			Add_Queue(row, (col+1) ); 
		} 
	   } 
	} 
	// add more neighbors to queue for 8-neighbor rule
	if(parameters.neighbor_rule == 8) {
		// Check NW neighbor  
		if(Matrix_3 [o_row - 1] [o_col-1] != 'X') {   //skip if it's cleared, missing, or border 
			if(Matrix_3 [o_row - 1] [o_col-1] != 'x') {   //skip if it's in the queue already 
			 if(Matrix_1 [row - 1] [col-1] == type) {   // check its type  
				Matrix_3 [o_row - 1] [o_col-1] = 'x';  // mark for later pass  
				Matrix_2 [row - 1] [col-1] = patch;    // assign patch number 
				Add_Queue( (row-1), (col-1) ); 
			 } 
			} 
		} 
		// Check SW neighbor  
		if(Matrix_3 [o_row + 1] [o_col-1] != 'X') {  //  see comments above  
			if(Matrix_3 [o_row + 1] [o_col-1] != 'x') {    
			 if(Matrix_1 [row + 1] [col-1] == type) { 
				Matrix_3 [o_row + 1] [o_col-1] = 'x'; 
				Matrix_2 [row + 1] [col-1] = patch; 
				Add_Queue( (row+1), (col-1) ); 
			 } 
			} 
		} 
		// Check NE neighbor  
		if(Matrix_3 [o_row - 1] [o_col+1] != 'X') {   //skip if it's cleared, missing, or border 
			if(Matrix_3 [o_row - 1] [o_col+1] != 'x') {   //skip if it's in the queue already 
			 if(Matrix_1 [row - 1] [col+1] == type) {   // check its type  
				Matrix_3 [o_row - 1] [o_col+1] = 'x';  // mark for later pass  
				Matrix_2 [row - 1] [col+1] = patch;    // assign patch number 
				Add_Queue( (row-1), (col+1) ); 
			 } 
			} 
		} 
		// Check SE neighbor  
		if(Matrix_3 [o_row + 1] [o_col+1] != 'X') {  //  see comments above  
			if(Matrix_3 [o_row + 1] [o_col+1] != 'x') {    
			 if(Matrix_1 [row + 1] [col+1] == type) { 
				Matrix_3 [o_row + 1] [o_col+1] = 'x'; 
				Matrix_2 [row + 1] [col+1] = patch; 
				Add_Queue( (row+1), (col+1) ); 
			 } 
			} 
		}
	}
	// Code this pixel as cleared  
	Matrix_3 [o_row] [o_col] = 'X'; 
return; 
} 
void Add_Queue(unsigned int arg_1, unsigned int arg_2) 
{ 
	// pointer to the next pixel in the queue  
	struct queue_1 *new_ptr;   // struct defined in header  
	new_ptr = (struct queue_1 *)malloc(sizeof(struct queue_1)); 
	if (new_ptr == NULL) { printf("\nError: Malloc in Add_Queue\n"); exit(8); } 
	(*new_ptr).Row = arg_1; 
	(*new_ptr).Col = arg_2; 
	(*new_ptr).next_ptr = first_ptr;   // initialized in header  
	first_ptr = new_ptr;               // and kept there between calls 
return; 
} 
 
/*   *********************
     Read_Parameter_File
     *********************
*/
int Read_Parameter_File(FILE *parfile)
{
	char ch;
	unsigned int value, flag;
	value = 99;
	flag = 1;
	while(flag == 1) {
		if (fscanf(parfile,"%s %d", &ch, &value) != 2) {
			flag = 0;
			continue;
		}
		if((ch == 'm') || (ch == 'M')) {
			parameters.missing_value_code = value;
			continue;
		}
		if((ch == 'z') || (ch == 'Z')) {
			parameters.recode = value;
			continue;
		}
		if((ch == 'n') || (ch == 'N')) {
			parameters.neighbor_rule = value;
			continue;
		}
		if((ch == 'x') || (ch == 'X')) {
			parameters.write_patch_map = value;
			continue;
		}
		if((ch == 'p') || (ch == 'P')) {
			parameters.patch_analyzer = value;
			continue;
		}
		return(1);
	}
	if(value == 99) {
		return(1);
	}
	return(0);
}
/* Pixel_Analyzer: CALCULATE very basic PIXEL STATISTICS   ** 4-neighbor ONLY */ 
int Pixel_Analyzer(unsigned int nrows, unsigned int ncols) 
{   unsigned int t, index;
	unsigned int r, c;
    // magic
	unsigned char (*Matrix_1)[ncols] = (unsigned char (*)[ncols])mat_in;
    // Count the number of pixels by type, and the total number of pixels  
    for (r = 0; r < nrows; r++) { 
		for (c = 0; c < ncols; c++) { 
			index = Matrix_1 [r] [c]; 
			matrix_stats.Pix_counts [index] ++; 
			matrix_stats.Num_pixels ++;    // for now, includes missing data 
		}                            
	} 
	 // Reduce by the number of pixels of missing data  
	 matrix_stats.Num_pixels =  matrix_stats.Num_pixels -  matrix_stats.Pix_counts [255]; 
	// Calculate the percentages by type 
	for (t = 0; t < matrix_stats.Num_types; t++) { 
		 matrix_stats.Pix_percents [t] = (1.0 * matrix_stats.Pix_counts [t] ) / (1.0 * matrix_stats.Num_pixels) ; 
	} 
return(0); 
} 

/* Patch_Analyzer_4neighbor:  */ 
void Patch_Analyzer_4neighbor(unsigned int nrows, unsigned int ncols) 
{ 
	int p, r, c, test;               // local indices and flag 
	int r1, c1;                      // temporary holders of row and col
	int o_row, o_col;                // local indices for Matrix_3 
	int r_org, c_org, r_max, c_max, r_min, c_min; 
	int out_edge, out_pix, sum_in_edges, edge_temp;  // counters 
	// the magic, see main 
	unsigned char (*Matrix_1)[ncols] = (unsigned char (*)[ncols])mat_in;
	unsigned int (*Matrix_2)[ncols] = (unsigned int (*)[ncols])mat_patnum;
	unsigned char (*Matrix_3)[ncols] = (unsigned char (*)[ncols])mat_util;
	// Initialize variables that throw warnings
	p=0; r=0; c = 0; r_org=0; c_org=0; c_max = 0; r_max=0; r_min=0; c_min=0; r1=0; c1=0;
	// Initialize the (*(patch_stats + patch)).Num_pixels counter  
	 for(p = 0; p <= matrix_stats.Num_patches; p++) { 
		 (*(patch_stats + p)).Num_pixels = 0; 
	 } 
	// Count the number of pixels in each patch, using patch number extracted 
	//   from Matrix_2 as subscript.  Note that any missing cells have a 0  
	//   patch number in Matrix_2, so they don't overrun the 'real' patches  
	//   for which the subscripts start at base 1.        
	for(r = 0; r < control.Rows; r++) { 
		for (c = 0; c < control.Cols; c++) { 
			(*(patch_stats + Matrix_2 [r] [c] )).Num_pixels++; 
		} 
	} 
	// The walkaround routine, counting edges per patch etc.      
	// Start by initializing Matrix_3, the utility matrix with extra rows&cols   
	for(r = 0; r <= control.Rows + 1; r++) { 
		for(c = 0; c <= control.Cols + 1; c++) { 
			Matrix_3 [r] [c] = '\0'; 
		} 
	} 
	// Loop through the patches, again note subscript base 1  
    for(p = 1; p <= matrix_stats.Num_patches; p++) { 
		// Initialize some counters, statistics, and flags  
		out_edge = 0; out_pix = 0;  sum_in_edges = 0; 
		(*(patch_stats + p)).In_edges = 0; 
		(*(patch_stats + p)).Touch_border = 0; 
		// Find the first member of patch p. For the first patch start at row 0. For subsequent
		//		patches start at the row of the preceding patch -- patch (p-1).
		//		Note that the row_org for patch p-1 will be set by the time patch p is checked.
		if(p > 1) {
			r = (*(patch_stats + p - 1)).Row_org; 
		} 
		else { // Use the first row for the first patch 
			r = 0; 
		}       
		while (r < control.Rows) {   // keep looping thru the rows 
			for(c = 0; c < control.Cols; c++) { 
				if(Matrix_2 [r] [c] != p) { 
					continue; 
				}
				else {
					// The upper left pixel of patch p was found  
					// Save the location and type  
					(*(patch_stats + p)).Row_org = r;  // This is also Row_min 
					(*(patch_stats + p)).Col_org = c; 
					(*(patch_stats + p)).Patch_type = Matrix_1 [r] [c]; 
					// Seed some values for later use in walkaround  
					r_org = r; 
					c_org = c; 
					r_min = r; 
					c_min = c; 
					r_max = r; 
					c_max = c; 
					goto breakout; // This throws a lot of warnings
				}
			} 
			r++;   // continue with the next row  
		} 
		breakout: 
		// If patch size is 1 or 2, the required statistics are easy and can skip   
		//  the rest of this routine by using the following few statements   
		if((*(patch_stats + p)).Num_pixels < 3)  { 
			(*(patch_stats + p)).In_edges = 0; 
			// Handle the patch size = 1 case  
			if((*(patch_stats + p)).Num_pixels == 1)  { 
				(*(patch_stats + p)).Out_edges = 4; 
				(*(patch_stats + p)).Out_pixels = 4; 
				(*(patch_stats + p)).Row_max = (*(patch_stats + p)).Row_org; 
				(*(patch_stats + p)).Col_max = (*(patch_stats + p)).Col_org; 
				(*(patch_stats + p)).Col_min = (*(patch_stats + p)).Col_org; 
				// See if the patch touches the border 
				if( r == 0 || r == (control.Rows - 1) || c == 0 ||  
					c == (control.Cols - 1) ) {  
					(*(patch_stats + p)).Touch_border = 1;  
				} 
				else {
					// See if the patch touches a missing value 
					if (Matrix_2 [r-1] [c] == 0 || Matrix_2 [r+1] [c] == 0 || 
						Matrix_2 [r] [c-1] == 0 || Matrix_2 [r] [c+1] == 0 ) { 
							(*(patch_stats + p)).Touch_border = 1;  
					}
				}
			} 
			// Handle the patch size = 2 case   
			if((*(patch_stats + p)).Num_pixels == 2)  { 
				(*(patch_stats + p)).Out_edges = 6; 
				(*(patch_stats + p)).Out_pixels = 6; 
				// See if the first pixel of the patch touches the border  
				if( r == 0 || r == (control.Rows - 1) || c == 0 ||  
					c == (control.Cols - 1) ) {  
					(*(patch_stats + p)).Touch_border = 1;  
				}
				else {
					// See if the first pixel of the patch touches a missing value   
					if (Matrix_2 [r-1] [c] == 0 || Matrix_2 [r+1] [c] == 0 || 
						 Matrix_2 [r] [c-1] == 0 || Matrix_2 [r] [c+1] == 0 ) { 
							(*(patch_stats + p)).Touch_border = 1;  
					}
				}
				// Now find the second pixel of the patch as r1 & c1  
				// It can't be in a previous row, or in a previous column 
				if(Matrix_2 [r+1] [c] == p) {   // check next row  
					r1 = r+1; c1 = c;         	// temp save of location 
					(*(patch_stats + p)).Row_max = r1; 
					(*(patch_stats + p)).Col_max = c; 
				} 
				if(Matrix_2 [r] [c+1] == p) {   // check next column 
					r1 = r; c1 = c+1; 
					(*(patch_stats + p)).Row_max = r; 
					(*(patch_stats + p)).Col_max = c1; 
				} 
				(*(patch_stats + p)).Col_min = c; 
				// Repeat the border & missing touch tests for pixel #2 
				if( r1 == 0 || r1 == (control.Rows - 1) || c1 == 0 ||  
					c1 == (control.Cols - 1) ) {  
					(*(patch_stats + p)).Touch_border = 1;  
				} 
				else { 
				   if (Matrix_2 [r1-1] [c1] == 0 || Matrix_2 [r1+1] [c1] == 0 || 
					 Matrix_2 [r1] [c1-1] == 0 || Matrix_2 [r1] [c1+1] == 0 ) { 
						(*(patch_stats + p)).Touch_border = 1;  
				   } 
				} 
			} // end of patch size 2 handling
		} // end of patch size <3 handling 
		else {   // if patch size is >= 3 then do the rest of this routine 
		//Walkaround this patch 
			test = 0; 
			while (test == 0) { 
				look_up: 
				// prepare to break out of while loop when walkabout is complete 
				if(out_edge != 0 && r == r_org && c == c_org) { 
					goto done;
				}
				// ensure not the matrix top edge  
				if(r > 0) { 
					if(Matrix_2 [r-1] [c] == p) { 
						// cell above is the same patch  
						// move on  
						r = r - 1; 
						r_min = min(r, r_min); 
						goto look_left; 
					} 
					else { 
						// cell above is not the same patch   
						out_edge++; 
						// note translation of cell to Matrix_3 base, r+1, c+1  
						if(Matrix_3 [r] [c+1] == '\0') { 
							 // count pixel and code it counted   
							 out_pix ++; 
							 Matrix_3 [r] [c+1] = 'C'; 
						} 
						// if the cell above is missing, code Touch_border */ 
						if(Matrix_2 [r-1] [c] == 0) { 
							(*(patch_stats + p)).Touch_border = 1; 
						} 
					} 
				}
				// handle the matrix top edge */ 
				else { 
					// a border pixel can be reached only once, so count/code it
					out_edge++; 
					out_pix++; 
					Matrix_3 [r] [c+1] = 'C'; 
					// Note the patches that touch the matrix border  
					(*(patch_stats + p)).Touch_border = 1; 
				} 
				look_right: 
				// ensure not the matrix right edge  
				if(c < control.Cols -1) { 
					if(Matrix_2 [r] [c+1] == p) { 
						// cell at right is the same patch  
						// move on   
						c = c + 1; 
						c_max = max(c, c_max); 
						goto look_up; 
					} 
					else { 
						// cell at right is not the same patch  
						out_edge++; 
						// note translation of cell to Matrix_3 base, r+1, c+1 
						if(Matrix_3 [r+1] [c+2] == '\0') { 
							// count pixel and code it counted 
							out_pix ++; 
							Matrix_3 [r+1] [c+2] = 'C'; 
						} 
						// if the cell at right is missing, code touch_border 
						if(Matrix_2 [r] [c+1] == 0) { 
						   (*(patch_stats + p)).Touch_border = 1; 
						} 
					} 
				} 
				// handle the matrix right edge  
				else { 
					// a border pixel can be reached only once, so count/code it 
					out_edge++; 
					out_pix++; 
					Matrix_3 [r+1] [c+2] = 'C'; 
					// Note the patches that touch the matrix border 
					(*(patch_stats + p)).Touch_border = 1; 
				} 
				look_down: 
				// ensure not the matrix bottom edge  
				if(r < control.Rows -1) { 
					if(Matrix_2 [r+1] [c] == p) { 
						// cell below is the same patch   
						// move on  
						r = r + 1; 
						r_max = max(r, r_max); 
						goto look_right; 
					} 
					else { 
						// cell below is not the same patch  
						out_edge++; 
						// note translation of cell to Matrix_3 base, r+1, c+1   
						if(Matrix_3 [r+2] [c+1] == '\0') { 
							// count pixel and code it counted   
							out_pix ++; 
							Matrix_3 [r+2] [c+1] = 'C'; 
						} 
						// if the cell below is missing, code touch_border   
						if(Matrix_2 [r+1] [c] == 0) { 
						   (*(patch_stats + p)).Touch_border = 1; 
						} 
					} 
				} 
				// handle the matrix bottom edge  
				else { 
					// a border pixel can be reached only once, so count/code it  
					out_edge++; 
					out_pix++; 
					Matrix_3 [r+2] [c+1] = 'C'; 
					// Note the patches that touch the matrix border   
					(*(patch_stats + p)).Touch_border = 1; 
				} 
				look_left: 
				// ensure not the matrix left edge  
				if(c > 0) { 
					if(Matrix_2 [r] [c-1] == p) { 
						// cell at left is the same patch   
						// move on  
						c = c - 1; 
						c_min = min(c, c_min); 
						goto look_down; 
					} 
					else { 
						// cell at left is not the same patch  
						out_edge++; 
						// note translation of cell to Matrix_3 base, r+1, c+1   
						if(Matrix_3 [r+1] [c] == '\0') { 
							 // count pixel and code it counted  
							 out_pix ++; 
							 Matrix_3 [r+1] [c] = 'C'; 
						} 
						// if the cell at left is missing, code touch_border  
						if(Matrix_2 [r] [c-1] == 0) { 
						   (*(patch_stats + p)).Touch_border = 1; 
						} 
					} 
				} 
				// handle the matrix left edge  
				else { 
					// a border pixel can be reached only once, so count/code it  
					out_edge++; 
					out_pix++; 
					Matrix_3 [r+1] [c] = 'C'; 
					// Note the patches that touch the matrix border   
					(*(patch_stats + p)).Touch_border = 1; 
				} 
				goto look_up; 
				done: 
				test = 1; 
			} // end of while loop  
			// Save the patch perimeter, out_pixels, and max cols and rows   
			(*(patch_stats + p)).Out_edges = out_edge; 
			(*(patch_stats + p)).Out_pixels = out_pix; 
			(*(patch_stats + p)).Row_max = r_max; 
			(*(patch_stats + p)).Col_max = c_max; 
			(*(patch_stats + p)).Col_min = c_min; 
			// Start peekinside routine  
			//  First, note that no patch with size < 8 can have an interior, so set  
			//    the statistics for those small ones here, and skip the peekinside routine 
			if( (*(patch_stats + p)).Num_pixels < 8) { 
				(*(patch_stats + p)).In_edges = 0; 
			} 
			else { 
				// Search for interior pixel: one enclosed by another patch   
				// Note that Matrix_3 has 'C' if it was an outside pixel, else has '\0'   
				// During the inside walkarounds, M_3 is set to 'Y' for pixels that have been dealt with 
				//   The Y's will be set back to '\0' at the end of this patch   
				// Code the pixels of patch "p" in Matrix_3 as 'X'   
				// Consider only the area within the min,max window for this patch 
				for(r = r_min; r <= r_max; r++) { 
					for(c = c_min; c <= c_max; c++) { 
						if(Matrix_2 [r] [c] == p) Matrix_3 [r+1] [c+1] = 'X'; 
					} 
				} 
				// Start the search for interior pixel  
				// Consider only the area within the min, max window for this patch   
				r = r_min - 1; 
				c = c_min - 1; 
				start_over:          // this routine is repeated for multiple inclusions  
				while(r <= r_max) {   //  search each row   
					test = 0; 
					while(c <= c_max) {  // search each column in the row  
						o_row = r + 1;    // translate to Matrix_3 base   
						o_col = c + 1; 
						if(Matrix_3 [o_row] [o_col] == '\0' && test == 0) { 
							// a non-patch, non-outside pixel found before the patch 
							goto next_col; 
						} 
						if(Matrix_3 [o_row] [o_col] == 'C' && test == 0) { 
							// this is first outside pixel found   
							test = 1; 
							goto next_col; 
						} 
						if(Matrix_3 [o_row] [o_col] == 'X' && test >= 1) { 
							// this is a pixel in the patch after an outside pixel found   
							test = 2; 
							goto next_col; 
						} 
						if(Matrix_3 [o_row] [o_col] == 'C' && test == 2) { 
							// a second outside pixel found before an inside pixel   
							test = 1; 
							goto next_col; 
						} 
						// deal with the case of finding the same inside pixel twice  
						// note that Matrix_3 was set to Y during the inside walkarounds   
						if(Matrix_3 [o_row] [o_col] == 'Y') { 
							// This pixel has already been dealt with, so increment columns   
							// until the next pixel in the patch is found, then continue search 
							// The next outside pixel is after the inclusion you want to ignore  
							//  and so you will pass Y and maybe \0 on the way to finding an X * 
							while(c <= c_max) { 
								c++; 
								if(Matrix_3 [r+1] [c+1] == 'X') { 
									test = 2; 
									goto next_col; 
								} 
							} 
							// the end of the row was reached before the next outside pixel, go to next row  
							goto next_row; 
						} 
						if(Matrix_3 [o_row] [o_col] == '\0' && test == 2) { 
							// inside pixel found  
							goto inside_found; 
						} 
						next_col: 
						c++; 
					} 
					// no inside found, now loop to search the next row  
					next_row: 
					r++; 
					// reset the column index  
					c = c_min - 1; 
				} 
				// end of search for interior pixels for this patch  
				goto reset; 
				inside_found: 
				// Set some flags and location markers   
				r_org = r; 
				c_org = c; 
				test = 0; 
				edge_temp = 0; 
				// code this inside pixel as having been found and dealt with  
				Matrix_3 [r+1] [c+1] = 'Y'; 
				// Walk-around the inside pixels   
				while(test == 0) { 
					look_up2: 
					if(edge_temp != 0 && r == r_org && c == c_org) { 
						// the transit is complete, check for next inclusion in this patch   
						goto next_incl;  
					} 
					if(Matrix_3 [o_row - 1] [o_col] != 'X') { 
						// pixel above is not in patch, move there   
						o_row = o_row - 1; 
						// and code it as having been dealt with   
						Matrix_3 [o_row] [o_col] = 'Y'; 
						r = r - 1; 
						goto look_left2; 
					} 
					else { 
						// pixel above is in patch, so count this edge  
						edge_temp++; 
						// and look right   
					} 
					look_right2: 
					if(Matrix_3 [o_row] [o_col + 1] != 'X') { 
						// pixel at right is not in patch, move there  
						o_col = o_col + 1; 
						// and code it as having been dealt with   
						Matrix_3 [o_row] [o_col] = 'Y'; 
						c = c + 1; 
						goto look_up2; 
					} 
					else { 
					// pixel at right is in patch, so count this edge   
					edge_temp++; 
					// and look down   
					} 
					look_down2: 
					if(Matrix_3 [o_row + 1] [o_col] != 'X') { 
						// pixel below is not in patch, move there   
						o_row = o_row + 1; 
						// and code it as having been dealt with   
						Matrix_3 [o_row] [o_col] = 'Y'; 
						r = r + 1; 
						goto look_right2; 
					} 
					else { 
					// pixel below is in patch, so count this edge   
						edge_temp++; 
					// and look left   
					} 
					look_left2: 
					if(Matrix_3 [o_row] [o_col - 1] != 'X') { 
						// pixel at left is not in patch, move there   
						o_col = o_col - 1; 
						// and code it as having been dealt with   
						Matrix_3 [o_row] [o_col] = 'Y'; 
						c = c - 1; 
						goto look_down2; 
					} 
					else { 
					// pixel at left is in patch, so count this edge   
						edge_temp++; 
					// and look up  
					} 
					goto look_up2; 
					next_incl: 
					// reset the flag to jump out of the while loop for this inclusion 
					test = 1; 
				}  // end of while inside walkaround
				//  Accumulate the counts for this patch, adding in the inclusion just completed  
				sum_in_edges = sum_in_edges + edge_temp; 
				// Update the number of inside pixels for patches of all sizes   
				(*(patch_stats + p)).In_edges = sum_in_edges; 
				// return to the peekinside routine to search for next inclusion   
				r = r_org;  // Go back to first column of current row   
				c = c_min-1; 
				goto start_over; 
			} // end search for interior pixel
			reset:  // landing zone for unsuccessful search for inside pixels  
			// End of ins and outs for this patch. Reset Matrix_3 for next patch   
			// Only need to reset within the min, max window for this patch.   
			// Again, note translation of indices to base for Matrix_3   
			r = r_min; 
			while (r <= r_max + 2) { 
				for (c = c_min; c <= c_max+2; c++) { 
					Matrix_3 [r] [c] = '\0'; 
				} 
				r++; 
			} 
		}  // Closing brace for if size >= 3   
		// That's all for this patch  
	} // closing brace for patches
	return; 
}  






/* Patch_Analyzer_8neighbor:  */ 
void Patch_Analyzer_8neighbor(unsigned int nrows, unsigned int ncols) 
{ 
	int p, r, c, cmin, cmax, rmin, rmax, test;
	// the magic, see main 
	unsigned char (*Matrix_1)[ncols] = (unsigned char (*)[ncols])mat_in;
	unsigned int (*Matrix_2)[ncols] = (unsigned int (*)[ncols])mat_patnum;
	cmin=0; cmax=0; rmax=0; rmin=0;
	// Initialize the (*(patch_stats + patch)).Num_pixels counter  
	for(p = 0; p <= matrix_stats.Num_patches; p++) { 
		(*(patch_stats + p)).Num_pixels = 0; 
	} 
	// Count the number of pixels in each patch, using patch number extracted 
	//   from Matrix_2 as subscript.  Note that any missing cells have a 0  
	//   patch number in Matrix_2, so they don't overrun the 'real' patches  
	//   for which the subscripts start at base 1.        
	for(r = 0; r < control.Rows; r++) { 
		for (c = 0; c < control.Cols; c++) { 
			(*(patch_stats + Matrix_2 [r] [c] )).Num_pixels++; 
		} 
	} 
	// find the patch type and origin row and column
	// Loop through the patches, again note subscript base 1  
    for(p = 1; p <= matrix_stats.Num_patches; p++) { 
		// Find the first member of patch p. For the first patch start at row 0. For subsequent
		//		patches start at the row of the preceding patch -- patch (p-1).
		//		Note that the row_org for patch p-1 will be set by the time patch p is checked.
		if(p > 1) {
			r = (*(patch_stats + p - 1)).Row_org; 
		} 
		else { // Use the first row for the first patch 
			r = 0; 
		}       
		while (r < control.Rows) {   // keep looping thru the rows 
			for(c = 0; c < control.Cols; c++) { 
				if(Matrix_2 [r] [c] != p) { 
					continue; 
				}
				else {
					// The upper left pixel of patch p was found  
					// Save the location and type  
					(*(patch_stats + p)).Row_org = r;  // This is also Row_min 
					(*(patch_stats + p)).Col_org = c; 
					(*(patch_stats + p)).Patch_type = Matrix_1 [r] [c]; 
					goto breakout; 
				}
			} 
			r++;   // continue with the next row  
		}
		breakout:
	}
	// find cmin, cmax, rmax for each patch
// OMP over patches		
#pragma omp parallel  for  	 private ( p, rmin, rmax, cmin, cmax, c, r, test)		
	for(p = 1; p <= matrix_stats.Num_patches; p++) { 
		// check the remaining columns of row = rmin = row_org
		rmin =(*(patch_stats + p)).Row_org; // row_org is by definition rmin
		rmax = rmin;
		cmin = (*(patch_stats + p)).Col_org; // Col_org is the cmin within the first row of the patch
		cmax = cmin;
		// check the rest of the row for cmax
		for(c = cmin+1; c < control.Cols; c++) {
			if(Matrix_2 [rmin] [c] == p) {
				cmin = min(c, cmin);
				cmax = max(c, cmax);
			}
		}
		// now check the remaining rows
		test = 0; 
		r = rmin;
		while(test == 0){
			r++;
			test = 1; // presuppose the row has none, in which case there can be none in the following rows either
			for(c = 0; c < control.Cols; c++) { 
				if(Matrix_2 [r] [c] == p) {
					test = 0; //oops the row has one, need to test at least one more row
					rmax = r;
					cmin = min(c, cmin);
					cmax = max(c, cmax);
				}	
			}
		}
		(*(patch_stats + p)).Row_max = rmax;
		(*(patch_stats + p)).Col_max = cmax;
		(*(patch_stats + p)).Col_min = cmin;	
	}
	return; 
}  
