/*
# Copyright 2007, 2008, 2009
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Moritz Gerstung,
# Lukas Geyrhofer,
# Osvaldo Zagordi,
# ETH Zurich

# This file is part of ShoRAH.
# ShoRAH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ShoRAH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ShoRAH.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DPM_SAMPLER_HPP
#define DPM_SAMPLER_HPP

/* variables and functions defined here
   for data structures see data_structures.hpp */

#include <map>
#include <string>
#include <utility>
#include "data_structures.hpp"

const unsigned int B = 5;  // characters in the alphabet
const int one_int = 153391689;
const int two_int = 306783378;
const int four_int = 613566756;
const int LIMIT = 100000;

//char* filein;
char** id;
unsigned short int** r;
crnode** readtable;
crnode** readtable2;
int** ftable;
int* ftable_sum;
// multimap<string,int> readmap;
// int** creads;
unsigned short int* h;
ssret* res;
int* res_dist;
int* cbase;
double* pbase;
double* log_pbase;

unsigned int n = 0;
unsigned int q = 0;
unsigned int totsites = 0;
unsigned int totbases = 0;
unsigned int hapbases = 0;
unsigned int lowest_free = 0;
unsigned int iter = 1000;
unsigned int MAX_K;
unsigned int J;
unsigned int K = 0;  // initial number of clusters
double avgNK = 0.0;  // average #reads in each startcluster, avgNK = n/K
double default_avgNK = 10.0;
unsigned int HISTORY = 100;
double theta = 0.90; // error rate
double eps1 = 0.985;
double eps2 = 0.001;
double gam = 0.90; // mutation rate
double alpha = 0.01; // probability of creating new classes
double g_noise = 0.0001;

double* P;
double* log_P;

/* Used to avoid double counting resulting in posteriors > 1 by keeping track of
 * which iteration/step a haplotype-count was last increased in. */
struct versioned_count {
    versioned_count() : count(0), last_seen_step(-1) {}

    int count;
    int last_seen_step;
};

typedef std::pair<std::string, int> sipair;
// bool comp_values (sipair e1, sipair e2);
typedef std::map<std::string, int> sup_map;
/* NOTE: 'sup_map' is used for read2hap, which appears to be used for output only.
 *       Since read2hap tracks read assignments, which are unique, no double counting
 *       issues can arise. */
typedef std::map<std::string, versioned_count> sup_map_versioned;
sup_map_versioned support; /*! How many times the haplotype/read is found in history*/

typedef std::map<std::string, int*> freq_map;
freq_map freq; /*! How many reads had the ahplotype assigned in history*/

typedef std::map<std::string, std::string*> ass_map;
ass_map assignment; /*! reads id assigned to the given haplotype*/
sup_map* read2hap;

cnode** cl_ptr;  // class local pointer, while c_ptr is global
cnode** c_ptr;
cnode* mxt = NULL;
hnode* hst = NULL;  // linked list for haplotype history, first element
// sample_rec* sample_hist;

FILE* assign;

hnode*** ass_hist;
unsigned int record = 0;

unsigned short int** haplotypes;
int* ch;  // count haplotypes

char* haplotype_output;

double double_threshold_min;  // will be set to = gsl_sf_log(DBL_MIN);
double double_threshold_max;  // will be set to = gsl_sf_log(DBL_MAX);
// where DBL_MIN is the smallest number != 0.0, a (double) random number generator
// can produce...
// needed for omitting underflow-errors of gsl-functions

// random number generator via gsl
unsigned long randseed = 0; /* random seed, if no command line parameter
                   -R given, set to current time */

bool read_data(const char* filein, std::ofstream& out_file);

void read_conversion(crnode* b, unsigned short int* a, int seq_length);

void conversion(int* b, unsigned short int* a, int seq_length);

int weight(const cnode* wn);

int weight_shift(const cnode* wn, unsigned int i, unsigned int removed);

void build_assignment(std::ofstream& out_file);

int isdna(char c);

unsigned int d2i(char c);

std::pair<int, int> seq_distance_new(int* A, crnode* B, int seq_length);

std::pair<int, int> seq_distance_rr(int* A, crnode* B, int seq_length);

std::pair<int, int> seq_distance(unsigned short int* a, unsigned short int* b, int seq_length);

ssret* sample_class(unsigned int i, unsigned int step);

void sample_hap(cnode* c);

double sample_ref();

void check_size(const cnode* cst, unsigned int n);

int count_classes(const cnode* cst);

void write_assignment(unsigned int it, unsigned int new_proposed, const cnode* tn);

void create_history(unsigned int k);

void record_conf(cnode* tn, unsigned int step);

void cleanup();

int compare(const void* a, const void* b);

int compare_hnss_seq(const void* a, const void* b);

int compare_hnss_count(const void* a, const void* b);

int isnewhaplotype(hnode** p, unsigned short int* h);

double setfinalhaplotype(unsigned int i);

void write_haplotype_frequencies(char* filename, unsigned int hcount);

void write_posterior_files(std::string instr);

void print_stats(std::ofstream& out_file, const cnode* cn, unsigned int J);

void reset_globals();

#endif
