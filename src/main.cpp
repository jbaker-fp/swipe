#include "swipe.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <climits>
#include <cstdbool>
#include <cstring>
#include <ctime>

#include <vector>

// feel free to change these defaults
#define ST       .1
#define DT       .1
#define MIN      4.0
#define MAX      43.0

#define VNUM    1.5 // current version

// main method, interfacing with user arguments
int main(int argc, char* argv[]) {
    char output[] = "OUTPUT:\npitch_0\ttime_0\npitch_1\ttime_1\n...\t...\
    \npitch_N\ttime_N\n\n"; 
    char header[] = "\nSWIPE' pitch tracker, implemented by Kyle Gorman \
<gormanky@ohsu.edu>, \nbased on: A. Camacho. 2007. A sawtooth \
waveform inspired pitch estimator\nfor speech and music. Doctoral \
dissertation, U of Florida.\n\n\
More information: <http://ling.upenn.edu/~kgorman/C/swipe/>\n\n";
    char synops[] = "SYNPOSIS:\n\n\
swipe [-i FILE] [-o FILE] [-b LIST] [-r MIN:MAX] [-s TS] [-t DT] [-mnhv]\n\
\nFLAG:\t\tDESCRIPTION:\t\t\t\t\tDEFAULT:\n\n\
-i FILE\t\tinput file\t\t\t\t\t<STDIN>\n\
-o FILE\t\toutput file\t\t\t\t\t<STDOUT>\n\
-b LIST\t\tbatch mode [LIST is a file containing\n\
\t\tone \"INPUT OUTPUT\" pair per line]\n\n\
-r MIN:MAX\tpitchrange in Hertz\t\t\t\t100:600\n\
-s TIMESTEP\ttimestep in seconds\t\t\t\t0.001\n\
-t THRESHOLD\tstrength threshold [0 <= x <= 1]\t\t0.300\n\n\
-m\t\tOutput Mel pitch\t\t\t\tno\n\
-n\t\tDon't output voiceless frames\t\t\tno\n\
-h\t\tDisplay this message, then quit\n\
-v\t\tDisplay version number, then quit\n\n";
    // all set by #defines
    double st = ST; 
    double dt = DT;
    double min = MIN;
    double max = MAX; 
    int ch;
    // some, but not all, compilers initialize char*s to be "\0"
    char* wav = NULL; 
    char* out = NULL;
    int needed;
    while ((ch = getopt(argc, argv, "i:f:o:r:s:t:b:mnhv")) != -1) {
        switch(ch) {
            case 'i':
                needed = (int) (strlen(optarg) + 1);
                if (needed > FILENAME_MAX) {
                    fprintf(stderr, "Filename too long, aborting.\n");
                    exit(EXIT_FAILURE);
                }
                wav = (char *) malloc(sizeof(char) * needed);
                strcpy(wav, optarg);
                break; 
            case 'o':
                needed = (int) (strlen(optarg) + 1);
                if (needed > FILENAME_MAX) {
                    fprintf(stderr, "Filename too long, aborting.\n");
                    exit(EXIT_FAILURE);
                }
                out = (char *) malloc(sizeof(char) * needed);
                strcpy(out, optarg);
                break;
            case 'r':
                min = atof(strtok(optarg, ":"));
                max = atof(strtok(NULL, ":")); 
                break;
            case 't':
                st = atof(optarg);
                break;
            case 's':
                dt = atof(optarg);
                break;
            case 'h':
                fprintf(stderr, "%s", header); 
                fprintf(stderr, "%s", synops); 
                fprintf(stderr, "%s", output);
                exit(EXIT_SUCCESS);
            case 'v':
                fprintf(stderr, "This is SWIPE', v. %1.1f.\n", VNUM); 
                exit(EXIT_SUCCESS);
            case '?':
            default: 
                fprintf(stderr, "%s", header);
                fprintf(stderr, "%s", synops);
                exit(EXIT_FAILURE);
            argc -= optind; 
            argv += optind;
        }
    }
    // santiny-check the args
    if (min < 1.) { 
        fprintf(stderr, "Min pitch < 1 Hz, aborting.\n"); 
        exit(EXIT_FAILURE);
    }
    if (max - min < 1.) {
        fprintf(stderr, "Max pitch <= min pitch, aborting.\n"); 
        exit(EXIT_FAILURE);
    } 
    if (st < 0. || st > 1.) { 
        fprintf(stderr, "Strength must be 0 <= x <= 1, set to %.3f.\n", ST); 
        st = ST;
    }
    if (dt < .001) {
        fprintf(stderr, "Timestep must be >= 0.001, set to %.3f.\n", DT);
        dt = DT;
    }

    FILE* in = fopen(wav, "r");
    if (in == NULL) {
        fprintf(stderr, "Reading from \"%s\" failed (try ", wav);
        fprintf(stderr, "specifying an input file with -i).\n");
        exit(EXIT_FAILURE);
    }

    // COunts lines
	int num_lines = 0;
	double tmp = 0;
	while(fscanf(in, "%lf\n", &tmp) != EOF){
		num_lines++;
	}
	rewind(in);
	

	std::vector<double> x(num_lines, 0);
	// in = fopen(argv[1], "r");
	for (int i = 0; i < x.size(); i++) {
		fscanf(in, "%lf\n", &(x[i]));
	}
	fclose(in);

    int fs = 4800;
    double derbs = .01;
    double dlog2p = 1./25.;
    int step = fs*dt;


    PitchVector p = swipe(x, fs, step, min, max, dlog2p, derbs, st);


    FILE* f_out = fopen(out, "w");
    fprintf(f_out, "{\"fo_strength\": [");
    for (int i = 0; i < p.strength.size(); i++) {
        fprintf(f_out, "%.40lf", p.strength[i]);
		if (i != p.strength.size()-1) {
			fprintf(f_out, ", ");
		}
    }


    fprintf(f_out, "], \"fo_times\": [");
    for (int i = 0; i < p.pitch.size(); i++) {
        fprintf(f_out, "%.1f", i*dt);
		if (i != p.pitch.size()-1) {
			fprintf(f_out, ", ");
		}
    }

    fprintf(f_out, "], \"fo\": [");
    for (int i = 0; i < p.pitch.size(); i++) {
        fprintf(f_out, "%.40lf", p.pitch[i]);
		if (i != p.pitch.size()-1) {
			fprintf(f_out, ", ");
		}
    }
    fprintf(f_out,"]}");
    fclose(f_out);

    free(out);
    free(wav);
    exit(EXIT_SUCCESS);
}