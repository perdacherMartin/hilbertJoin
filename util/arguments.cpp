

#include "arguments.h"


void parsing_args(int argc, char* argv[], size_t *n, double *epsilon, size_t *d, char *filename, int *activedims){
  char c;
  FILE *file;

  if ( argc < 5 ){
    fprintf (stderr, "There are obligatory parameters.\n");
    fprintf (stderr, "Usage: ./hilbertSelfJoinCardinality");

    fprintf(stderr, "Obligatory parameters: \n");
    fprintf(stderr, "n (number of objects )\ne (epsilon)\nd (dimensionality)\n");
    fprintf(stderr, "Optional parameters: \n\n");
    fprintf(stderr, "a number of acitve dimensions (default 3)\n");
    fprintf(stderr, "f (filename) if there is no filename we use random generated data [0.0, 1.0)\n");
    // fprintf(stderr, "b use the -b argument without options to specify that it is a binary file.\n");
    fprintf(stderr, "Example (with default values): ./hilbertSelfJoinCardinality -n 200000 -e 0.2 -d 64 -t 64\n");
    exit(1);
  }

  while ( (c = getopt(argc, argv, "n:e:d:t:f:k:s:") ) != -1) {

	if ( optarg ){
        switch(c){
            case 'n':
    		    *n = (size_t)atol(optarg);
    		    break;
            case 'd':
                *d = (size_t)atoi(optarg);
                break;
            case 'e':
                *epsilon = atof(optarg);
                break;
            case 's':
                *activedims = atoi(optarg);
                break;
            case 'f':
                strcpy(filename, optarg);
                break;
            case '?':
              if (optopt == 'c')
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
              else if (isprint (optopt)){
                  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                  exit(1);
              }
              else
                fprintf (stderr,
                         "Unknown option character `\\x%x'.\n",
                         optopt);
              exit(1);
    		default:
    			break;
        }
	}else{
        switch(c){
            // case 'b':
            //     *isBinary = true;
            //     break;
            case '?':
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                exit(1);
                break;
            default:
                break;
        }
    }
  }

  if ( *activedims < 0 || *activedims > 6 ){
    fprintf(stderr, "parameter active dimensions is typical in the range of [0-5]\n");
    exit(1);
  }

}


void parsing_args_join(int argc, char* argv[], size_t *n, size_t *m, double *epsilon, size_t *d, char *filename, char *filename2, int *activedims){
  char c;
  FILE *file;

  if ( argc < 5 ){
    fprintf (stderr, "There are obligatory parameters.\n");
    fprintf (stderr, "Usage: ./hilbertJoinCardinality");

    fprintf(stderr, "Obligatory parameters: \n");
    fprintf(stderr, "n (number of objects in set A)\nm (number of objects in set B)\ne (epsilon)\nd (dimensionality)\n");
    fprintf(stderr, "Optional parameters: \n\n");
    fprintf(stderr, "a number of active dimensions (default 3)\n");
    fprintf(stderr, "f (filename) if there is no filename we use random generated data [0.0, 1.0)\n");
    fprintf(stderr, "g (filename set B) if there is no filename we use random generated data [0.0, 1.0)\n");
    // fprintf(stderr, "b use the -b argument without options to specify that it is a binary file.\n");
    fprintf(stderr, "Example (with default values): ./hilbertJoinCardinality -n 200000 -m 200000 -e 0.2 -d 20 -t 64\n");
    exit(1);
  }

  while ( (c = getopt(argc, argv, "n:m:e:d:t:f:g:k:a:s:") ) != -1) {

	if ( optarg ){
        switch(c){
            case 'n':
    		    *n = (size_t)atol(optarg);
    		    break;
            case 'm':
    		    *m = (size_t)atol(optarg);
    		    break;
            case 'd':
                *d = (size_t)atoi(optarg);
                break;
            case 'e':
                *epsilon = atof(optarg);
                break;
            case 's':
                *activedims = atoi(optarg);
                break;
            case 'f':
                strcpy(filename, optarg);
                break;
            case 'g':
                strcpy(filename2, optarg);
                break;
            case '?':
              if (optopt == 'c')
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
              else if (isprint (optopt)){
                  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                  exit(1);
              }
              else
                fprintf (stderr,
                         "Unknown option character `\\x%x'.\n",
                         optopt);
              exit(1);
    		default:
    			break;
        }
	}else{
        switch(c){
            // case 'b':
            //     *isBinary = true;
            //     break;
            case '?':
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                exit(1);
                break;
            default:
                break;
        }
    }
  }

  if ( *activedims < 0 || *activedims > 6 ){
    fprintf(stderr, "parameter active dimensions is typical in the range of [0-5]\n");
    exit(1);
  }

}
