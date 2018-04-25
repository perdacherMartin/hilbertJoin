

#include "arguments.h"


void parsing_args(int argc, char* argv[], size_t *n, double *epsilon, size_t *d, size_t *threads, char *filename, bool *isBinary, int *KBLOCK, int *stripes){
  char c;
  FILE *file;

  if ( argc < 5 ){
    fprintf (stderr, "There are obligatory parameters.\n");
    fprintf (stderr, "Usage: ./egoHilb (or ./egoCano)");

    fprintf(stderr, "Obligatory parameters: \n");
    fprintf(stderr, "n (number of objects )\ne (epsilon)\nd (dimensionality)\n");
    fprintf(stderr, "Optional parameters: \n t number of threads\n\n");
    fprintf(stderr, "k KBLOCK (default 8)\n");
    fprintf(stderr, "s number of stripes (default 2)\n");
    fprintf(stderr, "f (filename) if there is no filename we use random generated data [0.0, 100.0)\n");
    fprintf(stderr, "b use the -b argument without options to specify that it is a binary file.\n");
    fprintf(stderr, "Example (with default values): ./blasJoin -n 200 -e 0.2 -d 20 -s 5000 -t 64\n");
    exit(1);
  }

  while ( (c = getopt(argc, argv, "n:e:d:t:f:k:s:b") ) != -1) {

	if ( optarg ){
        switch(c){
            case 'n':
    		    *n = (size_t)atol(optarg);
    		    break;
    		case 't':
    		    *threads = atoi(optarg);
    			break;
            case 'd':
                *d = (size_t)atoi(optarg);
                break;
            case 'e':
                *epsilon = atof(optarg);
                break;
            case 'k':
                *KBLOCK = atoi(optarg);
                break;
            case 's':
                *stripes = atoi(optarg);
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
            case 'b':
                *isBinary = true;
                break;
            case '?':
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                exit(1);
                break;
            default:
                break;
        }
    }
  }

  // if ( *blocksize > *n || *blocksize <= 1 ){
  //     fprintf (stderr, "Blocksize has to be greater than 1 and smaller or equal to N\n");
  //     printf("n:%d, blocksize: %d\n", *n, *blocksize);
  //     exit(1);
  // }

}
