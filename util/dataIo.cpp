
#include "dataIo.h"

void random_init(double *array, const int N, const int D){
    short unsigned seed[3];
    int i;

    #pragma omp parallel
    {
        int tId = omp_get_thread_num();
        seed[0]=(((tId*tId + 15) * 3)/7);
        seed[1]=(((tId*tId + 13) * 2)/3);
        seed[2]=tId;

        #pragma omp for
        for ( i=0 ; i < N * D ; i++ ){
            array[i] = erand48(seed) * 100.0;
        }
    }
}

void read_file(double *array, const int N, const int D, char filename[], const bool IS_BINARY){
    FILE *fp;
    size_t counts = 0;
    size_t i=0,j=0;
    char line[MAX_LINE_LENGTH];
    char *token=NULL;
    const char space[2] = " ";

    fp = fopen(filename,"r");

    if ( fp == NULL ){
        fprintf(stderr, "File '%s' does not exists!", filename);
        exit(1);
    }

    if ( IS_BINARY ){
        // printf("processing binary file!");fflush(stdout);
        // read binary file, everything at once
        counts = fread(array, N * D, sizeof(double), fp);
        // printf("%dx%d: %d readed\n", N, D, counts);
        if ( counts == 0 ) {
            fprintf(stderr, "Binary file '%s' could not be read. Wrong format.", filename);
            exit(1);
        }
    }else{
        // processing a text file
        // format: there are D double values each line. Each value is separated by a space character.
        // notice MAX_LINE_LENGTH = 2049
        // printf("processing text file!");fflush(stdout);
        i = 0;
        while ( fgets ( line, MAX_LINE_LENGTH, fp ) != NULL &&
                i < N ) {
            if ( line[0] != '%'){ // ignore '%' comment char
                token = strtok(line, space);
                j=0;

                while ( token != NULL &&
                        j < D ){
                    array[i*D + j] = atof(token); // 0.0 if no valid conversion
                    token = strtok(NULL, space);
                    j++;
                }
                i++;
            }
        }
    }

    fclose(fp);
}

void save_binary_file(double *array, const int N, const int D, char filename[]){
    FILE *fp=NULL;
    size_t counts = 0;

    fp = fopen(filename, "w");

    if ( fp == NULL ){
        fprintf(stderr, "Could not open file '%s'!", filename);
        exit(1);
    }

    counts = fwrite(array,sizeof(double) * N * D, 1, fp);

    if ( counts == 0 ){
        fprintf(stderr, "Error in writing file '%s'. Abort.", filename);
        exit(1);
    }

    fclose(fp);
}

void save_text_file(double *array, const int N, const int D, char filename[]){
    FILE *fp=NULL;
    size_t counts = 0;
    size_t i=0, j=0;
    char line[MAX_LINE_LENGTH];
    char strDouble[50];

    fp = fopen(filename, "w");

    if ( fp == NULL ){
        fprintf(stderr, "Could not open file '%s'!", filename);
        exit(1);
    }

    for ( i=0 ; i < N ; i++ ){
        strcpy(line, "");
        for ( j=0 ; j < D ; j++ ){
            strcpy(strDouble, "");
            sprintf(strDouble, "%f ", array[i*D + j]);
            strcat(line, strDouble);

        }
        fprintf(fp, "%s\n", line);

    }

    fclose(fp);
}
