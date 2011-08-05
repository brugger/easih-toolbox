/******************************************************************************
 * 
 * Extracts PF reads from illumina tiles and adjust the QV to the standard.
 * 
 * 
 *
 *
 * Kim Brugger (05 Aug 2011), contact: kim.brugger@easih.ac.uk
 *****************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <strings.h>
#include <malloc.h>


#define dprintf if (0) printf
#define uprintf if (0) printf



int main(int argc, char **argv){

  FILE *fh =  fopen (argv[1], "r");
  size_t lineln = 10;
  char *line = malloc(sizeof(char)*lineln);
  char data[11][256];

  printf("Hello world, looking in file '%s'\n", argv[1]);
  
  while ( getline (&line, &lineln, fh) > 0 ) {

  
    int i;
    int j = 0;
    int k = 0;
    for(i=0; i <  lineln; i++) {
      if ( line[i] == '\t' ) {
	data[j][k++] = 0x00;
	j++; 
	k=0;
	continue;
      }
      data[j][k++] = line[i];
      //    printf("%d%d == %c\n", j, k, line[i] );
    }
    
    //    push @read, "\@${instr}_$run_id:$lane:$tile:$x:$y/$read\n";
    
    if ( data[10][0] == '0')
      continue;

    //    printf("%s", line);

    for(i=0;data[8][i];i++) {
      if (data[8][i] == '.')
	data[8][i] = 'N';
      
    }

    for(i=0;data[9][i];i++) {
	data[9][i] = data[9][i]-33;
    }
    
    if (data[7][0] == '3' )
      data[7][0] = '2';

    printf("@%s_%s:%s:%s:%s:%s/%s\t%s\t%s\n", &data[0][0], &data[1][0], &data[2][0], &data[3][0], &data[4][0], &data[5][0], &data[7][0], &data[8][0], &data[9][0], &data[10][0]);
    //    break;
  }
  return(0);
}
