/******************************************************************************
 * 
 * Extracts PF reads from illumina tiles and adjust the QV prints result to stdout.
 * 
 *
 * Kim Brugger (05 Aug 2011), contact: kim.brugger@easih.ac.uk
 *****************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>


int main(int argc, char **argv){

  FILE *fh =  fopen (argv[1], "r");
  size_t lineln = 10;
  char *line = malloc(sizeof(char)*lineln);
  int index[20];
  
  while ( getline (&line, &lineln, fh) > 0 ) {

    int i;
    int j = 0;
    
    for(i=0; i <  lineln; i++) {
      if ( line[i] == '\t' ) {
	line[i] = 0x00;
	index[j++] = i+1;
      }
      else if ( line[i] == '\n' ) {
	line[i] = 0x00;
      }
    }

    //    if ( line[ index[9] ] == '0')
    //      continue;

    // Change read nr to 2 if it is 2 as read 2 is the barcode read.
    if (line[ index[6]] == '3' )
      line[ index[6]] = '2';

    for(i=0;line[ index[7] + i];i++) {
      if (line[ index[7] + i] == '.')
	line[ index[7] + i] = 'N';
    }

    int QV30 = 0;
    for(i=0;line[ index[8] + i];i++) {
	line[ index[8] + i] -= 31;
	if (line[ index[8] + i] >= 30 + 33)
	  QV30++;
    }


    printf("@%s_%s:%s:%s:%s:%s/%s\t%s\t%s\t%s\t%d\n", &line[0], &line[ index[0]], &line[ index[1]], &line[ index[2]], &line[ index[3]], &line[ index[4]], &line[ index[6]], &line[ index[7]], &line[ index[8]], &line[ index[9]], QV30);

  }
  return(0);
}


/*
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
*/
