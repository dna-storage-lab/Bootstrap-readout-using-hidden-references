#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void deletion_random(int cover, int num, const char *in_fastq, const char *outdir) {
    char buff[600] = {0};
    char subname[300] = {0};
    FILE *fseq = NULL;
    FILE *fsub = NULL;
    
    fseq = fopen(in_fastq, "rt");  
    if (fseq == NULL) {
        printf("Cannot open the file %s!\n", in_fastq);
        return;
    }
    
    int maxLines = 3000000;
    char **lines = malloc(maxLines * sizeof(char *));
    if (lines == NULL) {
        printf("Memory allocation failed!\n");
        fclose(fseq);
        return;
    }
    
    int lineCount = 0;
    while (fgets(buff, sizeof(buff), fseq) && lineCount < maxLines) {
        fgets(buff, sizeof(buff), fseq); 
        lines[lineCount] = malloc(strlen(buff) + 1);
        if (lines[lineCount] == NULL) {
            printf("Memory allocation for line failed at line %d\n", lineCount);
            break;
        }
        strcpy(lines[lineCount], buff);
        lineCount++;
        fgets(buff, sizeof(buff), fseq);  
        fgets(buff, sizeof(buff), fseq); 
    }
    
    fclose(fseq);
    
    srand(123);  
    
    for (int i = 0; i < num; i++) {
        sprintf(subname, "%s/sub_reads_%d.txt", outdir, i + 1);  
        fsub = fopen(subname, "wt");
        if (fsub == NULL) {
            printf("Cannot open file %s for writing!\n", subname);
            continue;
        }

        for (int j = 0; j < cover; j++) {
            int randomIndex = rand() % lineCount;
            fprintf(fsub, "%s", lines[randomIndex]);
        }
        
        fclose(fsub);
    }
    
    for (int i = 0; i < lineCount; i++) {
        free(lines[i]);
    }
    free(lines);
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        printf("Usage: %s <InputFastq> <NumReads> <ExpNum> <SubDir>\n", argv[0]);
        return 1;
    }
    
    const char *in_fastq = argv[1];   // Input FASTQ file
    int numReads = atoi(argv[2]);     // Number of reads to sample per subset
    int expNum = atoi(argv[3]);       // Number of subsets to generate
    const char *subDir = argv[4];     // Directory for output files
    
    deletion_random(numReads, expNum, in_fastq, subDir);
    return 0;
}
