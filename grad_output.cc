#include "include/grad_output.h"
#include <string.h>
#include <stdio.h>

void Grad_Output::parse_outputs(const char *out){
    int len = strlen(out);
    char outputs[len];
    strcpy(outputs,out);
    printf("%s",outputs);
    char *tok = NULL;
    tok = strtok(outputs, ",");
    while (tok){
      if (strcmp(tok,"test")==0){
            printf("output is %s\n", tok);
        }
        tok = strtok(NULL,",");
    }
}
