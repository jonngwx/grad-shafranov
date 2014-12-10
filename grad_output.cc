#include "include/grad_output.h"
#include <string.h>
#include <stdio.h>

void Grad_Output::parse_outputs(const char *out){
    int len = strlen(out);
    if (len < 1){
        return;
    }
    char outputs[len];
    strcpy(outputs,out);
    char *tok = NULL;
    tok = strtok(outputs, ",");
    while (tok){
      if (strcmp(tok,"j")==0){
          printf("output is %s\n", tok);
      }
      tok = strtok(NULL,",");
    }
}
