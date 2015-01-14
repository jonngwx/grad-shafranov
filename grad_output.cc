#include "include/grad_output.h"
#include <string.h>
#include <stdio.h>

void GradOutput::parse_outputs(const char *out){
    int len = strlen(out);
    if (len < 1){
        return;
    }
    char outputs[len];
    strcpy(outputs,out);
    char *tok = NULL;
    //    printf("%s\n",outputs);
    tok = strtok(outputs, ",");
    while (tok){
      if (strcmp(tok,"j")==0){
          printf("output is %s\n", tok);
          if (!find(CURRENT)) {
              output_list.push_back(CURRENT);
          } else {
              printf("already there!\n");
          }
      } else if (strcmp(tok,"bt") == 0){
          printf("bt\n");
          if (!find(TOROIDAL_FIELD)) {
              output_list.push_back(TOROIDAL_FIELD);
          } else {
              printf("bt already there!\n");
          }
      }
      tok = strtok(NULL,",");
    }
}
