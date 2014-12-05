#include "include/grad_output.h"
#include <string.h>
#include <stdio.h>

/**
 * \file Grad output base class implementation 
 * */



/**
 * @brief Parses the string of outputs to write to file and returns a list.
 * @param out comma separated string of outputs. currently allows p,bt,q
 * 
 */
void Grad_Output::parse_outputs(const char *out){
    int len = strlen(out);
    char outputs[len];
    strcpy(outputs,out);
    printf("%s",outputs);
    char *tok = NULL;
    tok = strtok(outputs, ",");
    while (tok){
        printf("output is %s\n", tok);
        tok = strtok(NULL,",");
    }
}
