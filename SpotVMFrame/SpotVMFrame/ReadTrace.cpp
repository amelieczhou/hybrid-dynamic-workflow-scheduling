#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#include "ReadTrace.h"

int ReadTrace::readfile(char* filename){	    
    if((fp=fopen(filename, "r"))==NULL){
	printf("Error: file not found: %s\n", filename);
	exit(1);
    }    

    return 0;
}

int ReadTrace::readfile(char* filename, int l){	    
    readfile(filename);
	if(fp != NULL){
		char str[1024]; // tmp buf for one line
		for(int index = 0; index < l; index ++){
			char* ptr = fgets(str, 1024, fp);
		}
	}

	return 0;
}

int ReadTrace::closefile(){
	if(fp != NULL)
		fclose(fp);
	return 0;
}

int ReadTrace::readprice(float* result){	
	int j;
	
    char str[1024]; // tmp buf for one line
    char buf[32]; // tmp buf for one digit
    char *ptr, *ptr2;
    
	if(fp == NULL)
		return NULL;
    
	ptr = fgets(str, 1024, fp);
	if(ptr==NULL){	    
	    // no more results, so just break this loop
	    return NULL;
	}

	for(j=0;j<42;++j){
	    ptr2 = strchr(ptr, ',');
	    if(ptr2 == NULL){
			ptr2 = strchr(ptr, '\n');
			if (ptr2 == NULL){		    
				// we don't have more, so just break this loop
				break;
			}
	    }
	    
	    // store the float value
	    memset(buf, 0, 32);
	    strncpy(buf, ptr, ptr2-ptr);
	    if(j == 3)
			result[0] = atof(buf); 
		else if(j == 38)
			result[1] = atof(buf); 
		else if(j == 30)
			result[2] = atof(buf); 
		else if(j == 4)
			result[3] = atof(buf); 
	    
	    // move pointer to the next position
	    ptr = ptr2 + 1; 
	}

	return 0;
}
