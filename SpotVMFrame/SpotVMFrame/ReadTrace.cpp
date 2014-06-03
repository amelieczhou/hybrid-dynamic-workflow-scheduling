#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <omp.h>

#include "ReadTrace.h"

int ReadTrace::readfile(char* filename){	    
    if((fp=fopen(filename, "r"))==NULL){
	printf("Error: file not found: %s\n", filename);
	exit(1);
    }    

    return 0;
}

int ReadTrace::readfile(char* filename, int l){//read the l-th line 
    readfile(filename);
	if(fp != NULL){
		char str[1024]; // tmp buf for one line
		for(int index = 0; index < l; index ++){
			char* ptr = fgets(str, 1024, fp);
		}	
	}

	return 0;
}

int ReadTrace::readtomem(char* filename, float* result){
	//readfile(filename);
	std::ifstream fin(filename);
	char str[14400];
	char buf[32];
	char *ptr, *ptr2;
	int lines,cols;
	if(strcmp(filename,"data1.csv")==0||strcmp(filename,"data2.csv")==0){
		lines = 35971;
		cols = 4;
	}
	else if(strcmp(filename,"tp.dat")==0||strcmp(filename,"newData/tp.dat")==0){
		lines = 800;
		cols = 440;
	}else if(strcmp(filename,"ffp_2.dat")==0||strcmp(filename,"newData/ffp_2.dat")==0){
		lines = 2520;//760;
		cols = 2400;
	}
	int i=0;
	while(fin.getline(str,14400)&&i<lines){
		ptr=str;
		for(int j=0; j<cols; j++){
			//if(j==2399)
				//printf("");
			ptr2 = strchr(ptr,',');
			if(ptr2 == NULL){
				ptr2 = strchr(ptr,'\n');
				if(ptr2 == NULL) {
					memset(buf,0,32);
					strncpy(buf,ptr,5);
					result[i*cols+j] = atof(buf);
					break;
				}
			}
			memset(buf,0,32);
			strncpy(buf,ptr,ptr2-ptr);
			result[i*cols+j] = atof(buf);
			ptr = ptr2 + 1;
			
		//	printf("%d\n",j);
		}
		i++;
	}
	//if(fp!=NULL){
	//	int setbufsize = setvbuf (fp , NULL , _IOFBF , 14400);
	//	char str[14400];
	//	char buf[32];
	//	char *ptr, *ptr2;
	//	int lines,cols;
	//	if(strcmp(filename,"data1.csv")==0||strcmp(filename,"data2.csv")==0){
	//		lines = 35971;
	//		cols = 4;
	//	}
	//	else if(strcmp(filename,"tp.dat")==0||strcmp(filename,"newData/tp.dat")==0){
	//		lines = 800;
	//		cols = 440;
	//	}else if(strcmp(filename,"ffp_2.dat")==0||strcmp(filename,"newData/ffp_2.dat")==0){
	//		lines = 760;
	//		cols = 2400;
	//	}
	//	for(int i=0; i<lines; i++){//totally 35971 lines in the trace
	//		ptr = fgets(str,14400,fp);
	//		for(int j=0; j<cols; j++){
	//			ptr2 = strchr(ptr,',');
	//			if(ptr2 == NULL){
	//				ptr2 = strchr(ptr,'\n');
	//				if(ptr2 == NULL) break;
	//			}
	//			memset(buf,0,32);
	//			strncpy(buf,ptr,ptr2-ptr);
	//			result[i*cols+j] = atof(buf);
	//			ptr = ptr2 + 1;

	//			printf("%d\n",j);
	//		}
	//	}
	//}
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
	    if(j == 3)//trace for small instance type
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
