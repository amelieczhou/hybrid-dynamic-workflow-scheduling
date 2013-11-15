#include "stdafx.h"
#include "PricingModel.h"
#include "ReadTrace.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

int PricingModel::init(){
	this->PonDemand[0] = 0.095;
	this->PonDemand[1] = 0.19;
	this->PonDemand[2] = 0.38;
	this->PonDemand[3] = 0.76;
	
	this->supportFile[0] = "Data/tp1.dat";
	this->supportFile[1] = "Data/tp2.dat";
	this->supportFile[2] = "Data/tp3.dat";
	this->supportFile[3] = "Data/tp4.dat";

	this->onDemandLag = 90;
	this->spotLag = 60;

	return 0;
}

int PricingModel::readFromFile(int fileindex, int dataindex, int storeindex){
	FILE* fp;
	if((fp=fopen(this->supportFile[fileindex], "r"))==NULL){
		printf("Error: file not found: %s\n", this->supportFile[fileindex]);
		exit(1);
    }    

//	printf("read from file %s", this->supportFile[fileindex]);
//	printf("read from data %d", dataindex);

	if(dataindex > 50) dataindex = 50; //boundary checking

	char str[1024]; // tmp buf for one line
    char buf[32]; // tmp buf for one digit
    char *ptr, *ptr2;
    
	int candp = 0;

	while(candp < 1000){		
		ptr = fgets(str, 1024, fp);
		if(ptr==NULL){	    
			// no more results, so just break this loop
			break;
		}

		for(int j=0;j<51;++j){
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
			if(j == dataindex){				
				this->configInfo[storeindex][candp] = atof(buf); 
			}
		    
			// move pointer to the next position
			ptr = ptr2 + 1; 
		}
		
		candp += 1;
	}

	if(fp != NULL)
		fclose(fp);

	return 0;
}

float PricingModel::getConfigInfo(int c, int p){
	return this->configInfo[c][p];
}

int PricingModel::getPricing(int* config, float* estimateT, float* result){ //config and estimateT should have the same length
	//reset status
	for(int i = 0; i<5; i++){
		for(int j = 0; j<1000; j++){
			this->configInfo[i][j] = 0.0;
		}
	}

	int configsize = 0;
	
	if(config != NULL && estimateT != NULL){
		//get configs

		int i = 0;
		while(config[i] != -1){
			if(config[i] >= 0){
				int fileIndex = config[i];
				int dataIndex = (estimateT[i]/180.0);//calculate every 3 minutes
				
				if(config[i+1] == -1){ //last one -- onDemand VM
					this->demandPrice = this->PonDemand[fileIndex];
				}else{//spot vm
					this->readFromFile(fileIndex, dataIndex, i);
				}
			}

			i += 1;
		}

		/*
		for(int m = 990; m<1000; m++){
			printf("%.3f\t", this->configInfo[0][m]);
		}
		printf("\n");
		*/
		
		//Search the pricing	
		configsize = i;
		int resultsize = i+1;

		for(int index = 0; index < resultsize-1; index++){
			result[index] = 0.0;
		}

		result[resultsize-1] = this->demandPrice;		

		//ondemand cost
		float currentcost = (estimateT[configsize-1]+this->onDemandLag) * (this->demandPrice) / 3600.0;

		for(int cindex = configsize -2; cindex >= 0; cindex--){
			if(config[cindex] >= 0){
			//forall bidding price
			float bp = 0.001;
			float candp = 0.0;
			
			for(int pindex = 999; pindex >=0; pindex--){
				float successp = this->configInfo[cindex][pindex]; 
				float newcost = 1.0 * ((estimateT[cindex] +this->spotLag)* bp/3600.0) + (1.0 - successp)*currentcost;//?
				//float newcost = successp * ((estimateT[cindex] +this->spotLag)* bp/3600.0) + (1.0 - successp)*currentcost;

				if(newcost < currentcost){
					candp = bp;
					currentcost = newcost;
				}
				bp += 0.001;
			}
			
			if(candp == 0.0)
				break;

			result[cindex+1] = candp;
			}
		}

		result[0] = currentcost;

		return 0;		
	}else{
		return 1;
	}
}

int PricingModel::finalize(){
	return 0;
}
int PricingModel::generateFile(int fileno){//index is the instance type
	//generate the success rate files for spot instances
	float bp = 0.001;
	float candp = 0.0;
	//first read all data in data.csv into memory
	float datatrace[143884];//35971*4
	ReadTrace* tr = new ReadTrace();
	tr->readtomem("data2.csv",datatrace);
	tr->closefile();
	delete tr;

	FILE* wFile;
	char line[4];
	wFile = fopen("newData/ffp.dat","w");
	if(wFile == NULL)
	{
		printf("cannot open ffp.dat file\n");
		exit(1);
	}
	//for each bidding price
	for(int j=0; j<760; j++) {
		//for each execution time
		float rate[200*4];
		for(int i=0; i<200*4; i++) rate[i] = 0;

		int totalsize = 1;
		for(int i=0; i<totalsize; i++){ //1000 times, calculate average
			int tracelag = rand()%35000;//0-34999, origin28800
			omp_set_num_threads(4);
			#pragma omp parallel	
			{
				#pragma omp for
				for(int index =0; index<4; index++){
					int time[] = {0,0,0,0};
					while(time[index]<200*3){
						if(datatrace[(tracelag+time[index])*4 + index]>bp){
						//	int id = omp_get_thread_num();
						//	printf("id:%d\n",id);
							rate[index*200+(int)floor((float)time[index]/3.0)] += 1; //fail at time
							break;
						}
						time[index] ++;
					}
				}
			}
		}
		for(int timeslot=0; timeslot<200; timeslot++){
			for(int index=0; index<4; index++)
				rate[index*200+timeslot] /= totalsize;
		}
		for(int i=0; i<4; i++)
			for(int k=0; k<200; k++){
				sprintf(line,"%.3f",rate[i*200+k]);
				fputs(line,wFile);
				if((i*200+k) == 799) fputs("\n",wFile);
				else fputs(",",wFile);
			}
		bp += 0.001;
	}

	fclose(wFile);
	return 0;
}