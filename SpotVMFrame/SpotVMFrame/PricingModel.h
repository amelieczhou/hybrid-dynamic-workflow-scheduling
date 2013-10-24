const float priceOnDemand[] = {0.095, 0.19, 0.38, 0.76};

class PricingModel{
public:
	int init();
	int getPricing(int* config, float* estimateT, float* result);
	int finalize();
	float getConfigInfo(int c, int p);
	int readFromFile(int fileindex, int dataindex, int storeindex);
	int generateFile(int i);
private:
	

	float PonDemand[4];
	char* supportFile[4];
	float onDemandLag;
	float spotLag;

	float configInfo[5][1000];	
	float demandPrice;
};