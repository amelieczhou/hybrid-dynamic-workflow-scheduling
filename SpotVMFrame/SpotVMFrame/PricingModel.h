const double priceOnDemand[] = {0.095, 0.19, 0.38, 0.76};

class PricingModel{
public:
	int init();
	int getPricing(int* config, double* estimateT, double* result);
	int finalize();
	double getConfigInfo(int c, int p);
	int readFromFile(int fileindex, int dataindex, int storeindex);
private:
	

	double PonDemand[4];
	char* supportFile[4];
	double onDemandLag;
	double spotLag;

	float configInfo[5][1000];	
	double demandPrice;
};