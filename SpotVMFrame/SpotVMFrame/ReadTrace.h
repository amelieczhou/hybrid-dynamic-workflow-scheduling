
class ReadTrace{
public:
	int readfile(char* filename);
	int readfile(char* filename, int l);
	int readprice(float* result);
	int readtomem(char* filename, float* result);
	int closefile();

	//istream implementation
	int isreadfile(char* filename, int l);
private:
	FILE* fp;
};