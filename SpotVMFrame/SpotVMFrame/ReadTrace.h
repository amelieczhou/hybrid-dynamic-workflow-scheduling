
class ReadTrace{
public:
	int readfile(char* filename);
	int readfile(char* filename, int l);
	int readprice(float* result);
	int closefile();
private:
	FILE* fp;
};