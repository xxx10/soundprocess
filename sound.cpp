#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <cmath>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <ctype.h>
#include <algorithm>
#ifdef _EiC
#define WIN32
#endif
#define CLEAR_FACEINFO(f) f.x=0;f.y=0;f.width=0;f.height=0;f.isSpeaking=0;
#define FRAME_NUM_MAX 30000
#define V_SOUND 340
void writeToFile();

using namespace std;
typedef struct __faceinfo {
	double x;
	double y;
	double width;
	double height;
	int isSpeaking;
}FaceInfo;
typedef struct _point3{
	double x;
	double y;
	double z;
}point3;
int pLeftThan(const void* p1, const void* p2)
{
	return  (*(point3*)p1).x - (*(point3*)p2).x;
}
double distanceP3(const point3 &p1, const point3 &p2)
{
	return sqrt((p1.x - p2.x) * (p1.x - p2.x) + \
		(p1.y - p2.y) * (p1.y - p2.y) + \
		(p1.z - p2.z) * (p1.z - p2.z) );
}

static FaceInfo result[FRAME_NUM_MAX][10];
static int totalFrame = 0;
static int dataLength = 0;//�������ݳ��ȣ���λs
static int frameNum = 0;
static int people = 2;
static int fps = 0;//֡��
static int fs = 0;//����������
static point3 micLoc[4];
static char validFile[128];
static char videoFile[128];
static char audioFile[4][128];
static point3 speakerLoc[5]; //���5��
static int audioDataNum = 0;//ÿ���ļ����������ݸ���
int is_speaking(int people ,FaceInfo result[FRAME_NUM_MAX][10]); 
void init_params(char* filename);
void writeToFile();
void readFromFile();
point3 getSourceLoc(point3 t_delay, int m);//���·���λ��speakerLoc[m]

int main(int argc, char** argv )
{
	init_params(argv[1]);
	readFromFile();
	is_speaking(people, result);
	writeToFile();
	return 0;
}

void init_params(char* filename)
{
	FILE* fileparam = fopen(filename, "r");
	if (filename == NULL)
	{
		printf("cannot open file : %s", filename);
		exit(0);
	}
	int num_video, num_audio;
	point3 camLoc;
	fscanf(fileparam, "%d%d%d", &num_video, &num_audio, &people);
	fscanf(fileparam, "%d", &dataLength);
	fscanf(fileparam, "%d%d", &fps, &fs);
	totalFrame = fps * dataLength;
	audioDataNum = fs * dataLength;
	fscanf(fileparam, "%s", validFile);
	fscanf(fileparam, "%s", videoFile);
	fscanf(fileparam, "%lf%lf%lf", &camLoc.x, &camLoc.y, &camLoc.z);
	for (int i = 0; i <= 3; i++)
	{
		fscanf(fileparam, "%s", audioFile[i]);
		fscanf(fileparam, "%lf%lf%lf", &micLoc[i].x, &micLoc[i].y, &micLoc[i].z);
	}
	for (int i = 0; i <= people - 1; i++)
	{
		fscanf(fileparam, "%lf%lf%lf", &speakerLoc[i].x, &speakerLoc[i].y, &speakerLoc[i].z);
	}
	qsort(speakerLoc, people, sizeof(speakerLoc[0]), pLeftThan);//�������
	fclose(fileparam);
}

void readFromFile()
{
	FILE* result_input;
	result_input = fopen("result_o.dat", "r");
	for(int i = 0; i < totalFrame ; i++)
	{
		for(int j = 0 ; j < people ; j++)
		{
			fscanf(result_input, "%lf%lf%lf%lf%d", &result[i][j].x, &result[i][j].y, &result[i][j].width, &result[i][j].height, &result[i][j].isSpeaking);
		}
	}
}

void writeToFile()
{
	FILE* result_output;
	result_output = fopen("result.dat", "w");
	for(int i = 0; i < totalFrame ; i++)
	{
		for(int j = 0 ; j < people ; j++)
		{
			fprintf(result_output, "%f\t%f\t%f\t%f\t%d\t", result[i][j].x, result[i][j].y, result[i][j].width, result[i][j].height, result[i][j].isSpeaking);
		}
		fprintf(result_output, "\n");
	}
}

double xcorr(int a[], int b[], int Length, int k)//����غ���
{
	double val = 0;
	if (k >= Length || k <= -Length)
	{
		val = -1;
	}
	else
	{
		for (int i = 0; i <= Length - 1 - abs(k); i ++)
		{
			if (k >= 0)
			{
				val = val + b[i] * a[i + k];
			}
			else
			{
				val = val + a[i] * b[i - k];
			}
		}
	}
	return val;
}

int is_speaking(int people ,FaceInfo result[FRAME_NUM_MAX][10])
{

	FILE* audio[4];
	for (int i = 0; i <= 3; i++)
	{
		audio[i] = fopen(audioFile[i], "r");
		if (audio[i] == NULL)
		{
			printf("cannot open file : %s", audioFile);
			exit(0);
		}
	}
	int data;
	double powerAvg = 0;	//��ƽ������
	for (int i = 0; i <= audioDataNum - 1; i ++)
	{
		fscanf(audio[0], "%d", &data);
		powerAvg = powerAvg + data * data;
	}
	powerAvg = powerAvg / audioDataNum;
	rewind(audio[0]);
	int jmax;
	int cnt = 0;
	double audioDataNumPerFrame = double(audioDataNum) / totalFrame;
	int** audioData = new int*[4];
	for (int i = 0; i <= 3; i ++)
	{
		audioData[i] = new int[int(audioDataNumPerFrame) + 1];
	}
	for (int i = 0; i <= totalFrame - 1; i ++)
	{
		if(i == 200)
		{
			int a = 1;
		}
		jmax = int(audioDataNumPerFrame);
		if (jmax + cnt < audioDataNumPerFrame * (i + 1) - 1) //������ܳ��ֵ�audioDataNumPerFrame���������⡣
		{
			jmax = jmax + 1;
		}
		for (int k = 0; k <= 3; k ++)
		{
			for (int j = 0; j <= jmax - 1; j ++)
			{
				fscanf(audio[k], "%d", &audioData[k][j]);
			}	
			cnt++;
		}
		if (jmax == int(audioDataNumPerFrame))
		{
			for (int k = 0; k <= 3; k ++)
			{
				audioData[k][jmax] = 0;
			}
		}
		double corrMax[3] = {0};
		int argMaxCorr[3] = {0};
		point3 t_delay;
		double powerFrameAvg = 0;//��֡��������
		for(int j = 0; j <= 2; j++)//��3�������ļ�ѭ��
		{	
			for (int k = - int(audioDataNumPerFrame) + 1; k <= int(audioDataNumPerFrame) - 1; k++)
			{
				double corr = xcorr(audioData[0], audioData[j + 1], int(audioDataNumPerFrame) + 1, k);
				if (corr > corrMax[j])
				{
					corrMax[j] = corr;
					argMaxCorr[j] = k;
				}
			}
			powerFrameAvg = powerFrameAvg + corrMax[j] / (jmax - argMaxCorr[j]); 
		}
		t_delay.x = - double(argMaxCorr[0]) / fs;
		t_delay.y = - double(argMaxCorr[1]) / fs;
		t_delay.z = - double(argMaxCorr[2]) / fs;
		powerFrameAvg = powerFrameAvg / 3;
		double thresholdPower = powerAvg / 100;
		double thresholdDelay = 15e-4;
		point3 t_delay_idea;
		point3 precisionLoc;
		if (powerFrameAvg >= thresholdPower)
		{
			bool finish = false;
			for (int j = 0; j <= people - 1; j++)
			{
				t_delay_idea.x = (distanceP3(speakerLoc[j], micLoc[1]) - distanceP3(speakerLoc[j], micLoc[0])) / V_SOUND; 
				t_delay_idea.y = (distanceP3(speakerLoc[j], micLoc[2]) - distanceP3(speakerLoc[j], micLoc[0])) / V_SOUND; 
				t_delay_idea.z = (distanceP3(speakerLoc[j], micLoc[3]) - distanceP3(speakerLoc[j], micLoc[0])) / V_SOUND; 
				if (distanceP3(t_delay, t_delay_idea) <= thresholdDelay && !finish)
				{
					result[i][j].isSpeaking = 1;
					precisionLoc = getSourceLoc(t_delay, j);
					finish = true;
				}
				else
				{
					result[i][j].isSpeaking = 0;
				}
			}
		}
		else
		{
			for (int j = 0; j <= people - 1; j++)
			{
				result[i][j].isSpeaking = 0;
			}
		}
	}
	for (int i = 0; i <= 3; i++)
	{
		delete[] audioData[i];
		fclose(audio[i]);
	}
	delete audioData;
	for (int i = 0; i <= totalFrame - 3; i ++)	//���010��101��ͻ������
	{
		for (int j = 0; j <= people - 1; j ++)
		{
			if (result[i][j].isSpeaking == 1 && result[i + 1][j].isSpeaking == 0 && result[i + 2][j].isSpeaking == 1)
			{
				result[i + 1][j].isSpeaking = 1;
			}
			else if (result[i][j].isSpeaking == 0 && result[i + 1][j].isSpeaking == 1 && result[i + 2][j].isSpeaking == 0)
			{
				result[i + 1][j].isSpeaking = 0;
			}
		}
	}
	return 0;
}

int findMin(double A[21],int range)
{
	int min = 0;
	double Min = A[0];
	for(int cnt=1;cnt<=range-1;cnt++)
	{
		if(Min > A[cnt])
		{
			Min = A[cnt];			//��BUG�������ܴ��ڶ����Сֵ���������
			min = cnt;
		}
	}
	return min;
}/*
point3 getSourceLoc(point3 t_delay, int m)
{
	//�����
	point3 t_delay_idea;
	//t_delay_idea.x = (distanceP3(speakerLoc[j], micLoc[1]) - distanceP3(speakerLoc[j], micLoc[0])) / V_SOUND; 
	//t_delay_idea.y = (distanceP3(speakerLoc[j], micLoc[2]) - distanceP3(speakerLoc[j], micLoc[0])) / V_SOUND; 
	//t_delay_idea.z = (distanceP3(speakerLoc[j], micLoc[3]) - distanceP3(speakerLoc[j], micLoc[0])) / V_SOUND; 
	double d21 = (t_delay.y - t_delay.x)*V_SOUND;
	double d31 = (t_delay.z - t_delay.x)*V_SOUND;
	double err[21];	//�����������λ�õĶ�Ӧ�ӳ���ת���ɾ���)
	double Fin[61 * 11];
	double X[60 * 10];
	double Y[60 * 10];
	double Z[60 * 10];
	double x,y,z,d1,d2,A;
	int cnt_Fin = 0;
	point3 temp;
	for(x=-3;x<=3;x=x+0.1)
	{
		for(y=3;y<=4;y=y+0.1)
		{
			int d[21]={0};
			int cnt_z=0;
			for(z=-1;z<=1;z=z+0.1)
			{
				temp.x=x;
				temp.y=y;
				temp.z=z;
				d1=distanceP3(temp,micLoc[0]);	//���������λ�ú͵�һ����˵ľ���
				d2=distanceP3(temp,micLoc[1]);	//���������λ�ú͵ڶ�����˵ľ���
				err[cnt_z] = fabs(d2-d1-d21);
				cnt_z=cnt_z+1;
			}
			A=findMin(err,21);
			if(0)
			{
				int Count = 1;
				if(Count <= 1)
				{
					cout<<A<<endl;
					for(int cn=0;cn<=20;cn++)
					{
						cout<<err[cn]<<"  ";
					}
					cout<<endl;

				}
				Count ++;
			}
			int ZZ = -1 + 0.1 * (A - 1);
			temp.z = ZZ;
			double d_1 = distanceP3(temp,micLoc[0]);
			double d_2 = distanceP3(temp,micLoc[1]);
			double d_3 = distanceP3(temp,micLoc[2]);
			double fin = abs(d_2 - d_1 - d21) + abs(d_3 - d_1 - d31);
			Fin[cnt_Fin] = fin;
			X[cnt_Fin] = x;
			Y[cnt_Fin] = y;
			Z[cnt_Fin] = z;
			//���Դ���
			//cout<<"x  "<<x<<"  "<<"y  "<<y<<"  "<<"z  "<<z<<"   "<<cnt_Fin<<endl;

			cnt_Fin=cnt_Fin+1;
		}
	}
	int pos1 = findMin(Fin,596);
	//���Դ���
	if(0)
	{
		cout<<pos1<<endl;
		for (int cnt2=0;cnt2<=10;cnt2++)
		{
			cout<<"Z["<<cnt2+"]="<<Z[cnt2]<<"    ";
			cout<<"Fin["<<cnt2+"]="<<Fin[cnt2]<<"    ";
			cout<<endl;

		}
	}

	point3 result;
	result.x=X[pos1];
	result.y=Y[pos1];
	result.z=Z[pos1];
	if(0)
	{
		cout<<pos1<<endl;
		cout<<X[pos1]<<endl;
		cout<<Y[pos1]<<endl;
		cout<<Z[pos1]<<endl;
	}
	return result;//���޸�
}*/
point3 getSourceLoc(point3 t_delay, int m)
{
	//�����
	point3 t_delay_idea;
	//t_delay_idea.x = (distanceP3(speakerLoc[j], micLoc[1]) - distanceP3(speakerLoc[j], micLoc[0])) / V_SOUND; 
	//t_delay_idea.y = (distanceP3(speakerLoc[j], micLoc[2]) - distanceP3(speakerLoc[j], micLoc[0])) / V_SOUND; 
	//t_delay_idea.z = (distanceP3(speakerLoc[j], micLoc[3]) - distanceP3(speakerLoc[j], micLoc[0])) / V_SOUND; 
	double d21 = t_delay.x*V_SOUND;
	double d31 = t_delay.y*V_SOUND;
	double d41 = t_delay.z*V_SOUND;
	double x_ref = speakerLoc[m].x;
	double y_ref = speakerLoc[m].y;
	double z_ref = speakerLoc[m].z;

	double err[81];	//�����������λ�õĶ�Ӧ�ӳ���ת���ɾ���)
	double Fin[100*100];
	const int jj = 21;
	double X[jj*jj];
	double Y[jj*jj];
	double Z[jj*jj];
	double x,y,z,d1,d2,A;
	int cnt_Fin = 0;
	point3 temp;
	for(x=x_ref-0.5;x<=x_ref+0.5;x=x+0.05)
	{
		for(y=y_ref-0.5;y<=y_ref+0.5;y=y+0.05)
		{
			int cnt_z=0;
			for(z=z_ref-0.2;z<=z_ref+0.2;z=z+0.005)	//�������þ�ȷЩ������Ӱ�����յĸ��Ӷȣ���
			{
				temp.x=x;
				temp.y=y;
				temp.z=z;
				d1=distanceP3(temp,micLoc[0]);	//���������λ�ú͵�һ����˵ľ���
				d2=distanceP3(temp,micLoc[1]);	//���������λ�ú͵ڶ�����˵ľ���
				err[cnt_z] = fabs(d2-d1-d21);
				cnt_z=cnt_z+1;
			}
			A=findMin(err,101);
			int ZZ = z_ref-0.5 + 0.01*A;
			temp.z = ZZ;
			double d_1 = distanceP3(temp,micLoc[0]);
			double d_2 = distanceP3(temp,micLoc[1]);
			double d_3 = distanceP3(temp,micLoc[2]);
			double d_4 = distanceP3(temp,micLoc[3]);
			double fin = abs(d_2 - d_1 - d21) + abs(d_3 - d_1 - d31) + abs(d_4 - d_1 - d41);
			Fin[cnt_Fin] = fin;
			//if(cnt_Fin >= 
			X[cnt_Fin] = x;
			Y[cnt_Fin] = y;
			Z[cnt_Fin] = z;
			cnt_Fin=cnt_Fin+1;
		}
	}
	int pos1 = findMin(Fin,(jj-1)*(jj-1));
	point3 result;
	result.x=X[pos1];
	result.y=Y[pos1];
	result.z=Z[pos1];
	if(0)
	{
		cout<<pos1<<endl;
		cout<<X[pos1]<<endl;
		cout<<Y[pos1]<<endl;
		cout<<Z[pos1]<<endl;
	}
	return result;
}
