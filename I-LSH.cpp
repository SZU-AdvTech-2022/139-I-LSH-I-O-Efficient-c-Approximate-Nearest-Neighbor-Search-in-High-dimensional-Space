#include <iostream>
#include <random>
#include <array>
#include <algorithm>
#include <cmath>
//#include "B+Tree.cpp"

#define DIMENSION 100 //���ݼ�ά��d��
#define MLSHFUNCTIONS 10 //m��LSH������
#define DATAAMOUNT 100 //�ݶ���100�����ݡ�

typedef struct {
	int dn; //datanum���ݱ�š�
	double dp; //dataprojection����ͶӰֵ��
	bool visited; //�Ƿ��ѱ����ʣ�
}nppair; //(num,projection)pair�ṹ�塣

std::array<std::array<double, DIMENSION>, MLSHFUNCTIONS> LSH; //һ����Ϊm��array��������һ����Ϊd��double��array��
std::array<std::array<double, DIMENSION>, DATAAMOUNT> Data; //һ����Ϊ���������ȵ�array��������һ����Ϊd��double��array��
std::array<std::array<nppair, DATAAMOUNT>, MLSHFUNCTIONS> Projection; //һ����Ϊm��array��������һ�����������ȵ�double��array��
std::array<double, DIMENSION> query; //һ����Ϊd��double�Ͳ�ѯ������
std::array<double, MLSHFUNCTIONS> qp; //һ����Ϊm��double�Ͳ�ѯ������ͶӰ������

std::random_device rd; //rd():�������������
std::default_random_engine randomgen{ rd() }; //randomgen():������ʼ����Ϊ���������α�����������
//std::default_random_engine rng1; //Ĭ�ϳ�ʼ���ӵ�α�����������
//std::default_random_engine rng2{100}; //������ʼ����Ϊ100��α�����������
std::normal_distribution<double> p2stable(0, 1); //p2stable(randomgen):��׼��̬�ֲ��������ȡֵ����Ҫ�������������

void InitialLSHFunction() {
	std::cout << "���������ʼ��LSH����......";
	for (int i = 0; i < MLSHFUNCTIONS; i++) {
		for (int j = 0; j < DIMENSION; j++) {
			LSH[i][j] = p2stable(randomgen);
			//std::cout << LSH[i][j] << std::endl;
		}
	}
	std::cout << "����ɣ�" << std::endl;
}//��ʼ��m��dά����������LSH��

void InitialData() {
	std::cout << "���������ʼ����������......";
	for (int i = 0; i < DATAAMOUNT; i++) {
		for (int j = 0; j < DIMENSION; j++) {
			Data[i][j] = randomgen()/1000000; //�������С�ݶ�
			//std::cout << ++testsum << std::endl;
		}
	}
	std::cout << "����ɣ�" << std::endl;
}//��ʼ��DATAAMOUNT��dά��������Ϊ�������ݼ���

double DotProduct(std::array<double, DIMENSION> a, std::array<double, DIMENSION> b) {
	double sum = 0;
	for (int d = 0; d < DIMENSION; d++) {
		sum += a[d] * b[d];
	}
	return sum;
}//dά�����ڻ���

void DataProject() {
	std::cout << "�����������������LSHӳ��......";
	for (int i = 0; i < MLSHFUNCTIONS; i++) {
		for (int j = 0; j < DATAAMOUNT; j++) {
			Projection[i][j].dn = j+1;
			Projection[i][j].dp = DotProduct(LSH[i], Data[j]);
			//Projection[i][j].visited = false;
			//std::cout << Projection[i][j].dn << std::endl;
		}
	}
	std::cout << "����ɣ�" << std::endl;
}//�����������m��ͶӰά�ȵ�ͶӰֵ��

void ResetDataProjectVisited() {
	for (int i = 0; i < MLSHFUNCTIONS; i++) {
		for (int j = 0; j < DATAAMOUNT; j++) {
			Projection[i][j].visited = false;
		}
	}
}//����Projection��visited��

bool nppcompair(nppair a, nppair b) {
	return a.dp < b.dp;
}//nppair�ıȽ�����

void ProjectionSort() {
	std::cout << "�����������˳��������......";
	for (int i = 0; i < MLSHFUNCTIONS; i++) {
		sort(Projection[i].begin(), Projection[i].end(), nppcompair);
	}
	std::cout << "����ɣ�" << std::endl;
}//������Projection�������

int QuickSearch(std::array<nppair, DATAAMOUNT> npp, int head, int rear, double qp) {
	if (head == rear - 1) return ((qp - npp[head].dp) < (npp[rear].dp - qp)) ? head : rear;
	int middle = (head + rear) / 2;
	return (npp[middle].dp > qp) ? QuickSearch(npp, head, middle, qp) : QuickSearch(npp, middle, rear, qp);
}//��nppair��һ��ά�ȵ�����ͶӰֵ���ֲ��ҡ�

void GetQuery() {
	std::cout << "���ڿ�ʼ�����ѯ������" << std::endl;
	for (int d = 0; d < DIMENSION; d++) {
		std::cout << "�������ѯ�����ĵ�" << d + 1 << "��������";
		std::cin >> query[d];
	}
	std::cout << "��ѯ��������������ɣ�������Ϊ��(" << query[0];
	for (int d = 1; d < DIMENSION; d++) {
		std::cout << ", " << query[d];
	}
	std::cout << ")��" << std::endl;
}//��ȡһ����ѯ������

void QueryProject() {
	std::cout << "������ɲ�ѯ������LSHӳ��......";
	for (int i = 0; i < MLSHFUNCTIONS; i++) {
		qp[i] = DotProduct(query, LSH[i]);
		//std::cout << qp[i];
	}
	std::cout << "����ɣ�" << std::endl;
}//��ɲ�ѯ������ͶӰ��

double Distance(std::array<double, DIMENSION> a, std::array<double, DIMENSION> b) {
	double sum = 0;
	for (int d = 0; d < DIMENSION; d++) {
		sum += pow(a[d], 2) * pow(b[d], 2);
	}
	return sqrt(sum);
}//dά�������롣

std::array<int, DATAAMOUNT> cn; //��¼ÿ��������ͶӰά�����ѯ������������Ĵ�����
int IncrementalGeto() {
	double mptemp=DBL_MAX; //��¼��ͶӰά���о�����Сֵ��
	int mtemp; //��¼��ͶӰά���о�����Сֵ����ά�ȱ�š�
	int otemp;//��¼��ͶӰά���о�����Сֵ�ĵ��š�
	int ltemp, rtemp;
	int o;
	for (int m = 0; m < MLSHFUNCTIONS; m++) {
		o = QuickSearch(Projection[m], 0, DATAAMOUNT-1, qp[m]);
		if (Projection[m][o].visited == true) {
			ltemp = rtemp = o;
			while (Projection[m][ltemp].visited == true) {
				ltemp--;
				if (ltemp == -1) break;
			}
			while (Projection[m][rtemp].visited == true) {
				rtemp++;
				if (rtemp == DATAAMOUNT) break;
			}
			if ((ltemp == -1) && (rtemp == DATAAMOUNT)) {
				std::cout << "��" << m << "ά�ȵ�ͶӰֵ�ѱ���ѯ��ϡ�" << std::endl;
				continue;
			}
			else if (ltemp == -1) {
				o = rtemp;
			}
			else if (rtemp == DATAAMOUNT) {
				o = ltemp;
			}
			else {
				o = ((qp[m] - Projection[m][ltemp].dp) < (Projection[m][rtemp].dp - qp[m])) ? ltemp : rtemp;
			}
		}
		if (fabs(Projection[m][o].dp - qp[m]) < mptemp) {
			mptemp = Projection[m][o].dp - qp[m];
			mtemp = m;
			otemp = o;
		}
	}
	Projection[mtemp][otemp].visited = true;
	cn[otemp]++;
	return otemp;
}//������ѯ

int testsum = 0;

int main(void) {
	std::cout << "I-LSH��������" << std::endl;
	InitialLSHFunction();
	InitialData();
	DataProject();
	ProjectionSort();
	//for (int i = 0; i < DATAAMOUNT; i++) {
	//	std::cout << "(" << Projection[0][i].dn << "," << Projection[0][i].dp << ")" << std::endl;
	//}
	char queryagain; //�Ƿ�Ҫ������ѯ��
	int ncan; //��ѡ����������
	double alpha, beta, c; //����������
	double dmin; //���ѯ�����ľ�����Сֵ��
	int onum; //��ǰ������ѯ�������ı�š�
	//int mnum; //��ǰ������ѯ������ͶӰά�ȱ�š�
	int omin=DATAAMOUNT; //���ѯ����������С��������š�
	while (1) {
		queryagain = 'n';
		GetQuery();
		//for (int k = 0; k < DIMENSION; k++) query[k] = 1;//�������룬�����á�
		QueryProject();
		ResetDataProjectVisited();
		ncan = 0;
		std::cout << "���������ֵ,0<��<1��";
		std::cin >> alpha;
		while ((alpha <= 0) || (alpha >= 1)) {
			std::cout << "�������������ֵ��ע��0<��<1��";
			std::cin >> alpha;
		}
		std::cout << "������µ�ֵ,0<��<1��";
		std::cin >> beta;
		while ((beta <= 0) || (beta>= 1)) {
			std::cout << "����������µ�ֵ��ע��0<��<1��";
			std::cin >> beta;
		}
		std::cout << "������c��ֵ,c>1��";
		std::cin >> c;
		while (c <= 1) {
			std::cout << "����������c��ֵ��ע��c>1��";
			std::cin >> c;
		}
		dmin = DBL_MAX;
		for (int d = 0; d < DATAAMOUNT; d++) cn[d] = 0;
		while (ncan < beta * DATAAMOUNT) {
			onum = IncrementalGeto();
			if (cn[onum] == floor(alpha * MLSHFUNCTIONS)) {
				if (Distance(Data[onum], query) < dmin) {
					dmin = Distance(Data[onum], query);
					omin = onum;
				}
				ncan++;
			}
		}//ѭ�����Һ�ѡ��
		std::cout << "��ѯ������c-ANN�ǵ�" << omin + 1 << "�������������������ǣ�" << std::endl << "(" << Data[omin][0];
		for (int d = 1; d < DIMENSION; d++) {
			std::cout << ", " << Data[omin][d];
		}
		std::cout << ")" << std::endl << "�Ƿ�Ҫ������ѯ�����ǣ�������y���س���" << std::endl;
		std::cin >> queryagain;
		if (queryagain != 'y') break;
	}//ѭ����ѯ����

	std::cout << "��л����ʹ�ã�ף��������죡" << std::endl;
	return 0;
}