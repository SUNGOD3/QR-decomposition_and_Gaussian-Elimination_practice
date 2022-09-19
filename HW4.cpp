#include <iostream>
#include <iomanip>
#include <cmath> 
#include <GL/glut.h>
using namespace std;
int width = 850, height = 850, typ = 7;
double x[21], y[21], y2[21], ty[21], A[21][21], A2[21][21], At[21][21], B[21][21], B2[21][21], B3[21][21], c2[21], c3[21], d[21], d2[21], d3[21], c[21];
double dis[11][3], nans[15][3][3],ba=1;
double px(double v, int deg) {
	double rt = 1;
	for (int i = 0; i < deg; ++i) {
		rt *= v;
		rt += 1.0;
	}
	return rt;
}
void gauss_elm(int n) {
	for (int i = 0; i < n - 1; i++) {
		// Partial pivoting
		double maxEntry = fabs(B[i][i]);
		int p = i;
		for (int k = i; k < n; k++)
			if (fabs(B[k][i]) > maxEntry) {
				p = k;
				maxEntry = fabs(B[k][i]);
			}
		if (p != i) {
			for (int j = i; j < n; j++)
				swap(B[p][j], B[i][j]);
			swap(d[p], d[i]);
		}
		//Forward elimination.
		for (int k = i + 1; k < n; k++) {
			if (B[k][i] == 0.0) continue;//lucky you
			double r = B[k][i] / B[i][i];
			for (int j = i; j < n; j++)
				B[k][j] = B[k][j] - r * B[i][j];
			d[k] = d[k] - r * d[i];
		}
	}
}
double  inner_product(double* a, double* b, int n) {
	int   i;
	double sum;
	sum = 0.0;
	for (i = 0; i < n; i++)
		sum += a[i] * b[i];
	return (sum);
}
void QR_reflect(int n) {
	for (int j = 0; j <= n - 2; ++j) {
		//creater vector v(cv)
		double cv[20] = { 0 }, t[20] = { 0 };
		for (int i = j; i <= n - 1; ++i)
			cv[i] = B2[i][j];
		double ttt = inner_product(cv, cv, n), as = sqrt(ttt);
		if (cv[j] >= 0)
			cv[j] += as;
		else
			cv[j] -= as;
		double vtv = inner_product(cv, cv, n);
		for (int k = j; k <= n - 1; ++k) {
			for (int i = j; i <= n - 1; ++i)
				t[i] = B2[i][k];
			double vtt = inner_product(cv, t, n);
			for (int i = j; i <= n - 1; ++i)
				B2[i][k] = B2[i][k] - 2.0 * (vtt / vtv) * cv[i];
		}
		for (int i = j; i <= n - 1; ++i)
			t[i] = d2[i];
		double vtb = inner_product(cv, d2, n);
		for (int k = j; k <= n - 1; ++k)
			d2[k] = d2[k] - 2.0 * (vtb / vtv) * cv[k];
	}
}
void QR_decomposition_reflect(int n, int m) {
	for (int j = 0; j <= n - 1; ++j) {
		//creater vector v(cv)
		double cv[20] = { 0 }, t[20] = { 0 };
		for (int i = j; i <= m - 1; ++i)
			cv[i] = A2[i][j];
		double ttt = inner_product(cv, cv, m), as = sqrt(ttt);
		if (cv[j] >= 0)
			cv[j] += as;
		else
			cv[j] -= as;
		double vtv = inner_product(cv, cv, m);
		for (int k = j; k <= n - 1; ++k) {
			for (int i = j; i <= m - 1; ++i)
				t[i] = A2[i][k];
			double vtt = inner_product(cv, t, m);
			for (int i = j; i <= m - 1; ++i)
				A2[i][k] = A2[i][k] - 2.0 * (vtt / vtv) * cv[i];
		}
		for (int i = j; i <= m - 1; ++i)
			t[i] = y2[i];
		double vtb = inner_product(cv, y2, m);
		for (int k = j; k <= m - 1; ++k)
			y2[k] = y2[k] - 2.0 * (vtb / vtv) * cv[k];
	}
}
void back_substitute(int n, int ty) {
	for (int i = n - 1; i >= 0; i--) {
		if (ty == 0)c[i] = d[i] / B[i][i];
		else if (ty == 1) c2[i] = d2[i] / B2[i][i];
		else c3[i] = y2[i] / A2[i][i];
		for (int j = i - 1; j >= 0; j--) {
			if (ty == 0)d[j] = d[j] - B[j][i] * c[i];
			else if (ty == 1) d2[j] = d2[j] - B2[j][i] * c2[i];
			else y2[j] = y2[j] - A2[j][i] * c3[i];
		}
	}
}
void init(int deg, int prit) {
	for (int i = 0; i < 21; ++i) {
		x[i] = y[i] = c[i] = c2[i] = c3[i] = y2[i] = ty[i] = d[i] = d2[i] = d3[i] = 0.0;
		for (int j = 0; j < 21; ++j) {
			A[i][j] = A2[i][j] = At[i][j] = B[i][j] = B2[i][j] = B3[i][j] = 0.0;
		}
	}
	for (int i = 0; i <= 10; ++i) {//sample points
		x[i] = 2.0 + i * 0.2;
		y2[i] = y[i] = px(x[i], deg);
		if (prit)cout << fixed << setprecision(15) << x[i] << " " << y[i] << "\n";
	}
	for (int i = 0; i <= 10; ++i) {
		for (int j = 0; j <= deg; ++j) {
			A2[i][j] = A[i][j] = (j == 0) ? 1.0 : A[i][j - 1] * x[i];
			At[j][i] = A[i][j];
			d3[j] = d2[j] = d[j] += At[j][i] * y[i];
		}
	}
	for (int i = 0; i <= deg; ++i)
		for (int j = 0; j <= deg; ++j)
			for (int k = 0; k <= 10; k++)
				B3[i][j] = B2[i][j] = B[i][j] += At[i][k] * A[k][j];
}
void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	double k = -4;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			glLineWidth(20.0);
			glBegin(GL_LINES);
			if(i==0)glColor3f(0, 1, 0);
			else if(i==1)glColor3f(1, 0, 0);
			else glColor3f(0, 0, 1);
			glVertex2d((double)(i + k) / 10.0, -0.75);
			glVertex2d((double)(i +k) / 10.0, nans[typ][i][j] * (1000.0*ba) -0.75);
			glEnd();
			++k;
		}
	}
	glFlush();
}
void keyboard(unsigned char key, int x, int y) {
	if (key >= '3' && key <= '9') {
		typ = key - '0';
	}
	else if (key == 'T' || key == 't') {
		typ = 10;
	}
	else if (key == '+') {
		ba *= 10;
	}
	else if (key == '-') {
		ba /= 10;
	}
	display();
}
void get_dis(int deg) {
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j)
			nans[deg][i][j] = 0.0;
	}
	for (int ca = 0, n = deg + 1, m = 11; ca < 3; ++ca) {
		double norm_1 = 0, norm_2 = 0, if_norm = 0;
		if (ca == 2) {
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j)
					dis[i][ca] += A[i][j] * c3[j];
				dis[i][ca] -= y[i];
				dis[i][ca] = abs(dis[i][ca]);
				norm_1 += dis[i][ca];
				norm_2 += dis[i][ca] * dis[i][ca];
				if_norm = max(if_norm, dis[i][ca]);
			}
		}
		else {
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					if (ca == 0)
						dis[i][ca] += A[i][j] * c[j];
					else if (ca == 1)
						dis[i][ca] += A[i][j] * c2[j];
				}
				dis[i][ca] -= y[i];
				dis[i][ca] = abs(dis[i][ca]);
				norm_1 += dis[i][ca];
				norm_2 += dis[i][ca] * dis[i][ca];
				if_norm = max(if_norm, dis[i][ca]);
			}
		}
		nans[deg][ca][0] = norm_1;
		nans[deg][ca][1] = sqrt(norm_2);
		nans[deg][ca][2] = if_norm;
	}
	cout << "     1-norm" << setw(23) << "2-norm   " << setw(23) << "infiniti-norm\n";
	for (int i = 0; i < 3; ++i) {
		cout << i + 1 << " :";
		for (int j = 0; j < 3; ++j)
			cout << nans[deg][i][j] << " ";
		cout << '\n';
	}
}
int main(int argc, char** argv) {
	cout << "Q1:\n" << setw(5) << "x" << setw(18) << "y" << "\n";
	init(7, 1);
	gauss_elm(8);
	back_substitute(8, 0);
	cout << "\nQ2:\nGauss Eliminating :\n";
	for (int i = 0; i <= 7; ++i)
		cout << "c" << i << ": " << c[i] << "\n";
	cout << "\n";
	QR_reflect(8);//A = B2 , b = d2
	// Solve the upper triangular system by using backward substitution.
	back_substitute(8, 1);
	cout << "QR Eliminating :\n";
	for (int i = 0; i <= 7; ++i)
		cout << "c" << i << ": " << c2[i] << "\n";
	QR_decomposition_reflect(8, 11);
	back_substitute(8, 2);
	cout << "\nQR Eliminating (Original system) :\n";
	for (int i = 0; i <= 7; ++i)
		cout << "c" << i << ": " << c3[i] << "\n";
	cout << "\nQ3:\n";
	get_dis(7);
	cout << "\nQ4:\n";
	for (int i = 3; i <= 10; ++i) {
		if (i != 7) {
			cout << "degree = " << i << "\n";
			init(i, 0);
			gauss_elm(i + 1);
			back_substitute(i + 1, 0);
			//cout << "\nQ2:\nGauss Eliminating :\n";
			//for (int j = 0; j <= i; ++j) 
			//	cout << "c" << j << ": " << c[j] << "\n";
			QR_reflect(i + 1);//A = B2 , b = d2
			// Solve the upper triangular system by using backward substitution.
			back_substitute(i + 1, 1);
			//cout << "\nQR Eliminating :\n";
			//for (int j = 0; j <= i; ++j) 
			//	cout << "c" << j << ": " << c2[j] << "\n";
			QR_decomposition_reflect(i + 1, 11);
			back_substitute(i + 1, 2);
			//cout << "\nQR Eliminating (Original system) :\n";
			//for (int j = 0; j <= i; ++j) 
			//	cout << "c" << j << ": " << c3[j] << "\n";
			get_dis(i);
		}
	}
	cout << "\nQ5:\n在任何情況下，直接使用QR分解原始矩陣都會是最準的；\n在degree小的情況下，對於處理過的matrix，高斯分解會比QR分解準，除了degree 8\n";

	glutInit(&argc, argv);
	glutInitWindowPosition(400, 500);	//窗口初始位置
	glutInitWindowSize(width, height);	//窗口初始大小
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);	//設定模式

	glutCreateWindow("window");
	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT);
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	display();
	glutMainLoop();
}
