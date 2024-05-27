#include <bits/stdc++.h>

using namespace std;

const int nj = 50, ni = 100;
double H = 1.0, L = 3.0;
double eta, xi, p, q, alpha;
double d_xi = 1.0 / (ni - 1);
//double d_xi1 = 1.0 / (ni1 - 1);
double d_eta = 1.0 / (nj - 1);
double X[ni + 1][nj + 1], Y[ni + 1][nj + 1];

void TaoBien() {
	// Block 1
	for (int i = 1; i <= ni; i++) {
		xi = d_xi * (i - 1);
		X[i][1] = L * xi;
		Y[i][1] = 0;
		X[i][nj] = X[i][1];
		Y[i][nj] = H / 2 + H/2 * xi;
	}


}

void chialuoideu() {
	d_eta = 1.0 / (nj - 1);
	for (int i = 1; i <= ni; i++) {
		double dx = d_xi * (i - 1);

		for (int j = 1; j <= nj; j++) {
			double dy = d_eta * (j - 1);
			X[i][j] = (X[i][nj] - X[i][1])*dx + X[i][1];
			Y[i][j] = (Y[i][nj] - Y[i][1])*dy + Y[i][1];
		}
	}
}

void nenbienTraihoacPhai() {
	p = 0.1; // Nén biên trái
	//p = 1.8; // Nén biên phải
	q = 2.0;
	for (int j = 1; j <= nj; j++) {

		for (int i = 2; i < ni; i++) {
			eta = d_xi * (i - 1);
			double dx = p * eta + (1 - p) * (double)(1 - tanh(q * (1 - eta)) / tanh(q));
			X[i][j] = (X[ni][j] - X[1][j]) * dx + X[1][j];
			Y[i][j] = (Y[ni][j] - Y[1][j]) * dx + Y[1][j]; // ?
		}
	}
}
void nencatraiphai() {
	double xi0 = 0.5; // Chọn vị trí nén dần về 2 biên
	alpha = 3;
	for (int j = 1; j <= nj; j++) {

		for (int i = 1; i < ni; i++) {
			xi = d_xi * (i - 1);
			double dx = 0;
			if (xi >= 0 && xi <= xi0) {
				dx = xi0 * (exp(alpha * xi / xi0) - 1) / (exp(alpha) - 1);
			}
			else if (xi > xi0 && xi < 1) {
				dx = 1 - (1 - xi0) * (exp(alpha * (xi - xi0) / (1 - xi0)) - 1) / (exp(alpha) - 1);
			}
			X[i][j] = (X[ni][j] - X[1][j]) * dx + X[1][j];
			Y[i][j] = (Y[ni][j] - Y[1][j]) * dx + Y[1][j];
		}
	}
}


void nenbienTrenhoacDuoi() {
	//p = 0.1; // Nén biên dưới
	p = 1.8; // Nén biên trên
	q = 2.0;
	for (int i = 1; i <= ni; i++) {

		for (int j = 2; j < nj; j++) {
			eta = d_eta * (j - 1);
			double dy = p * eta + (1 - p) * (double)(1 - tanh(q * (1 - eta)) / tanh(q));
			X[i][j] = (X[i][nj] - X[i][1]) * dy + X[i][1];
			Y[i][j] = (Y[i][nj] - Y[i][1]) * dy + Y[i][1];
		}
	}
}
void nencatrenduoi() {
	double eta0 = 0.5; // Chọn vị trí nén dần về 2 biên
	alpha = 3;
	for (int i = 1; i <= ni; i++) {

		for (int j = 1; j < nj; j++) {
			eta = d_eta * (j - 1);
			double dy = 0;
			if (eta >= 0 && eta <= eta0) {
				dy = eta0 * (exp(alpha * eta / eta0) - 1) / (exp(alpha) - 1);
			}
			else if (eta > eta0 && eta < 1) {
				dy = 1 - (1 - eta0) * (exp(alpha * (eta - eta0) / (1 - eta0)) - 1) / (exp(alpha) - 1);
			}
			X[i][j] = (X[i][nj] - X[i][1]) * dy + X[i][1];
			Y[i][j] = (Y[i][nj] - Y[i][1]) * dy + Y[i][1];
		}
	}
}


void lamtronluoi() {
	int n = 0;
	// Elip
	while (n < 100) { // Độ cong lưới tỉ lệ thuận số lần lặp
		for (int j = 2; j <= nj - 1; j++) {

			for (int i = 2; i <= ni - 1; i++) {
				// Làm tròn lưới
				double alpha = pow(X[i][j + 1] - X[i][j - 1], 2) / 4 + pow(Y[i][j + 1] - Y[i][j - 1], 2) / 4;
				double beta = (X[i + 1][j] - X[i - 1][j]) * (X[i][j + 1] - X[i][j - 1]) / 4 + (Y[i][j + 1] - Y[i][j - 1]) * (Y[i + 1][j] - Y[i - 1][j]) / 4;
				double gam_mar = pow(X[i + 1][j] - X[i - 1][j], 2) / 4 + pow(Y[i + 1][j] - Y[i - 1][j], 2) / 4;
				double A = X[i + 1][j + 1] + X[i - 1][j - 1] - X[i - 1][j + 1] - X[i + 1][j - 1];
				double B = Y[i + 1][j + 1] + Y[i - 1][j - 1] - Y[i - 1][j + 1] - Y[i + 1][j - 1];
				// Rút gọn
				double k = alpha * (X[i + 1][j] + X[i - 1][j]) + gam_mar * (X[i][j + 1] + X[i][j - 1]);
				double l = alpha + gam_mar;
				double m = alpha * (Y[i + 1][j] + Y[i - 1][j]) + gam_mar * (Y[i][j + 1] + Y[i][j - 1]);

				X[i][j] = k / (2 * l) - beta * A / (4 * l);
				Y[i][j] = m / (2 * l) - beta * B / (4 * l);
			}
		}
		n++;
	}
}
void XuatTec() {
	ofstream outfiletec("Luoithang.tec");
	outfiletec << "VARIABLES = X Y " << endl;
	outfiletec << "ZONE I = " << ni << ", J = " << nj << ", F = POINT" << endl;
	for (int j = 1; j <= nj; j++) {

		for (int i = 1; i <= ni; i++)
		{
			outfiletec << X[i][j] << "\t\t" << Y[i][j] << endl;
		}
	}
	outfiletec.close();
}
int main() {
	TaoBien();
	chialuoideu();
	nenbienTraihoacPhai();
	//nencatraiphai();
	//nenbienTrenhoacDuoi();
	//nencatrenduoi();
	//lamtronluoi();
	XuatTec();
	return 0;
}