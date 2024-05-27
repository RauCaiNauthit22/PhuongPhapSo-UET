#include <bits/stdc++.h>

using namespace std;

const int nj = 50; // Số đường lưới theo trục y
const int ni = 100; // Số đường lưới theo trục x
double H = 1.0, L = 3.0; // Độ cao chiều dài lưới
double eta, xi, p, q, alpha;
double d_xi = 1.0 / (ni - 1); // Khai báo d_xi và d_eta 
double d_eta = 1.0 / (nj - 1); 
double X[ni + 1][nj + 1], Y[ni + 1][nj + 1];

//--------Tạo biên----------//
void TaoBien() {
	for (int i = 1; i <= ni; i++) {
		xi = d_xi * (i - 1);
		X[i][1] = L * xi;
		Y[i][1] = 0;
		X[i][nj] = X[i][1];
		Y[i][nj] = H;
	}
}

//---------Chia lưới đều-------//
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

//----------Nén Biên trái hoặc biên bên phải---------//
void nenbienTraihoacPhai() {
	p = 0.1;   // nén biên trái
	//p = 1.8;   // nén biên phải
	q = 2.0;
	for (int j = 1; j <= nj; j++) {

		for (int i = 2; i < ni; i++) {
			eta = d_xi * (i - 1);
			double dx = p * eta + (1 - p) * (double)(1 - tanh(q * (1 - eta)) / tanh(q));
			X[i][j] = (X[ni][j] - X[1][j]) * dx + X[1][j];
			Y[i][j] = (Y[ni][j] - Y[1][j]) * dx + Y[1][j];
		}
	}
}

//-----------Nén biên trên hoặc biên dưới-------------//
void nenbienTrenhoacDuoi() {
	p = 0.1;     // nén biên dưới
	//p = 1.8;   // nén biên trên
	q = 2.0;
	for (int i = 1; i <= ni; i++) {

		for (int j = 2; j < nj;j++) {
			eta = d_eta * (j - 1);
			double dy = p * eta + (1 - p) * (double)(1 - tanh(q * (1 - eta)) / tanh(q));
			X[i][j] = (X[i][nj] - X[i][1]) * dy + X[i][1];
			Y[i][j] = (Y[i][nj] - Y[i][1]) * dy + Y[i][1];
		}
	}
}

//--------------Nén cả hai biên trái phải---------------//
void nencatraiphai() {
	double xi0 = 0.5; // chọn vị trí nén dần về 2 biên
	alpha = 3;
	for (int j = 1; j <= nj; j++) {

		for (int i = 1; i < ni; i++) {
			xi = d_xi * (i - 1); // Tinh khoảng cách đường lưới theo trục x
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

//--------------Nén hai biên trên dưới---------------//
void nencatrenduoi() {
	double eta0 = 0.5; // Chọn vị trí nén dần về 2 biên
	alpha = 3;
	for (int i = 1; i <= ni; i++) {

		for (int j = 1; j < nj; j++) {
			eta = d_eta * (j - 1); // Tính khoảng cách đường lưới theo trục x
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

//--------------------Nén tại x và y----------------//
void nentaiXvaY() {
	double x0 = 0.00001;
	double y0 = 0.5;
	double dy0 = 1.0 * y0 / H;
	double dx0 = 1.0 * x0 / L;
	alpha = 3;
	for (int i = 1; i <= ni; i++) {
		double dxi = d_xi * (i-1); // Vị trí điểm lưới theo chiều x

		for (int j = 1; j <= nj; j++) {
			if (dxi <= dx0) { // Với các điểm <= dx0
				double m = exp(alpha) - exp(alpha * (1 - dxi / dx0));
				double n = exp(alpha) - 1;
				double dx = m / n;
				X[i][j] = L * dx0 * dx;
			}
			else {
				double m = (1 - dx0) * (exp(alpha * (dxi - dx0) / (1 - dx0)) - 1);
				double n = exp(alpha) - 1;
				double dx = m / n;
				X[i][j] = x0 + L * dx;
			}
			double eta = d_eta * (j-1); // Vị trí điểm lưới theo chiều y
			if (eta <= dy0) { // Với các điểm <= dy0
				// Hàm biến đổi nền trên
				double m = exp(alpha) - exp(alpha * (1 - eta / dy0));
				double n = exp(alpha) - 1;
				double dy = m / n;
				Y[i][j] = H * dy * dy0;
			}
			else {
				// Hàm biến đổi nền dưới
				double m = (1 - dy0) * (exp(alpha * (eta - dy0) / (1 - dy0)) - 1);
				double n = exp(alpha) - 1;
				double dy = m / n;
				Y[i][j] = H * dy + y0;
			}
		}
	}
}

//-------------Xuất file----------------//
void XuatTec() {
	ofstream outfiletec("Luoivuong.tec");
	outfiletec << "VARIABLES = X Y " << endl;
	outfiletec << "ZONE I = " << ni << ", J = " << nj << ", F = POINT" << endl;

	for (int j = 1; j <= nj; j++)
	{
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
	//nenbienTrenhoacDuoi();
	//nencatraiphai();
	//nencatrenduoi();
	//nentaiXvaY();
	XuatTec();
	return 0;
}
