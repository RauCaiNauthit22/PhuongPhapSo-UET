#include<bits/stdc++.h>

using namespace std;

const int nj = 50, ni = 100, ni1 = 20, ni2 = 30, ni3 = 50;
double H = 1.0, L = 3.0, L1 = 0.5, L2 = 0.5, L3 = 2.0;
double eta, d_eta, xi, p, q, alpha;
double d_xi1 = 1.0 / (ni1 - 1),
       d_xi2 = 1.0 / (ni2 - 1),
       d_xi3 = 1.0 / (ni3 - 1);
double X[ni+1][nj+1], Y[ni+1][nj+1];

//--------Tạo biên----------//
void TaoBien()
{
    // Block1
    for(int i = 1;i <= ni1;i++) {

        xi = d_xi1 * (i-1);
        // Biên dưới Block 1 
        X[i][1] = L1 * xi ;
        Y[i][1] = X[i][1]; 
        // Biên trên Block 1
        X[i][nj] = X[i][1]; 
        Y[i][nj] = H;
    }

    // Block 2
    for(int i = ni1;i <= ni1 + ni2;i++) {
        xi = d_xi2 * (i-ni1);
        // Biên dưới Block 2
		X[i][1] = L1 + L2 * xi;
        Y[i][1] = -X[i][1] + 1.0;
        // Biên trên Block 2
		X[i][nj] = X[i][1];
        Y[i][nj] = H; 
    }

    // Block 3
	for(int i= ni1 + ni2;i <= ni1 + ni2 + ni3;i++)
	{
		xi = d_xi3 * (i - ni1 - ni2);
        // Biên dưới Block 3
		X[i][1] = L1 + L2 + L3 * xi;
		Y[i][1] = 0;
        // Biên trên Block 3
		X[i][nj] = X[i][1];
		Y[i][nj] = H;
	}

}

//---------Chia lưới đều-------//
void ChiaDeu()
{
    d_eta = 1.0 / (nj-1); 
	for (int i = 1; i <= ni; i++) {

        for (int j = 1; j < nj; j++)
        {
			eta = (double)d_eta * (j - 1);
			double dy = eta;
            X[i][j] = (X[i][nj] - X[i][1]) * dy + X[i][1];
			Y[i][j] = (Y[i][nj] - Y[i][1]) * dy + Y[i][1];
		}
    }
}

//----------Nén biên trên hoặc dưới-----------//
void NenBienDuoi() 
{
    p = 0.1; // Nén biên dưới
    //p = 1.8; // Nén biên trên
    q = 2.0;
    for(int i = 1;i <= ni;i++)
    {
        for(int j = 2;j <= nj-1;j++)
        {
            eta = (double)d_eta * (j - 1); // Tính khoảng cách đường lưới theo trục y
            double dy = p * eta + (1 - p) * (double)(1 - tanh(q *(1 - eta)) / tanh(q));

            X[i][j] = (X[i][nj] - X[i][1]) * dy + X[i][1];
            Y[i][j] = (Y[i][nj] - Y[i][1]) * dy + Y[i][1];
        }
    }
}

//------------------Nén biên trái hoặc phải------------//
void NenPhai_bl1() {
    //p = 0.1 // Nén biên trái
	p = 1.8;// Nén biên phải
    q = 2.0;
    for(int j = 1; j <= nj;j++) {

        for(int i = 2; i < ni1; i++) {
            xi = d_xi1 *(i-1); // Tinh Khoảng cách đường lưới theo trục x
            double dx = p*xi+(1-p)*(1 - tanh(q*(1-xi))/tanh(q));
            X[i][j] = (X[ni1][j] - X[1][j])*dx + X[1][j]; 
            Y[i][j] = (Y[ni1][j] - Y[1][j])*dx + Y[1][j];
        }
    }
}

//-------------------Nén hai biên trái phải------------//
void Nen2BienTraiPhai() {
    double xi0 = 0.5;
    alpha = 3; 
    for(int j = 1;j <= nj;j++) {

        for(int i = ni1 + 1; i < ni1 + ni2; i++) {
            xi = d_xi2 * (i - ni1); // Tính khoảng cách đường lưới theo trục x
            double dx = 0;
            if( xi >= 0 && xi <= xi0 ) {
                dx = L2 - xi0 * ( exp(alpha) - exp( alpha * (1 - xi / xi0) ) )/( exp(alpha) - 1 );
            }
            else if ( xi > xi0 && xi < 1) {
                dx = 1 - (1-xi0) * ( exp(alpha * (xi - xi0)/(1 - xi0)) - 1 )/(exp(alpha) - 1);
            }
            X[i][j] = (X[ni1+ni2][j] - X[ni1][j]) * dx + X[ni1][j]; 
            Y[i][j] = (Y[ni1+ni2][j] - Y[ni1][j]) * dx + Y[ni1][j];
        }
    }
}

void NenTrai_bl3() {
	p = 0.1;
    q = 2.0;
    for(int j = 1; j <= nj; j++) {

        for(int i = ni1 + ni2 + 1; i < ni1 + ni2 + ni3; i++) {
            xi = d_xi3 * (i - ni1 - ni2); // Tinh buoc nhay theo chieu x
            double dx = p * xi+(1 - p) * (1 - tanh(q * (1 - xi)) / tanh(q));
            X[i][j] = (X[ni1 + ni2 + ni3][j] - X[ni1+ni2][j]) * dx + X[ni1 + ni2][j];
            Y[i][j] = (Y[ni1 + ni2 + ni3][j] - Y[ni1+ni2][j]) * dx + Y[ni1 + ni2][j];
        }
    }
}

//----------------Làm tròn lưới-------------//
void lamtronluoi() {
    int n = 0;
    // Elip 
    while (n < 1000) { // Đường cong tỉ lệ thuận số lần lặp

        for (int j = 2; j <= nj - 1; j++) {

            for (int i = 2; i <= ni - 1; i++) {
                // Làm tròn lưới
                double alpha = pow(X[i][j + 1] - X[i][j - 1], 2) / 4 + pow(Y[i][j + 1] - Y[i][j - 1], 2) / 4;
                double beta = (X[i + 1][j] - X[i - 1][j]) * (X[i][j + 1] - X[i][j - 1]) / 4 + (Y[i][j + 1] - Y[i][j - 1]) * (Y[i + 1][j] - Y[i - 1][j]) / 4;
                double gam_mar = pow(X[i + 1][j] - X[i - 1][j], 2) / 4 + pow(Y[i + 1][j] - Y[i - 1][j], 2) / 4;
                double A = X[i + 1][j + 1] + X[i - 1][j - 1] - X[i - 1][j + 1] - X[i + 1][j - 1];
                double B = Y[i + 1][j + 1] + Y[i - 1][j - 1] - Y[i - 1][j + 1] - Y[i + 1][j - 1];
                // Thu gọn 
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

void XuatTec()
{
	ofstream outfiletec("zanhhhhh3bl.tec");
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

int main(){
    TaoBien();
    ChiaDeu();
    //NenBienDuoi();
    //NenPhai_bl1();
    //Nen2BienTraiPhai();
    //NenTrai_bl3();
    lamtronluoi();
    XuatTec();
    return 0;
}