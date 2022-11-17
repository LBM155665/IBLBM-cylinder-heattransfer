#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include <vector>   // vector containers

#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code
#define Pi 3.1415926

using namespace std;
const int q = 9;  //D2Q9模型 
const int NX = 300;   //X方向 
const int NY = 100;     //Y方向 
const double U = 0.1;     //速度 
double cr = 5;    //圆的半径 修改


int e[q][2] = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 1, 1 }, { -1, 1 }, { -1, -1 }, { 1, -1 } };   //Ck ：沿着格子迁移方向的单位矢量 
double w[q] = { 4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36 };   // Wk 
double rho[NX + 1][NY + 1], u[NX + 1][NY + 1][2], uf[NX + 1][NY + 1][2], f[NX + 1][NY + 1][q], f0[NX + 1][NY + 1][q], F[NX + 1][NY + 1][q], velocity[NX + 1][NY + 1];
double Force_x[NX + 1][NY + 1], Force_y[NX + 1][NY + 1];
double T[NX + 1][NY + 1], tg[NX + 1][NY + 1], g0[NX + 1][NY + 1][q], g[NX + 1][NY + 1][q], G[NX + 1][NY + 1][q], Q[NX + 1][NY + 1];
int i, j, k, ip, jp, t;
double c, Re, Td, dx, dy, Lx, Ly, cx, cy, a, b, dt, rho0, Tin, rhow, tau_f, tau_g, niu, temw,alpha,Cp, Pr;
double ds;

double Fd, Cd, Fl, Cl;

int IBM_stencil = 4; //插值模板，可用值: 2, 4
const int num_nodes = 100; // 每面墙上IBM节点的数量,原95
double x[num_nodes], y[num_nodes], vel_x[num_nodes], vel_y[num_nodes], tem[num_nodes], node_force_x[num_nodes], node_force_y[num_nodes], node_terper[num_nodes];
double Ub_x[num_nodes], Ub_y[num_nodes], Tb[num_nodes];

// *****************
// DECLARE FUNCTIONS
// *****************

void init();
double feq(int k, double rho, double u[2]);

double geq(int k, double T, double u[2]);

double force(int k, double u[2], double Force_x, double Force_y);

void output(int m);
//void Error();
double stencil(double); // 计算IBM插值模板

void first_forcing_step();
void unforced_velocity_interpolation_on_boundary();
void boundary_force_evalution_on_boundary();
void force_spread(); // 将节点力扩展到流体晶格
void update_of_velocity();
void evolution();
void update_particle_position();
void mac_quantities();
void output_cd();



int main()
{
	using namespace std;
	init();
	for (t = 0; t <= 80000; t++)
	{

		first_forcing_step();
		unforced_velocity_interpolation_on_boundary();//边界上的非节点速度插值
		boundary_force_evalution_on_boundary();//边界力的计算
		force_spread();// 将节点力扩展到流体晶格
		update_of_velocity();//更新速度
		evolution();  //3-4-5:计算未修正速度——计算LBM——修正速度 
		update_particle_position(); // 7.更新粒子的位置
		mac_quantities();  //8.计算宏观量 

		if (t >= 100)
		{
			if (t % 1000 == 0)
			{
				printf("t=%d,Cd=%lf,Cl=%lf\n", t, Cd, Cl);
				output_cd();
			}
			if (t % 1000 == 0)
			{
				output(t);
			}

		}
	}
	return 0;
}

void init()
{
	dx = 1.0;
	dy = 1.0;
	Lx = dx * double(NY);
	Ly = dy * double(NX);
	dt = dx;

	c = dx / dt;
	cx = 10 * cr;//原0.6
	cy = 0.5 * NY;//原0.5
	a = cx - cr;    //圆最左端点的坐标 
	b = cy;
	Fd = 0.;
	Fl = 0.;

	rho0 = 1.0;
	//Tin = 1.0;
	Re = 40;//修改，原20
	//Td = 0.2;
	Pr = 1.0;
	Cp = 1.0;
	niu = U * 2 * cr / Re;
	alpha = niu*Cp / Pr;
	tau_f = 3.0 * niu + 0.5;
	tau_g = 3.0 * alpha + 0.5;
	//tau_g = 1.5 * Td + 0.5;
	std::cout << "tau_f=" << tau_f << ",tau_g=" << tau_g << endl;     //std::cout<<"Hello world!!!"<<std::endl;是标准输出格式 

	for (int n = 0; n < num_nodes; ++n) {
		//圆形粒子边界    
		// C++中cos,sin,asin,acos这些三角函数操作的是弧度,而非角度.你需要把角度转化为弧度. 弧度=角度*Pi/180;//

		x[n] = a + cr - cr * cos((360.0 / num_nodes) * (Pi / 180.0) * n);
		y[n] = b + cr * sin((360.0 / num_nodes) * (Pi / 180.0) * n);

		node_force_x[n] = 0.0;
		node_force_y[n] = 0.0;
		node_terper[n] = 0.0;

		Ub_x[n] = 0.0;
		Ub_y[n] = 0.0;
		Tb[n] = 1.0;

	}
	ds = (2 * Pi * cr) / num_nodes;

	for (i = 0; i <= NX; i++)
	for (j = 0; j <= NY; j++)
	{
		u[i][j][0] = 0;
		u[i][j][1] = 0;
		T[i][j] = 0;
		uf[i][j][0] = 0;
		uf[i][j][1] = 0;
		tg[i][j] = 0;
		Force_x[i][j] = 0;
		Force_y[i][j] = 0;
		Q[i][j] = 0;
		rho[i][j] = rho0;
		//T[0][j] = Tin;
		u[0][j][0] = U;         //!!!!!!!!!!!!!!!!!!!!!!!!!!
		for (k = 0; k < q; k++)
		{
			F[i][j][k] = feq(k, rho[i][j], u[i][j]);
			G[i][j][k] = geq(k, T[i][j], u[i][j]);
		}
	}
}

double feq(int k, double rho, double u[2])     //计算平衡态分布函数 
{
	double eu, uv, feq;
	eu = (e[k][0] * u[0] + e[k][1] * u[1]);
	uv = (u[0] * u[0] + u[1] * u[1]);
	feq = w[k] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
	return feq;
}

double geq(int k, double T, double u[2]) //平衡态温度
{
	double eu, uv, geq;
	eu = (e[k][0] * u[0] + e[k][1] * u[1]);
	uv = (u[0] * u[0] + u[1] * u[1]);
	geq = w[k] * T * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
	return geq;
}

double force(int k, double u[2], double Force_x, double Force_y) //强迫项fk，采用高阶形式
{
	double fx, fy, eu;
	eu = (e[k][0] * u[0] + e[k][1] * u[1]);
	fx = w[k] * (1.0 - 0.5 / tau_f) * (3.0 * (e[k][0] - u[0]) + 9.0 * eu * e[k][0]) * Force_x;
	fy = w[k] * (1.0 - 0.5 / tau_f) * (3.0 * (e[k][1] - u[1]) + 9.0 * eu * e[k][1]) * Force_y;
	return (fx + fy);
}

double qk(int k, double Q, double rho) //源项qk
{
	double qk;
	qk = w[k] * (1.0 - 0.5 / tau_g) * rho * Q;
	return qk;
}

void first_forcing_step()
{

	for (i = 0; i <= NX; i++)
	for (j = 0; j <= NY; j++)
	{
		rho[i][j] = 0;
		u[i][j][0] = 0;
		u[i][j][1] = 0;
		T[i][j] = 0;
		for (k = 0; k < q; k++)
		{
			f[i][j][k] = F[i][j][k];
			g[i][j][k] = G[i][j][k];
			rho[i][j] += f[i][j][k];//a+=b--->a=a+b
			u[i][j][0] += e[k][0] * f[i][j][k];
			u[i][j][1] += e[k][1] * f[i][j][k];
			T[i][j] += g[i][j][k];
		}
		u[i][j][0] /= rho[i][j];//a/=b--->a=a/b
		u[i][j][1] /= rho[i][j];
		T[i][j] /= rho[i][j];
	}

	return;
}

void unforced_velocity_interpolation_on_boundary()
{
	int X_s, X_e, Y_s, Y_e;

	for (int n = 0; n < num_nodes; ++n) {
		// 首先重置节点速度，因为使用了'+='。
		vel_x[n] = 0;
		vel_y[n] = 0;
		tem[n] = 0;

		// 确定插值范围内的最低流体点阵节点(参见展开)。
		int x_int = (int)(x[n]);
		int y_int = (int)(y[n]);

		if (IBM_stencil == 2) {
			X_s = x_int;
			X_e = x_int + 1;
			Y_s = y_int;
			Y_e = y_int + 1;
		}
		else if (IBM_stencil == 4) {
			X_s = x_int - 1;
			X_e = x_int + 2;
			Y_s = y_int - 1;
			Y_e = y_int + 2;
		}

		// 检查所有邻近的流体节点。
		// 在两点插值的情况下，它是2x2流体节点。
		for (int X = X_s; X <= X_e; ++X) {
			for (int Y = Y_s; Y <= Y_e; ++Y) {

				// 计算对象节点与流体网格节点之间的距离。
				const double dist_x = x[n] - X;
				const double dist_y = y[n] - Y;

				// 基于距离计算x和y方向的插值权值。
				const double weight_x = stencil(dist_x);
				const double weight_y = stencil(dist_y);

				// 计算节点速度。
				vel_x[n] += (u[X][Y][0] * weight_x * weight_y);
				vel_y[n] += (u[X][Y][1] * weight_x * weight_y);

				//节点温度
				tem[n] += (T[X][Y] * weight_x * weight_y);
			}
		}
	}

	return;
}


void boundary_force_evalution_on_boundary() {

	for (int n = 0; n < num_nodes; ++n) {
		node_force_x[n] = 0.;
		node_force_y[n] = 0.;
		node_terper[n] = 0.;
	}

	for (int n = 0; n < num_nodes; ++n) {
		node_force_x[n] = 2 * rho0 * (Ub_x[n] - vel_x[n]) / dt;
		node_force_y[n] = 2 * rho0 * (Ub_y[n] - vel_y[n]) / dt;
		node_terper[n] = 2 * (Tb[n] - tem[n]) / dt;
		//node_terper[n] = (Tb[n] - tem[n]) / dt;
	}
	return;
}


void force_spread()
{

	int X_s, X_e, Y_s, Y_e;

	// 重置force，因为使用了'+='。
	for (int X = 0; X < NX; ++X) {
		for (int Y = 0; Y < NY; ++Y) {
			Force_x[X][Y] = 0.;
			Force_y[X][Y] = 0.;
			Q[X][Y] = 0.;
		}
	}

	//运行所有边界节点。
	for (int n = 0; n < num_nodes; ++n) {
		//确定插值范围内最低的流格节点。'Lowest'的意思是:它的x值和y值最小。
		// 范围内的其他流体节点具有坐标(x_int + 1, y_int), (x_int, y_int + 1), and (x_int + 1, y_int + 1).
		int x_int = (int)(x[n]);
		int y_int = (int)(y[n]);
		if (IBM_stencil == 2) {
			X_s = x_int;
			X_e = x_int + 1;
			Y_s = y_int;
			Y_e = y_int + 1;
		}
		else if (IBM_stencil == 4) {
			X_s = x_int - 1;
			X_e = x_int + 2;
			Y_s = y_int - 1;
			Y_e = y_int + 2;
		}

		// 检查所有邻近的流体节点。
		// 在两点插值的情况下，它是2x2流体节点。
		for (int X = X_s; X <= X_e; ++X) {
			for (int Y = Y_s; Y <= Y_e; ++Y) {

				// 计算对象节点与流体网格节点之间的距离。
				const double dist_x = x[n] - X;
				const double dist_y = y[n] - Y;

				// 基于距离计算x和y方向的插值权值。
				const double weight_x = stencil(dist_x);
				const double weight_y = stencil(dist_y);

				// 计算晶格的力。
				Force_x[X][Y] += (node_force_x[n] * weight_x * weight_y * ds);
				Force_y[X][Y] += (node_force_y[n] * weight_x * weight_y * ds);
				Q[X][Y] += (node_terper[n] * weight_x * weight_y * ds);

			}
		}
	}
	return;
}

void update_of_velocity() {

	for (i = 0; i <= NX; i++)
	for (j = 0; j <= NY; j++)
	{
		uf[i][j][0] = u[i][j][0] + 0.5 * dt * Force_x[i][j] / rho[i][j];
		uf[i][j][1] = u[i][j][1] + 0.5 * dt * Force_y[i][j] / rho[i][j];
		tg[i][j] = T[i][j] + 0.5 * dt * Q[i][j];
	}
	return;
}

void evolution()
{

	//演化
	for (i = 1; i < NX; i++)
	{
		for (j = 1; j < NY; j++)
		{
			for (k = 0; k < q; k++)
			{
				ip = i - e[k][0];
				jp = j - e[k][1];   //流动

				f0[i][j][k] = f[ip][jp][k];//只迁移，没有碰撞的f函数 
				g0[i][j][k] = g[ip][jp][k];

				F[i][j][k] = f0[i][j][k] + (feq(k, rho[ip][jp], uf[ip][jp]) - f0[i][j][k]) / tau_f + force(k, uf[ip][jp], Force_x[ip][jp], Force_y[ip][jp]);
				G[i][j][k] = g0[i][j][k] + (geq(k, T[ip][jp], uf[ip][jp]) - g0[i][j][k]) / tau_g + qk(k, Q[ip][jp], rho[ip][jp]);

			}
		}
	}

	//边界处理
	for (j = 1; j < NY; j++)                     //方腔左边界
	{
		rhow = (F[0][j][0] + F[0][j][2] + F[0][j][4] + 2 * (F[0][j][3] + F[0][j][6] + F[0][j][7])) / (1.0 - U);
		F[0][j][1] = F[0][j][3] + 2 * rhow * U / 3.0;
		F[0][j][5] = F[0][j][7] - 0.5 * (F[0][j][2] - F[0][j][4]) + rhow * U / 6.0;
		F[0][j][8] = F[0][j][6] + 0.5 * (F[0][j][2] - F[0][j][4]) + rhow * U / 6.0;

		//-------------------------------------
		temw = (G[0][j][0] + G[0][j][2] + G[0][j][4] + 2 * (G[0][j][3] + G[0][j][6] + G[0][j][7])) / (1.0 - U);//不确定U================
		G[0][j][1] = G[0][j][3] + 2 * temw * U / 3.0;
		G[0][j][5] = G[0][j][7] - 0.5 * (G[0][j][2] - G[0][j][4]) + temw * U / 6.0;
		G[0][j][8] = G[0][j][6] + 0.5 * (G[0][j][2] - G[0][j][4]) + temw * U / 6.0;
	}

	{     //入口角节点
		F[0][0][1] = F[0][0][3];
		F[0][0][2] = F[0][0][4];
		F[0][0][5] = F[0][0][7];
		rhow = (F[0][1][0] + F[0][1][2] + F[0][1][4] + 2 * (F[0][1][3] + F[0][1][6] + F[0][1][7])) / (1.0 - U);
		F[0][0][6] = 0.5 * (rhow - (F[0][0][0] + F[0][0][1] + F[0][0][2] + F[0][0][3] + F[0][0][4] + F[0][0][5] + F[0][0][7]));
		F[0][0][8] = 0.5 * (rhow - (F[0][0][0] + F[0][0][1] + F[0][0][2] + F[0][0][3] + F[0][0][4] + F[0][0][5] + F[0][0][7]));

		F[0][NY][4] = F[0][NY][2];
		F[0][NY][1] = F[0][NY][3];
		F[0][NY][8] = F[0][NY][6];
		rhow = (F[0][NY - 1][0] + F[0][NY - 1][2] + F[0][NY - 1][4] + 2 * (F[0][NY - 1][3] + F[0][NY - 1][6] + F[0][NY - 1][7])) / (1.0 - U);
		F[0][NY][5] = 0.5 * (rhow - (F[0][NY][0] + F[0][NY][1] + F[0][NY][2] + F[0][NY][3] + F[0][NY][4] + F[0][NY][6] + F[0][NY][8]));
		F[0][NY][7] = 0.5 * (rhow - (F[0][NY][0] + F[0][NY][1] + F[0][NY][2] + F[0][NY][3] + F[0][NY][4] + F[0][NY][6] + F[0][NY][8]));

		//---------------------------------------------
		G[0][0][1] = G[0][0][3];
		G[0][0][2] = G[0][0][4];
		G[0][0][5] = G[0][0][7];
		temw = (G[0][1][0] + G[0][1][2] + G[0][1][4] + 2 * (G[0][1][3] + G[0][1][6] + G[0][1][7])) / (1.0 - U);
		G[0][0][6] = 0.5 * (temw - (G[0][0][0] + G[0][0][1] + G[0][0][2] + G[0][0][3] + G[0][0][4] + G[0][0][5] + G[0][0][7]));
		G[0][0][8] = 0.5 * (temw - (G[0][0][0] + G[0][0][1] + G[0][0][2] + G[0][0][3] + G[0][0][4] + G[0][0][5] + G[0][0][7]));

		G[0][NY][4] = G[0][NY][2];
		G[0][NY][1] = G[0][NY][3];
		G[0][NY][8] = G[0][NY][6];
		temw = (G[0][NY - 1][0] + G[0][NY - 1][2] + G[0][NY - 1][4] + 2 * (G[0][NY - 1][3] + G[0][NY - 1][6] + G[0][NY - 1][7])) / (1.0 - U);
		G[0][NY][5] = 0.5 * (temw - (G[0][NY][0] + G[0][NY][1] + G[0][NY][2] + G[0][NY][3] + G[0][NY][4] + G[0][NY][6] + G[0][NY][8]));
		G[0][NY][7] = 0.5 * (temw - (G[0][NY][0] + G[0][NY][1] + G[0][NY][2] + G[0][NY][3] + G[0][NY][4] + G[0][NY][6] + G[0][NY][8]));

	}
	for (i = 1; i < NX; i++)              //下边界
	{
		F[i][0][5] = F[i][0][7];
		F[i][0][2] = F[i][0][4];
		F[i][0][6] = F[i][0][8];
		G[i][0][5] = G[i][0][7];
		G[i][0][2] = G[i][0][4];
		G[i][0][6] = G[i][0][8];
	}

	for (i = 1; i < NX; i++)                //上边界
	{
		F[i][NY][4] = F[i][NY][2];
		F[i][NY][7] = F[i][NY][5];
		F[i][NY][8] = F[i][NY][6];
		G[i][NY][4] = G[i][NY][2];
		G[i][NY][7] = G[i][NY][5];
		G[i][NY][8] = G[i][NY][6];
	}



	for (j = 0; j <= NY; j++)             //右边界
	for (k = 0; k < q; k++)
	{
		F[NX][j][k] = F[NX - 1][j][k];
		G[NX][j][k] = G[NX - 1][j][k];
	}


	return;
}


// *******************
// IBM计算模板
// *******************

// 计算IBM插值模板。

double stencil(double dist) {

	double stencil = 0.;

	if (IBM_stencil == 2) {
		stencil = 1. - abs(dist);
	}
	else if (IBM_stencil == 4) {
		if (abs(dist) < 1.) {
			stencil = 0.125 * (3. - 2 * abs(dist) + sqrt(1. + 4 * abs(dist) - 4 * SQ(dist)));
		}
		else if (abs(dist) < 2.) {
			stencil = 0.125 * (5. - 2 * abs(dist) - sqrt(-7. + 12 * abs(dist) - 4 * SQ(dist)));
		}
	}

	return stencil;
}

// *********************
// UPDATE NODE POSITIONS
// *********************

// 根据粒子节点的速度更新粒子节点的位置。
// 新的节点位置是它的旧位置加上它的当前速度(欧拉积分)。


void update_particle_position()
{
	Fd = 0.0; Fl = 0.0; Cd = 0.0; Cl = 0.0;
	for (int n = 0; n < num_nodes; ++n) {
		Fd += (-node_force_x[n] * ds);
		Fl += (-node_force_y[n] * ds);
	}

	Cd = Fd * 2 / (rho0 * U * U * 2 * cr);
	Cl = Fl * 2 / (rho0 * U * U * 2 * cr);
	return;
}

void mac_quantities()
{

	for (j = 0; j <= NY; j++)             //计算宏观变量velocity 
	for (i = 0; i <= NX; i++)
	{
		velocity[i][j] = sqrt(uf[i][j][0] * uf[i][j][0] + uf[i][j][1] * uf[i][j][1]);
	}
	return;
}

void output(int m)  //输出
{
	ostringstream name;
	name << "cavity_" << m << ".dat";
	ofstream out(name.str().c_str());
	out << "Title=\"LBM Lid Driven Flow\"\n" << "VARIABLES = \"X\",\"Y\",\"U\",\"V\",\"velocity\",\"temperature\"\n" << "ZONE T=\"BOX\",I=" << NX + 1 << ",J=" << NY + 1 << ".F=POINT" << endl;

	for (j = 0; j <= NY; j++)
	for (i = 0; i <= NX; i++)
	{
		out << double(i) << " "
			<< double(j) << " " << uf[i][j][0] << " " << uf[i][j][1] << " " << velocity[i][j] << " " << tg[i][j] << endl;
	}


}

void output_cd()
{
	ofstream outfile;
	outfile.open("cd_cl.txt", ofstream::app);
	outfile << t << " " << Cd << " " << Cl << endl;
	outfile.close();
}