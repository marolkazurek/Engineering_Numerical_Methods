import math

nx = 400
ny = 90
i1 = 200
i2 = 210
j1 = 50
delta = 0.01
sigma = 10 * delta
xA = 0.45
yA = 0.45
D=0.
IT_MAX = 20000

u0 = [0 for i in range (nx + 1)]
u1 = [0 for i in range (nx + 1)]
Vx = [0 for i in range (nx + 1)]
Vy = [0 for i in range (nx + 1)]
psi= [0 for i in range (nx + 1)]

for i in range(nx):
	u0[i] = [0 for i in range(ny + 1)]
	u1[i] = [0 for i in range(ny + 1)]
	Vx[i] = [0 for i in range(ny + 1)]
	Vy[i] = [0 for i in range(ny + 1)]
	psi[i]= [0 for i in range(ny + 1)]
		
Vx[0][0] = 0.
Vy[0][0] = 0.
	

x, y=0,0
ps=0.
while (fscanf(file, "%d %d %lf", &x, &y, &ps) == 3) {
		psi[x][y] = ps;
}


for (int i = 1; i < n_x; i++) {
	for (int j = 1; j < n_y; j++) {
		if (i >= i_1 && i <= i_2 && j <= j_1) {}
		else {
			Vx[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2 * dt);
			Vy[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * dt);
		}
	}
}
for (int i = 1; i < n_x; i++) {
		Vx[i][0] = Vy[i][n_y] = 0.;
}
for (int j = 0; j <= n_y; j++) {
	Vx[0][j] = Vx[1][j];
	Vx[n_x][j] = Vx[n_x - 1][j];
}

for (int i = 0; i <= n_x; i++) {
		for (int j = 0; j <= n_y; j++) {
			f8<<i<<"\t"<<j<<"\t"<<Vx[i][j]<<std::endl;
			f9<<i<<"\t"<<j<<"\t"<<Vy[i][j]<<std::endl;
		}
		f8<<std::endl;
		f9<<std::endl;
}

Vmax = 0.
V=0;


for (int i = 0; i <= n_x; i++) {
	for (int j = 0; j <= n_y; j++) {
		V = sqrt(pow(Vx[i][j], 2) + pow(Vy[i][j], 2));
		if (V > Vmax) {
			Vmax = V;
		}
	}
}
	double tdt = dt / (4 * Vmax);


	for (int i = 0; i <= n_x; i++) {
		for (int j = 0; j <= n_y; j++) {
			u0[i][j] = exp(-(pow(dt * i - x_A, 2) + pow(dt * j - y_A, 2)) / 2. / pow(sigma, 2)) / (2 * M_PI * pow(sigma, 2));
			u1[i][j] = u0[i][j];
		}
	}

	for (int it = 1; it <= IT_MAX; it++) {

		double c = 0.;
		double xsr = 0.;
		for (int k = 0; k < 20; k++) {
			for (int i = 0; i <= n_x; i++) {
				for (int j = 1; j < n_y; j++) {
					if (i >= i_1 && i <= i_2 && j <= j_1) {}
					else{

						int inext,iprev;
						inext=i+1;
						iprev=i-1;
						if(i==n_x)
              inext=0;
						if(i==0)
              iprev=n_x;
						u1[i][j]=(1.0 / (1.0 + ((2.0 * D * tdt) / dt*dt)))*(u0[i][j]-tdt/2*Vx[i][j]*((u0[inext][j]-u0[iprev][j])/2/dt+(u1[inext][j]-u1[iprev][j])/2/dt )-tdt/2*Vy[i][j]*((u0[i][j+1]-u0[i][j-1])/2/dt+(u1[i][j+1]-u1[i][j-1])/2/dt )+tdt*D/2*((u0[inext][j]+u0[iprev][j]+u0[i][j+1]+u0[i][j-1]-4*u0[i][j])/dt/dt+(u1[inext][j]+u1[iprev][j]+u1[i][j+1]+u1[i][j-1])/dt/dt ));
					}
				}
			}
		}
		for (int i = 0; i <= n_x; i++) {
			for (int j = 0; j <= n_y; j++) {
				u0[i][j] = u1[i][j];
				c += u0[i][j] * dt*dt;
				xsr += i * dt * u0[i][j] * dt*dt;
			}
		}
		f1<<it<<"\t"<<c<<std::endl;
		f2<<it<<"\t"<<xsr<<std::endl;

		if(D != 0){

			if(it == 1000){
				for (int i = 0; i <= n_x; i++) {
					for (int j = 0; j <= n_y; j++) {
						f3<<i<<"\t"<<j<<"\t"<<u1[i][j]<<std::endl;
					}
					f3<<std::endl;
				}
			}
			if(it == 2000){
				for (int i = 0; i <= n_x; i++) {
					for (int j = 0; j <= n_y; j++) {
						f4<<i<<"\t"<<j<<"\t"<<u1[i][j]<<std::endl;
					}
					f4<<std::endl;
				}
			}
			if (it == 3000){
				for (int i = 0; i <= n_x; i++) {
					for (int j = 0; j <= n_y; j++) {
						f5<<i<<"\t"<<j<<"\t"<<u1[i][j]<<std::endl;
					}
					f5<<std::endl;
				}
			}
			if (it == 4000){
				for (int i = 0; i <= n_x; i++) {
					for (int j = 0; j <= n_y; j++) {
						f6<<i<<"\t"<<j<<"\t"<<u1[i][j]<<std::endl;
					}
					f6<<std::endl;
				}
			}
			if (it == 5000){
				for (int i = 0; i <= n_x; i++) {
					for (int j = 0; j <= n_y; j++) {
						f7<<i<<"\t"<<j<<"\t"<<u1[i][j]<<std::endl;
					}
					f7<<std::endl;
				}
			}
		}

	}





