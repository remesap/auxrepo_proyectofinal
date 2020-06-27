#include <iostream>
#include <vector>
#include <cmath>
#include <random>


//ecuación diferencia por resolver d²x/dt² + 2*(da/dt)*(a) = -Gm/a³
//posición acual y posición anterior

struct galaxy
{
  double mass; //masa de cada galaxia
  double r[3], v[3], f[3]; 
};


//constantes para la ecuacion de Friedmann
const double omega_A =0.73;
const double omega_k = 0;
const double omega_R = 0;
const double omega_M = 0.23;
const double H0 = 0.5; 
float DT = 0.01;
float TA = 0;
float TB = 5;
double TSTEPS = (TB - TA)/DT;
int n = 1; //orden de la ec diferencial
double a0 = 1;

//constantes para la ecuación de movimiento

const int N = 1000; //número de galaxias
const double G = 1; //constante de gravitación universal
const double e = 0.1; //smoothing parameter
const double L = 1.0;

//constantes para la función g(r)
const int NREGIONS = 5; //"cascarones" en los cuales dividimos el espacio
const double DR = 2.0; //grosor del cascaron
const double PI = M_PI;

//funciones
void initial_conditions(std::vector<double> & a);
void initial_conditions_galaxies(std::vector<galaxy> & galaxies);
void timestep_euler(std::vector<double> & a, double dt);
void timestep_rk4(std::vector<double> & a, double dt);
void timestep_motion(std::vector<galaxy> & galaxies, double dt);
void print_system(std::vector<galaxy> & galaxies);
double fe(std::vector<double> & a);
void forces(std::vector<galaxy> & galaxies,std::vector<double> & a, double y);
void start_timeintegration(std::vector<galaxy> & galaxies, double dt);

//función g:
double volume(double x, double dx);
int n_rdistance(double x, double dx, std::vector<galaxy> & galaxies, double i); //retorna el número de galaxias en ese radio
double g(double l, double nx, double x, double dx);
  
int main(void)
{
  //double e = 0.1/a
  std::vector<double> y(n);
  std::vector<galaxy> galaxies(N);
  std::vector<double> nr_histogram(NREGIONS, 0.0);
  initial_conditions(y);
  initial_conditions_galaxies(galaxies);
  forces(galaxies,y,y[0]);
  start_timeintegration(galaxies, DT);

  //print_system(galaxies);
  
  for(int i = 0; i < TSTEPS; ++i)
    {
      double t = TA + i*DT;
      //timestep_euler(y,DT);
      timestep_rk4(y, DT);
      timestep_motion(galaxies, DT); // actualiza r y v -> actualizar f
      forces(galaxies,y,y[0]);
      }

  for (int ii = 0; ii < NREGIONS; ++ii)
    {
      double r = 1.0 + ii*DR; int igalax = 34;
      for (int igalax = 0; igalax < N; ++igalax)  //for para correrlo por todas las galaxias (se necesita promediar) diviendo por N
	nr_histogram[ii] += 1.0*n_rdistance(r, DR, galaxies, igalax);
      //promediando: ....
      nr_histogram[ii] = nr_histogram[ii]/N;

      std::cout << r << "  " << nr_histogram[ii] << std::endl;
    }

    
    //print_system(galaxies);

  return 0;
}

void initial_conditions(std::vector<double> & a)
{
  a[0] = a0;
}

void timestep_euler(std::vector<double> & a,double dt)
{
  a[0] = a[0] + dt*fe(a);
}

void timestep_rk4(std::vector<double> & a,double dt)
{
  std::vector<double> stmp = a;
  std::vector <double> k1(N), k2(N), k3(N), k4(4);

  //compute k1
  for (int ii = 0; ii < n; ++ii)
    {
      k1[ii] = dt*fe(a);
    }

  //compute k2
  for (int ii = 0; ii < n; ++ii)
    {
      stmp[ii] = a[ii] + k1[ii]/2;
    }
  
  for (int ii = 0; ii < n; ++ii)
    {
      k2[ii] = dt*fe(stmp);
    }

  //compute k3
  for (int ii = 0; ii < n; ++ii)
    {
      stmp[ii] = a[ii] + k2[ii]/2;
    }

  for (int ii = 0; ii < n; ++ii)
    {
      k3[ii] = dt*fe(stmp);
    }

  //compute k4
  for (int ii = 0; ii < n; ++ii)
    {
      stmp[ii] = a[ii] + k3[ii];
    }

  for (int ii = 0; ii < n; ++ii)
    {
      k4[ii] = dt*fe(stmp);
    }

  //new state
  for (int ii = 0; ii < n; ++ii)
    {
      a[ii] = a[ii] + 1.0/6*(k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii]);
    }
  
}

double fe(std::vector<double> & a)
{
  return a[0]*H0*std::pow((omega_A + omega_k*(std::pow(a0,2)/std::pow(a[0],2)) + omega_M*(pow(a0,3)/pow(a[0],3)) + omega_R*(std::pow(a0,4)/std::pow(a[0],3))),1/2);
}


void initial_conditions_galaxies(std::vector<galaxy> & galaxies)
{
  std::mt19937 gen(1); //declarando el generador
  std::uniform_real_distribution<double> dis0(-L/2, L/2);
  std::uniform_real_distribution<double> dis1(-L/2, L/2);
  std::uniform_real_distribution<double> dis2(-L/2, L/2);
  
  
  
  for(int ii = 0; ii < N; ++ii)
    {
      galaxies[ii].mass = 1;

      galaxies[ii].r[0] = dis0(gen);
      galaxies[ii].r[1] = dis1(gen);
      galaxies[ii].r[2] = dis2(gen);
      
      galaxies[ii].v[0] = 0;
      galaxies[ii].v[1] = 0;
      galaxies[ii].v[2] = 0;
      
    }
  
}

void timestep_motion(std::vector<galaxy> & galaxies, double dt)
{
  for (auto & galaxia : galaxies) {
    for(int ii = 0; ii <= 2; ++ii) {
      galaxia.v[ii] += dt*galaxia.f[ii]/galaxia.mass;
      galaxia.r[ii] += galaxia.v[ii]*dt;
    }
  }
}

void start_timeintegration(std::vector<galaxy> & galaxies, double dt)
{
  for (auto & galaxia : galaxies) {
    for(int ii = 0; ii <= 2; ++ii) {
      galaxia.v[ii] = galaxia.v[ii] - dt*galaxia.f[ii]/(2*galaxia.mass);
    }
  }
}
void forces(std::vector<galaxy> & galaxies,std::vector<double> & a, double y)
{
  
  for (auto & galaxia : galaxies) {
    galaxia.f[0] = galaxia.f[1] = galaxia.f[2] = 0.0;
  }
  //Fuerza de gravedad 
  for(int i = 0; i <= N; ++i){
    for(int j = 0; j <= N; ++j){
      if(i==j) continue;
      else {
	double xij0 = galaxies[i].r[0] - galaxies[j].r[0];
	double xij1 = galaxies[i].r[1] - galaxies[j].r[1];
	double xij2 = galaxies[i].r[2] - galaxies[j].r[2];
    
    	galaxies[i].f[0] -= (G*galaxies[i].mass)/std::pow(y, 3)*( xij0 /( std::pow(xij0, 2) + std::pow(e/y, 2) ));
    	galaxies[i].f[1] -= (G*galaxies[i].mass)/std::pow(y, 3)*( xij1 /( std::pow(xij1, 2) + std::pow(e/y, 2) ));
    	galaxies[i].f[2] -= (G*galaxies[i].mass)/std::pow(y, 3)*( xij2 /( std::pow(xij2, 2) + std::pow(e/y, 2) ));
      }
    }
  }
  
  //Fuerza de amortiguamiento 
  for (auto & galaxia : galaxies)
    {   
      galaxia.f[0] = galaxia.f[0] - (fe(a)/y)*galaxia.v[0]; 
      galaxia.f[1] = galaxia.f[1] - (fe(a)/y)*galaxia.v[1]; 
      galaxia.f[2] = galaxia.f[2] - (fe(a)/y)*galaxia.v[2];
    }
}
  

void print_system(std::vector<galaxy> & galaxies)
{
  std::cout.precision(15); std::cout.setf(std::ios::fixed);

  for(int ii = 0; ii < N; ++ii)
    {
      std::cout << galaxies[ii].r[0] << "  " << galaxies[ii].r[1] << "  " <<  galaxies[ii].r[2] << "\n";
    }  
}


//para la función g:
double volume(double x, double dx)
{
  return 4/3*PI*(std::pow(x+dx,3) - std::pow(x,3)); 
}


int n_rdistance(double x, double dx, std::vector<galaxy> & galaxies, double i)
{
  int totalp = 0;
  for(int j = 0; j <= N; ++j){
    if(i == j) continue;
    else {
      double xij0 = galaxies[i].r[0] - galaxies[j].r[0];
      double xij1 = galaxies[i].r[1] - galaxies[j].r[1];
      double xij2 = galaxies[i].r[2] - galaxies[j].r[2];
      
      double d = std::sqrt(xij0*xij0 + xij1*xij1 + xij2*xij2);
      if ( (d > x - dx/2) && (d <= x + dx/2) )
	{
	  totalp = totalp + 1; 
	}
    }
  }
  return totalp;
}


double g(double l, double nx, double x, double dx)
{
  return (2*l/N*(N-1))*(nx/(4*PI*x*x*dx));
}

