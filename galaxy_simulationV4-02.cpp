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
const double L = 1.0; //tamaño inicial del universo

//constantes para la función g
const double PI = std::acos(-1);
const double rmax = 15.0;
const double DR = 0.1; //diferencial en el cascarón esférico
const int NREGIONS = rmax/DR + 1;


//funciones

//funciones para la ecuación de Friedmann
void initial_conditions(double & a);
void timestep_euler(double & a, double dt);
void timestep_rk4(double & a, double dt);
double fe(double & a);

//funciones para resolver la ecuación de movimiento implementando Leap-Frog
void initial_conditions_galaxies(std::vector<galaxy> & galaxies);
void timestep_motion(std::vector<galaxy> & galaxies, double dt);
void forces(std::vector<galaxy> & galaxies, double & a);
void start_timeintegration(std::vector<galaxy> & galaxies, double dt);
//funciones para hallar g
double volume(double x, double dx);
int n_rdistance(double x,double dx,std::vector<galaxy> & galaxies, double i);
double g(double l, double nx, double x, double dx);

//función que imprime
void print_system(std::vector<galaxy> & galaxies);


int main(void)
{
  //double e = 0.1/a
  double y = 0;
  std::vector<galaxy> galaxies(N);
  std::vector<double> nr_histogram(NREGIONS, 0.0); int NTOTAL = 0.0;
  initial_conditions(y);
  initial_conditions_galaxies(galaxies);
  start_timeintegration(galaxies, DT);
  forces(galaxies, y);
    
  for(int i = 0; i < TSTEPS; ++i)
    {
      double t = TA + i*DT;
      timestep_rk4(y,DT);
      timestep_motion(galaxies, DT); // actualiza r y v -> actualizar f
      forces(galaxies, y);
    }
  
  for (int ii = 0; ii < NREGIONS; ++ii)
    {
      double r = 1.0 + ii*DR;
      for (int igalax = 0; igalax < N; ++igalax)  //for para correrlo por todas las galaxias (se necesita promediar) diviendo por N
	nr_histogram[ii] += 1.0*n_rdistance(r, DR, galaxies, igalax);
      //promediando: ....
      nr_histogram[ii] = nr_histogram[ii]/N;
      NTOTAL += nr_histogram[ii];
      double vol = volume(r, DR); 
      //std::cout << r << "  " << nr_histogram[ii] << std::endl;
      std::cout << r << "  " << g(vol, nr_histogram[ii], r, DR) << std::endl;
      //double NX = n_rdistance(r,DR,galaxies,igalax);
      //std::cout << NX << std::endl;

    }

  std::cout << NTOTAL << std::endl;
  // print_system(galaxies);
    
  return 0;
}

void initial_conditions(double & a)
{
  a = a0;
}

void timestep_euler(double & a,double dt)
{
  a = a + dt*fe(a);
}

void timestep_rk4(double & a,double dt)
{
  double stmp = a;
  double k1, k2, k3, k4;
  
  //compute k1
  k1 = dt*fe(a);

  //compute k2
  stmp = a + k1/2;
  k2 = dt*fe(stmp);
    

  //compute k3
  stmp = a + k2/2;
  k3 = dt*fe(stmp);

  //compute k4
  stmp = a + k3;
  k4 = dt*fe(stmp);

  //new state
  a = a + 1.0/6*(k1 + 2*k2 + 2*k3 + k4);
  
}

double fe(double & a)
{
  return a*H0*std::pow((omega_A + omega_k*(std::pow(a0,2)/std::pow(a,2)) + omega_M*(pow(a0,3)/pow(a,3)) + omega_R*(std::pow(a0,4)/std::pow(a,3))),1/2);
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
void forces(std::vector<galaxy> & galaxies, double & a)
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
    
    	galaxies[i].f[0] -= (G*galaxies[i].mass)/std::pow(a, 3)*( xij0 /( std::pow(xij0, 2) + std::pow(e/a, 2) ));
    	galaxies[i].f[1] -= (G*galaxies[i].mass)/std::pow(a, 3)*( xij1 /( std::pow(xij1, 2) + std::pow(e/a, 2) ));
    	galaxies[i].f[2] -= (G*galaxies[i].mass)/std::pow(a, 3)*( xij2 /( std::pow(xij2, 2) + std::pow(e/a, 2) ));
      }
    }
  }
  
  //Fuerza de amortiguamiento 
  for (auto & galaxia : galaxies)
    {   
      galaxia.f[0] = galaxia.f[0] - (fe(a)/a)*galaxia.v[0]; 
      galaxia.f[1] = galaxia.f[1] - (fe(a)/a)*galaxia.v[1]; 
      galaxia.f[2] = galaxia.f[2] - (fe(a)/a)*galaxia.v[2];
    }
}

double volume(double x, double dx)
{
  return 4/3*PI*(std::pow(x+dx,3) - std::pow(x,3)); 
}


int n_rdistance(double x,double dx, std::vector<galaxy> & galaxies, double i)
{
  int totalp = 0;
  for(int j = 0; j <= N; ++j){
    if(i==j) continue;
    else {
      double xij0 = galaxies[i].r[0] - galaxies[j].r[0];
      double xij1 = galaxies[i].r[1] - galaxies[j].r[1];
      double xij2 = galaxies[i].r[2] - galaxies[j].r[2];

      double d = std::sqrt(xij0*xij0 + xij1*xij1 + xij2*xij2);
      if ( d > x - dx && d <= x + dx )
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


void print_system(std::vector<galaxy> & galaxies)
{
  std::cout.precision(15); std::cout.setf(std::ios::fixed);

  for(int ii = 0; ii < N; ++ii)
    {
      std::cout << galaxies[ii].r[0] << "  " << galaxies[ii].r[1] << "  " <<  galaxies[ii].r[2] << "\n";
    }  
}
