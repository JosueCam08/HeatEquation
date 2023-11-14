#include <iostream>
#include <cmath>
#include <time.h>

using namespace std;

int main(){

    // <><><><><> Declaracion de variables <><><><><> //

    int Nx, Ny, Nt;                     // Numero de divisiones del tiempo y puntos x y
    double kx, ky;                      // Coeficiente de difusibilidad
    double x_ini, x_fin, y_ini, y_fin;  // Inicio y fin del dominio x y
    double t_ini, t_fin;                // Variables de tiempo
    double dx, dy, dt;                  // Diferencial de x, y, t

    double **u;                         // Matriz iteracion n + 1
    double **v;                         // Matriz iteracion n
    double **aux;                       // Puntero auxiliar

    double *x;                          // Valores de puntos en eje x
    double *y;                          // Valores de puntos en eje y

    double rx, ry;                      // Constantes de la ecuacion
    double sum;                         // Verificacion de resultados

    // <><><><><> Inicializacion <><><><><> //

    // Parametros
    kx = ky = 1.0;

    x_ini = y_ini = 0.0;
    x_fin = y_fin = 1.0;
    t_ini = 0.0;
    t_fin = 0.5;

    Nt = 100000;
    Nx = Ny = 200;

    // Discretizacion
    dx = (x_fin - x_ini) / (Nx - 1);
    dy = (y_fin - y_ini) / (Ny - 1);
    dt = (t_fin - t_ini) / (Nt - 1);
    rx = kx * (dt / (dx * dx));
    ry = ky * (dt / (dy * dy));

    // Impresion
    cout << "Division del dominio [0,1] x [0, 1] en (" << Nx << "," << Ny << ").\n";
    cout << "Criterio CFL " << rx + ry << " < 1/2\n";

    // <><><><><> Reserva de memoria <><><><><> //

    // Iniciamos a contabilizar el tiempo
    t_ini = clock();

    x = new double [Nx + 1];
    y = new double [Ny + 1];
    
    u = new double* [Nx + 1];
    v = new double* [Nx + 1];

    for(int i = 1; i <= Nx; i++){
        u[i] = new double [Ny + 1];
        v[i] = new double [Ny + 1];
    }

    // <><><><><> Condiciones iniciales <><><><><> //

    // Calculo de puntos de x
    for(int i = 1; i <= Nx; i++)
        x[i] = x_ini + (i - 1) * dx;

    // Calculo de puntos de y
    for(int j = 1; j <= Ny; j++)
        y[j] = y_ini + (j - 1) * dy;

    // Condiciones iniciales interior
    for(int i = 1; i <= Nx; i++)
        for(int j = 1; j<= Ny; j++)
            v[i][j] = u[i][j] = sin(x[i] + y[j]) * sin(x[i] + y[j]);

    // Condiciones de frontera
    for(int j = 1; j <= Ny; j++)
        u[1][j] = u[Nx][j] = 1.0;
    
    for(int i = 1; i <= Nx; i++)
        u[i][1] = u[i][Ny] = 1.0;

    // <><><><><> Programa principal <><><><><> //
    
    for(int k = 1; k <= Nt; k++){
        for(int i = 2; i <= (Nx - 1); i++){
            for(int j = 2; j <= (Ny - 1); j++){
                double c = 1.0 - (2.0 * (rx + ry));
                u[i][j] = (c * v[i][j])
                        + (rx * v[i - 1][j])
                        + (rx * v[i + 1][j])
                        + (ry * v[i][j - 1])
                        + (ry * v[i][j + 1]);
            }
        }
        // Actualizacion
        if(k == 1){
            for(int j = 1; j <= Ny; j++)
                v[1][j] = v[Nx][j] = 1.0;
    
            for(int i = 1; i <= Nx; i++)
                v[i][1] = v[i][Ny] = 1.0;
        }
        swap(u, v);
    }
    swap(u, v);

    // Finalizamos el conteo de tiempo
    t_fin = clock();

    // <><><><><> Verificacion de resultados <><><><><> //

    sum = 0.0;
    for(int i = 1; i <= Nx; i++)
        for(int j = 1; j <= Ny; j++)
            sum += u[i][j];

    cout << "Suma total: " << sum << '\n';
    cout << "Tiempo: " << (t_fin - t_ini) / CLOCKS_PER_SEC << '\n';

    // <><><><><> Liberacion de memoria <><><><><> //

    delete [] x;
    delete [] y;

    for(int i = 1; i <= Nx; i++){
        delete [] u[i];
        delete [] v[i];
    }

    delete [] u;
    delete [] v;


    return 0;
}