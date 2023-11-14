#include <iostream>
#include <cmath>
#include <mpi.h>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{

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

    // <><><><><> Declaracion de variables MPI <><><><><> //

    int numtasks;                       // Numero de procesos total
    int taskid;                         // ID del proceso actual

    int ndims = 2;                      // Dimesion del Grid
    int dims_vec[2];                    // Tamanio de cada dimension por entrada
    int periodicite[2];                 // Vector booleanos indicando la periodicidad por dimension
    int reorganisation;                 // Valor booleano sobre la reorganizacion del grid

    int coords[2];                      // Coordenadas de la vision de la topologia cartesiana

    int vecino[2];                      // Vector que indica el vecino proveniente y destino del proceso actual
    int direction;                      // x-direction = 0 | y-direction = 1
    int displasment;                    // Numero de shifts

    int NxG;                            // Tamanio global de divisiones sobre eje x
    int NyG;                            // Tamanio global de divisiones sobre eje y
    int NN;                             // Divisiones sobre el eje a paralelizar

    int it_ini, it_fin;                 // Variables que definen rangos de evaluacion segun el bloque

    int *index_global;               // Mapeo entre indices globales y locales
    double *xG;                         // Puntos x globales
    double *yG;                         // Puntos y globales

    int etiqueta = 888;                 // Etiqueta para comunicacion

    double suma_global;                         // Verificacion de resultados

    // <><><><><> Inicializacion <><><><><> //

    // Parametros
    kx = ky = 1.0;

    x_ini = y_ini = 0.0;
    x_fin = y_fin = 1.0;
    t_ini = 0.0;
    t_fin = 0.5;

    Nt = 100000;
    NxG = NyG = 200;

    // Discretizacion
    dx = (x_fin - x_ini) / (NxG - 1);
    dy = (y_fin - y_ini) / (NyG - 1);
    dt = (t_fin - t_ini) / (Nt - 1);
    rx = kx * (dt / (dx * dx));
    ry = ky * (dt / (dy * dy));

    // Impresion
    cout << "Division del dominio [0,1] x [0, 1] en (" << NxG << "," << NyG << ").\n";
    cout << "Criterio CFL " << rx + ry << " < 1/2\n";

    // <><><><><> Componentes MPI <><><><><> //

    // Inicializacion
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

    // Topologia cartesiana
    dims_vec[0] = numtasks; dims_vec[1] = 1;
    periodicite[0] = periodicite[1] = 0;
    reorganisation = 0;
    MPI_Comm comm2D;

    MPI_Cart_create(MPI_COMM_WORLD,
                    ndims,
                    dims_vec, 
                    periodicite,
                    reorganisation, 
                    &comm2D);
                
    // Obtenemos las coordenadas (identificador) del hilo
    MPI_Cart_coords(comm2D,
                    taskid,
                    ndims,
                    coords);

    cout << " Task ID: " << taskid << " Coordenada: (" << coords[0] << "," << coords[1] << ").\n";

    // MPI Vecinos
    vecino[0] = vecino[1] = MPI_PROC_NULL;
    direction = 0;
    displasment = 1;

    MPI_Cart_shift(comm2D,
                   direction,
                   displasment,
                   &vecino[0],
                   &vecino[1]);

    cout << " Task ID: " << taskid << " Vecinos: (" << vecino[0] << "," << vecino[1] << ").\n";

    // <><><><><> Division del dominio <><><><><> //

    // Divsion en eje X

    NN = floor(1.0 * NxG / numtasks);

    if(numtasks == 1)                           // Caso serial
        Nx = NxG;
    else 
        if(taskid == 0)                         // Proceso 1
            Nx = NN;
        else
            if(taskid == numtasks - 1)          // Proceso n - 1
                Nx = NxG - NN * taskid + 2;
            else                                // Proceso 0 < i < numtasks - 1
                Nx = NN + 2;
     
    // Division en eje Y
    Ny = NyG;

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

    index_global = new int [Nx + 1];
    xG = new double [NxG + 1];
    yG = new double [NyG + 1];

    // <><><><><> Indices locales y globales <><><><><> //
    /* Busca posicionar el intervalo correspondiente del proceso en el
    intervalo [2, numtasks - 1] */

    if(numtasks == 1)
        for(int i = 1; i <= Nx; i++)
            index_global[i] = i;
    else
        if(taskid == 0)                         // Caso proceso cero
            for(int i = 1; i <= Nx; i++)
                index_global[i] = i;
        else                                    // Caso 0 < i <= numtasks
            for(int i = 1; i <= Nx; i++)
                index_global[i] = (taskid * NN) - 1 + (i - 1);

    cout << " Task ID: " << taskid << " Nx Local: " << Nx << " Indices globales: (" << index_global[1] << "," << index_global[Nx] << ").\n";

    // <<< MPI Tipo de vectores (Comunicaciones) >>> //
    MPI_Datatype type_row;
    MPI_Type_vector(Ny, 1, 1, MPI_DOUBLE, &type_row);
    MPI_Type_commit(&type_row);

    // <><><><><> Condiciones iniciales <><><><><> //

    // Calculo de puntos globales x
    for(int i = 1; i <= NxG; i++)
        xG[i] = x_ini + (i - 1) * dx;

    // Calculo de puntos globales y
    for(int j = 1; j <= NyG; j++)
        yG[j] = y_ini + (j - 1) * dy;            

    // Calculo de puntos de x
    for(int i = 1; i <= Nx; i++)
        x[i] = xG[index_global[i]];

    // Calculo de puntos de y
    for(int j = 1; j <= Ny; j++)
        y[j] = yG[j];

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
        
        // Comunicacion
        if(numtasks > 1){
            MPI_Sendrecv(&u[2][1], 1, type_row, vecino[0], etiqueta,
                        &u[Nx][1], 1, type_row, vecino[1], etiqueta,
                        comm2D, MPI_STATUS_IGNORE);

            MPI_Sendrecv(&u[Nx - 1][1], 1, type_row, vecino[1], etiqueta,
                        &u[1][1], 1, type_row, vecino[0], etiqueta,
                        comm2D, MPI_STATUS_IGNORE);
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

    // Obtenemos rango a considerar
    if(numtasks == 1){
        it_ini = 1;
        it_fin = Nx;
    }else{
        if(taskid == 0){
            it_ini = 1;
            it_fin = Nx - 1;
        }else if(taskid == numtasks - 1){
            it_ini = 2;
            it_fin = Nx;
        }else{
            it_ini = 2;
            it_fin = Nx - 1;
        }
    }

    // Realizamos la suma
    sum = 0.0;
    for(int i = it_ini; i <= it_fin; i++)
        for(int j = 1; j <= Ny; j++)
            sum += u[i][j];

    // Realizamos la reduccion de todos los procesos
    MPI_Allreduce(&sum, &suma_global, 1, MPI_DOUBLE, MPI_SUM, comm2D);
    sum = suma_global;

    cout << " TaskID = " << taskid << " Suma total: " << sum << '\n';
    cout << "Tiempo: " << (t_fin - t_ini) / CLOCKS_PER_SEC << '\n';

    // <><><><><> Liberacion de memoria <><><><><> //

    // Liberacion de memoria
     MPI_Type_free(&type_row);

     // Finalizacion de la zona paralela
     MPI_Finalize();

    delete [] x;
    delete [] y;
    delete [] index_global;
    delete [] xG;
    delete [] yG;

    for(int i = 1; i <= Nx; i++){
        delete [] u[i];
        delete [] v[i];
    }

    delete [] u;
    delete [] v;


    return 0;
}