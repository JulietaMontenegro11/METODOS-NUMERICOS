#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

using namespace std;

class SistemaLineal {
public:
    int n;
    int fixDecimales;
    vector<vector<double>> A;
    vector<double> b;
    bool diagonalDominante = false;

    SistemaLineal() {
        n = 0;
        fixDecimales = 5;
    }

    bool esDiagonalDominante(const vector<vector<double>> &mat) {
        for(int i=0;i<n;i++){
            double suma = 0;
            for(int j=0;j<n;j++){
                if(i!=j) suma += abs(mat[i][j]);
            }
            if(abs(mat[i][i]) <= suma) return false;
        }
        return true;
    }

    void acomodarDiagonalDominante() {
        vector<int> filas(n);
        for(int i=0;i<n;i++) filas[i]=i;

        vector<vector<double>> mejorA = A;
        vector<double> mejorB = b;
        bool encontrada = false;

        do {
            vector<vector<double>> tempA(n, vector<double>(n));
            vector<double> tempB(n);
            for(int i=0;i<n;i++){
                tempA[i] = A[filas[i]];
                tempB[i] = b[filas[i]];
            }
            if(esDiagonalDominante(tempA)){
                mejorA = tempA;
                mejorB = tempB;
                encontrada = true;
                break;
            }
        } while(next_permutation(filas.begin(), filas.end()));

        A = mejorA;
        b = mejorB;

        if(encontrada){
            cout << "\n                                                      Sistema reacomodado con diagonal dominante estricta.\n";
            diagonalDominante = true;
        } else {
            cout << "\n                                                      No es posible reorganizar para diagonal dominante estricta.\n"
                 << "\n                                                      Se mantiene el sistema tal como fue ingresado.\n";
            diagonalDominante = false;
        }
    }

    void ingresarSistema() {
        cout << "\n                                                                   INGRESO DE SISTEMA \n\n";
        do {
            cout <<"                                                      Ingrese numero de incognitas (2-5): ";
            cin >> n;
        } while (n < 2 || n > 5);

        A.assign(n, vector<double>(n));
        b.assign(n, 0.0);

        for (int i = 0; i < n; i++) {
            cout << "\n                                                      Ecuacion " << i+1 << ":\n";
            for (int j = 0; j < n; j++) {
                cout << "                                                      Coeficiente de x" << j+1 << ": ";
                cin >> A[i][j];
            }
            cout << "                                                      Termino independiente: ";
            cin >> b[i];
        }

        do {
            cout << "\n                                                      Ingrese numero de decimales (0-15): ";
            cin >> fixDecimales;
            if(fixDecimales < 0) fixDecimales = 0;
            if(fixDecimales > 15){
                cout << "                                                      Se ajusta a 15 decimales maximo.\n";
                fixDecimales = 15;
            }
        } while(fixDecimales <0 || fixDecimales >15);

        acomodarDiagonalDominante();
        mostrarSistema();
    }

    void mostrarSistema() {
        cout << "\n                                                                   SISTEMA ACTUAL \n";

        int anchoTerminal = 140; // ancho estimado de tu terminal
        int nColumnas = n;      // número de incógnitas
        int colAncho = 5;      // espacio para cada columna (coeficiente + xN)
        int totalAncho = nColumnas * colAncho + 5; // 5 para el "= b"
        int margenIzq = max(0, (anchoTerminal - totalAncho)/2);

        for(int i=0;i<n;i++){
            cout << string(margenIzq, ' '); // margen izquierdo para centrar la fila
            for(int j=0;j<n;j++){
                cout << setw(colAncho) << A[i][j] << "x" << j+1;
            }
            cout << " = " << b[i] << "\n";
        }
    }

    void mostrarFormulas() {
        cout << "\n                                                                   FORMULAS DE JACOBI \n";
        int anchoTerminal = 80;         // ancho estimado de la terminal
        int nColumnas = n;              // número de incógnitas
        int colAncho = 20;              // espacio reservado para cada término (coeficiente*x)
        int totalAncho = nColumnas * colAncho + 10; // +10 para paréntesis y división
        int margenIzq = max(0, (anchoTerminal - totalAncho)/2);

        for(int i=0;i<n;i++){
            cout << string(margenIzq, ' '); // margen izquierdo para centrar la fórmula
            cout << "x" << i+1 << " = ( " << b[i] << " - (";

            for(int j=0;j<n;j++){
                if(j!=i){
                    cout << A[i][j] << "*x" << j+1;
                    if(j < n-1) cout << " + ";
                }
            }
            cout << ") ) / " << A[i][i] << "\n";
        }
    }
};

class JacobiSolver {
public:
    SistemaLineal *S;
    vector<vector<double>> tabla;
    vector<double> resultadoFinal;
    bool ejecutado = false;

    JacobiSolver(SistemaLineal *s){S=s;}

    void resolver() {
        int n = S->n;
        vector<double> x(n,0.0);      // vector con valores iniciales de las incógnitas (iteración 0)
        vector<double> x_new(n,0.0);  // vector para almacenar los nuevos valores de cada iteración
        ejecutado = false;

        // Verificar que la diagonal no tenga ceros
        for(int i=0;i<n;i++){
            if(S->A[i][i]==0){
                cout << "\n                                                      Error: hay cero en diagonal. Jacobi no puede continuar.\n";
                return;
            }
        }

        // Verificar si la matriz es estrictamente diagonal dominante
        if(!S->esDiagonalDominante(S->A)){
            cout << "\n                                                      ADVERTENCIA: La matriz no es estrictamente diagonal dominante.\n"
                 << "\n                                                      Jacobi puede no converger o generar valores infinitos.\n";
        }

        ejecutado = true;
        tabla.clear();
        tabla.push_back(x); // guardamos iteración 0

        int iter=0;
        double tolerancia = pow(10,-S->fixDecimales); // tolerancia según el número de decimales
        double error = tolerancia+1;

        // Ciclo principal de Jacobi
        while(error>tolerancia){
            iter++;
            for(int i=0;i<n;i++){
                double suma=0;
                for(int j=0;j<n;j++)
                    if(j!=i) suma += S->A[i][j]*x[j];
                x_new[i]=(S->b[i]-suma)/S->A[i][i];
            }

            // calcular error máximo
            error = 0;
            for(int i=0;i<n;i++){
                double dif = abs(x_new[i]-x[i]);
                if(dif>error) error=dif;
            }

            x = x_new;
            tabla.push_back(x); // guardar resultados de la iteración
        }

        resultadoFinal = x;

        // ---- Imprimir resultados centrados ----
        int anchoTerminal = 140;       // ancho estimado de la terminal
        int nColumnas = n;            // número de incógnitas
        int colAncho = 20;            // espacio reservado para cada valor (x1, x2, ...)
        int totalAncho = nColumnas * colAncho + 15; // +15 para "Iteraciones necesarias" y "="
        int margenIzq = max(0, (anchoTerminal - totalAncho)/2);

        cout << "\n" << string(margenIzq, ' ') << "                                                                   RESULTADO FINAL \n";
        cout << string(margenIzq, ' ') << "Iteraciones necesarias: " << iter << "\n";

        for(int i=0;i<n;i++){
            cout << string(margenIzq, ' '); // aplica margen izquierdo
            cout << "x" << i+1 << " = " << x[i] << "\n";
        }
    }

    void mostrarTablaCompleta() {
        if(S->n == 0){
            cout << "\n                                                      Debe ingresar un sistema primero.\n";
            return;
        }

        if(tabla.empty()){
            cout << "\n                                                      ADVERTENCIA: Jacobi no se ha ejecutado. Los valores se muestran como ceros iniciales.\n";
        }

        int n = S->n;
        int anchoTerminal = 140;      // ancho estimado de la terminal
        int colAncho = 20;           // espacio para cada columna
        int totalAncho = (n + 1) * colAncho; // +1 para la columna de iteración
        int margenIzq = max(0, (anchoTerminal - totalAncho)/2);

        cout << "\n" << string(margenIzq, ' ') << "                                                                   TABLA DE ITERACIONES \n";

        // Cabecera
        cout << string(margenIzq, ' ') << setw(colAncho) << "Iteracion";
        for(int i=0;i<n;i++)
            cout << setw(colAncho) << ("x" + to_string(i+1));
        cout << "\n";

        // Filas de la tabla
        for(int k=0;k<tabla.size();k++){
            cout << string(margenIzq, ' ') << setw(colAncho-2) << k << ":";
            for(int i=0;i<n;i++)
                cout << setw(colAncho) << fixed << setprecision(S->fixDecimales) << tabla[k][i];
            cout << "\n";
        }
    }

    void mostrarTablaRango(int ini,int fin){
        if(S->n == 0){
            cout << "                                                      Debe ingresar un sistema primero.\n";
            return;
        }

        if(tabla.empty()){
            cout << "                                                      ADVERTENCIA: Jacobi no se ha ejecutado. Los valores se muestran como ceros iniciales.\n";
        }

        if(ini<0 || fin>=tabla.size() || ini>fin){
            cout << "                                                      ADVERTENCIA: Rango invalido. Iteraciones disponibles: 0 a " << tabla.size()-1 << "\n";
            return;
        }

        int n = S->n;
        int anchoTerminal = 140;      // ancho estimado de la terminal
        int colAncho = 20;           // espacio para cada columna
        int totalAncho = (n + 1) * colAncho; // +1 para la columna de iteración
        int margenIzq = max(0, (anchoTerminal - totalAncho)/2);

        cout << "\n" << string(margenIzq, ' ') << "                                                                   TABLA DE ITERACIONES EN RANGO\n";

        // Cabecera
        cout << string(margenIzq, ' ') << setw(colAncho) << "Iteracion";
        for(int i=0;i<n;i++)
            cout << setw(colAncho) << ("x" + to_string(i+1));
        cout << "\n";

        // Filas del rango
        for(int k=ini;k<=fin;k++){
            cout << string(margenIzq, ' ') << setw(colAncho-2) << k << ":";
            for(int i=0;i<n;i++)
                cout << setw(colAncho) << fixed << setprecision(S->fixDecimales) << tabla[k][i];
            cout << "\n";
        }
    }

    void guardarArchivo(){
        ofstream archivo("resultado.txt");
        if(!archivo.is_open()){cout<<"                                                      Error al crear archivo.\n"; return;}
        archivo << fixed << setprecision(S->fixDecimales);
        archivo << "                                                                   SISTEMA CON DIAGONAL DOMINANTE \n";

        for(int i=0;i<S->n;i++){
            for(int j=0;j<S->n;j++)
                archivo << S->A[i][j] << " x" << j+1 << "  ";
            archivo << "= " << S->b[i] << "\n";
        }

        if(ejecutado){
            archivo << "\n                                                                   RESULTADO FINAL \n\n";
            for(int i=0;i<S->n;i++)
                archivo << " x" << i+1 << " = " << resultadoFinal[i] << "\n";
            archivo << "\n                                                                   TABLA DE ITERACIONES \n\n";
            for(int k=0;k<tabla.size();k++){
                archivo << "                                                      Iteracion " << k << ": ";
                for(int i=0;i<S->n;i++)
                    archivo << "x" << i+1 << "=" << tabla[k][i] << "  ";
                archivo << "\n";
            }
        } else {
            archivo << "\n                                                      No es posible generar resultados ni tabla dado a que el método no puede converger.\n";
        }
        archivo.close();
        cout << "\n                                                      Archivo 'resultado.txt' creado correctamente.\n";
    }
};

// Función para limpiar pantalla
void limpiarPantalla(){
#ifdef _WIN32
    system("cls");
#else
    system("clear");
#endif
}

int main() {
    SistemaLineal S;
    JacobiSolver solver(&S);
    bool haySistema=false;
    int opcion;

    do {
        cout << "\n\n                                                                   MENU \n\n";
        cout << "                                                  1. Ingresar sistema y reacomodar diagonal dominante\n";
        cout << "                                                  2. Formulas de Jacobi\n";
        cout << "                                                  3. Calcular Jacobi hasta resultado final\n";
        cout << "                                                  4. Tabla de iteraciones\n";
        cout << "                                                  5. Guardar resultado\n";
        cout << "                                                  6. Salir\n";
        cout << "                                                  Seleccione: ";
        cin >> opcion;
        limpiarPantalla();

        switch(opcion){
            case 1:
                S.ingresarSistema();
                haySistema=true;
                break;
            case 2:
                if(haySistema) S.mostrarFormulas();
                else cout << "\n                                                      Debe ingresar un sistema primero.\n";
                break;
            case 3:
                if(haySistema) solver.resolver();
                else cout << "\n                                                      Debe ingresar un sistema primero.\n";
                break;
            case 4:
                if(!haySistema){cout<<"\n                                                      Debe ingresar un sistema primero.\n"; break;}
                int tipo; cout<<"\n                                                      1. Tabla completa\n\n                                                      2. Tabla por rango\n\n                                                      Seleccione: "; cin>>tipo;
                if(tipo==1) solver.mostrarTablaCompleta();
                else{
                    int i,f; cout<<"\n                                                      Iteracion inicial: "; cin>>i; cout<<"\n                                                      Iteracion final: "; cin>>f;
                    solver.mostrarTablaRango(i,f);
                }
                break;
            case 5:
                if(!haySistema) cout<<"                                                      Debe ingresar un sistema primero.\n";
                else solver.guardarArchivo();
                break;
            case 6: cout<<"                                                      Fin\n"; break;
            default: cout<<"                                                      Opcion invalida.\n";
        }

    } while(opcion!=6);

    return 0;
}




