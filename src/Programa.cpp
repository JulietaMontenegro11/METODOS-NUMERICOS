#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

class SistemaLineal {
public:
    int n; // número de incógnitas
    int fixDecimales;
    vector<vector<double>> A;
    vector<double> b;

    SistemaLineal() {
        n = 0;
        fixDecimales = 5;
    }

    void ingresarSistema() {
        cout << "\n--- INGRESO DE SISTEMA ---\n";
        do {
            cout << "Ingrese numero de incognitas (2-5): ";
            cin >> n;
        } while (n < 2 || n > 5);

        A.assign(n, vector<double>(n));
        b.assign(n, 0.0);

        cout << "\nIngrese cada ecuacion en formato coeficientes x1..xn y termino independiente.\n";
        for (int i = 0; i < n; i++) {
            cout << "\nEcuacion " << i + 1 << ":\n";
            for (int j = 0; j < n; j++) {
                cout << "  Coeficiente de x" << j + 1 << ": ";
                cin >> A[i][j];
            }
            cout << "  Termino independiente: ";
            cin >> b[i];
        }

        do {
            cout << "\nIngrese numero de decimales (0-15): ";
            cin >> fixDecimales;
            if(fixDecimales < 0) fixDecimales = 0;
            if(fixDecimales > 15) {
                cout << "Se ajusta automáticamente a 15 decimales (maximo seguro para double).\n";
                fixDecimales = 15;
            }
        } while (fixDecimales < 0 || fixDecimales > 15);

        acomodarDiagonalDominante();
        cout << "\nSistema reacomodado con diagonal dominante.\n";
        mostrarSistema();
    }

    void mostrarSistema() {
        cout << "\n--- SISTEMA ACTUAL ---\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                cout << A[i][j] << "x" << j + 1 << "  ";
            cout << "= " << b[i] << "\n";
        }
    }

    void acomodarDiagonalDominante() {
        for (int i = 0; i < n; i++) {
            int mejorFila = i;
            double maxValor = abs(A[i][i]);
            for (int f = i + 1; f < n; f++) {
                if (abs(A[f][i]) > maxValor) {
                    mejorFila = f;
                    maxValor = abs(A[f][i]);
                }
            }
            if (mejorFila != i) {
                swap(A[i], A[mejorFila]);
                swap(b[i], b[mejorFila]);
            }
        }
    }

    void mostrarFormulas() {
        cout << "\n--- FORMULAS DE JACOBI ---\n";
        for (int i = 0; i < n; i++) {
            cout << "x" << i + 1 << " = ( " << b[i] << " - (";
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    cout << A[i][j] << " * x" << j + 1;
                    if (j < n - 1) cout << " + ";
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

    JacobiSolver(SistemaLineal *s) {
        S = s;
    }

    // Opción 3: resultado final
    void resolver() {
        int n = S->n;
        vector<double> x(n, 0.0);
        vector<double> x_new(n, 0.0);

        // Verificar ceros en diagonal
        for(int i=0; i<n; i++){
            if(S->A[i][i] == 0){
                cout << "Error: la matriz tiene 0 en la diagonal. Jacobi no puede continuar.\n";
                return;
            }
        }

        // Verificar diagonal dominante estricta
        bool converge = true;
        for(int i=0; i<n; i++){
            double suma = 0;
            for(int j=0; j<n; j++){
                if(i != j) suma += abs(S->A[i][j]);
            }
            if(abs(S->A[i][i]) <= suma){
                converge = false;
                break;
            }
        }
        if(!converge){
            cout << "ADVERTENCIA: La matriz no es estrictamente diagonal dominante.\n";
            cout << "Jacobi puede no converger o generar valores infinitos.\n";
            return;
        }

        tabla.clear();
        tabla.push_back(x); // Iteración 0

        int iter = 0;
        double tolerancia = pow(10, -S->fixDecimales);
        double error = tolerancia + 1;

        while (error > tolerancia) {
            iter++;
            for (int i = 0; i < n; i++) {
                double suma = 0;
                for (int j = 0; j < n; j++) {
                    if (j != i) suma += S->A[i][j] * x[j];
                }
                x_new[i] = (S->b[i] - suma) / S->A[i][i];
            }

            error = 0;
            for (int i = 0; i < n; i++) {
                double dif = abs(x_new[i] - x[i]);
                if (dif > error) error = dif;
            }

            x = x_new;
            tabla.push_back(x);
        }

        resultadoFinal = x;

        cout << fixed << setprecision(S->fixDecimales);
        cout << "\n--- RESULTADO FINAL ---\n";
        cout << "Iteraciones necesarias: " << iter << "\n";
        for (int i = 0; i < n; i++) cout << "x" << i + 1 << " = " << x[i] << "\n";
    }

    void mostrarTablaCompleta() {
        int n = S->n;
        cout << fixed << setprecision(S->fixDecimales);
        cout << "\n--- TABLA COMPLETA ---\n";
        for (int k = 0; k < tabla.size(); k++) {
            cout << "Iteracion " << k << ": ";
            for (int i = 0; i < n; i++) cout << "x" << i + 1 << "=" << tabla[k][i] << "  ";
            cout << "\n";
        }
    }

    void mostrarTablaRango(int ini, int fin) {
        int n = S->n;
        cout << fixed << setprecision(S->fixDecimales);

        if (ini < 0 || fin >= tabla.size() || ini > fin) {
            cout << "ADVERTENCIA: El rango solicitado no es válido.\n";
            cout << "Iteraciones disponibles: 0 a " << tabla.size() - 1 << "\n";
            return;
        }

        cout << "\n--- TABLA EN RANGO ---\n";
        for (int k = ini; k <= fin; k++) {
            cout << "Iteracion " << k << ": ";
            for (int i = 0; i < n; i++) cout << "x" << i + 1 << "=" << tabla[k][i] << "  ";
            cout << "\n";
        }
    }

    void guardarArchivo() {
        ofstream archivo("resultado.txt");
        if(!archivo.is_open()){
            cout << "Error al crear el archivo.\n";
            return;
        }

        archivo << fixed << setprecision(S->fixDecimales);
        archivo << "--- SISTEMA CON DIAGONAL DOMINANTE ---\n";
        for(int i=0;i<S->n;i++){
            for(int j=0;j<S->n;j++)
                archivo << S->A[i][j] << "x" << j+1 << "  ";
            archivo << "= " << S->b[i] << "\n";
        }

        archivo << "\n--- RESULTADO FINAL ---\n";
        for(int i=0;i<S->n;i++)
            archivo << "x" << i+1 << " = " << resultadoFinal[i] << "\n";

        archivo << "\n--- TABLA DE ITERACIONES ---\n";
        for(int k=0;k<tabla.size();k++){
            archivo << "Iteracion " << k << ": ";
            for(int i=0;i<S->n;i++)
                archivo << "x" << i+1 << "=" << tabla[k][i] << "  ";
            archivo << "\n";
        }

        archivo.close();
        cout << "Archivo 'resultado.txt' creado correctamente.\n";
    }
};

int main() {
    SistemaLineal S;
    JacobiSolver solver(&S);

    bool haySistema = false;
    int opcion;

    do {
        cout << "\n\n====== MENU ======\n";
        cout << "1. Ingresar sistema y reacomodar diagonal dominante\n";
        cout << "2. Mostrar formulas de Jacobi\n";
        cout << "3. Calcular con Jacobi hasta tolerancia (resultado final)\n";
        cout << "4. Mostrar tabla de iteraciones\n";
        cout << "5. Guardar resultado en archivo txt\n";
        cout << "6. Salir\n";
        cout << "Seleccione: ";
        cin >> opcion;

        switch (opcion) {
        case 1:
            S.ingresarSistema();
            haySistema = true;
            break;

        case 2:
            if (haySistema) S.mostrarFormulas();
            else cout << "Debe ingresar un sistema primero.\n";
            break;

        case 3:
            if (haySistema) solver.resolver();
            else cout << "Debe ingresar un sistema primero.\n";
            break;

        case 4:
            if (!haySistema) {
                cout << "Debe ingresar un sistema primero.\n";
                break;
            }
            int tipo;
            cout << "1. Tabla completa\n2. Tabla por rango\nSeleccione: ";
            cin >> tipo;
            if (tipo == 1) solver.mostrarTablaCompleta();
            else {
                int i, f;
                cout << "Iteracion inicial: ";
                cin >> i;
                cout << "Iteracion final: ";
                cin >> f;
                solver.mostrarTablaRango(i, f);
            }
            break;

        case 5:
            if (!haySistema) cout << "Debe ingresar un sistema primero.\n";
            else solver.guardarArchivo();
            break;

        case 6:
            cout << "Saliendo...\n";
            break;

        default:
            cout << "Opcion invalida.\n";
        }

    } while (opcion != 6);

    return 0;
}
