#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
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

    // Ingresar sistema 2-5 variables
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
            cout << "\nIngrese numero de decimales (0-20): ";
            cin >> fixDecimales;
        } while (fixDecimales < 0 || fixDecimales > 20);

        acomodarDiagonalDominante();
        cout << "\nSistema reacomodado con diagonal dominante.\n";
        mostrarSistema();
    }

    // Mostrar sistema
    void mostrarSistema() {
        cout << "\n--- SISTEMA ACTUAL ---\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                cout << A[i][j] << "x" << j + 1 << "  ";
            cout << "= " << b[i] << "\n";
        }
    }

    // Reacomodo para diagonal dominante
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

    // Mostrar fórmulas de Jacobi
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

//--------------------------------------------
//  JACOBI SOLVER
//--------------------------------------------
class JacobiSolver {
public:
    SistemaLineal *S;
    vector<vector<double>> tabla;

    JacobiSolver(SistemaLineal *s) {
        S = s;
    }

    // Método Jacobi: iterar hasta error 0
    void resolver() {
        int n = S->n;
        vector<double> x(n, 0.0);
        vector<double> x_new(n, 0.0);

        tabla.clear();
        tabla.push_back(x); // Iteración 0

        double error = 1;
        int iter = 0;

        while (error > 0) {
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

        cout << fixed << setprecision(S->fixDecimales);
        cout << "\n--- RESULTADO FINAL ---\n";
        cout << "Iteraciones necesarias: " << iter << "\n\n";

        for (int i = 0; i < n; i++)
            cout << "x" << i + 1 << " = " << x[i] << "\n";
    }

    // Mostrar tabla completa
    void mostrarTablaCompleta() {
        int n = S->n;
        cout << fixed << setprecision(S->fixDecimales);

        cout << "\n--- TABLA COMPLETA ---\n";
        for (int k = 0; k < tabla.size(); k++) {
            cout << "Iteracion " << k << ": ";
            for (int i = 0; i < n; i++)
                cout << "x" << i + 1 << "=" << tabla[k][i] << "  ";
            cout << "\n";
        }
    }

    // Tabla por rango
    void mostrarTablaRango(int ini, int fin) {
        int n = S->n;
        cout << fixed << setprecision(S->fixDecimales);

        cout << "\n--- TABLA EN RANGO ---\n";
        for (int k = ini; k <= fin && k < tabla.size(); k++) {
            cout << "Iteracion " << k << ": ";
            for (int i = 0; i < n; i++)
                cout << "x" << i + 1 << "=" << tabla[k][i] << "  ";
            cout << "\n";
        }
    }
};

//--------------------------------------------
//  PROGRAMA PRINCIPAL
//--------------------------------------------
int main() {
    SistemaLineal S;
    JacobiSolver solver(&S);

    bool haySistema = false;
    int opcion;

    do {
        cout << "\n\n====== MENU ======\n";
        cout << "1. Ingresar sistema y reacomodar diagonal dominante\n";
        cout << "2. Mostrar formulas de Jacobi\n";
        cout << "3. Calcular con Jacobi hasta error 0\n";
        cout << "4. Mostrar tabla de iteraciones\n";
        cout << "5. Salir\n";
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
            cout << "Saliendo...\n";
            break;

        default:
            cout << "Opcion invalida.\n";
        }

    } while (opcion != 5);

    return 0;
}
