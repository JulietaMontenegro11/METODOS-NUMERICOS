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

        sort(filas.begin(), filas.end());
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

        int anchoTerminal = 140;
        int nColumnas = n;
        int colAncho = 5;
        int totalAncho = nColumnas * colAncho + 5;
        int margenIzq = max(0, (anchoTerminal - totalAncho)/2);

        for(int i=0;i<n;i++){
            cout << string(margenIzq, ' ');
            for(int j=0;j<n;j++){
                cout << setw(colAncho) << A[i][j] << "x" << j+1;
            }
            cout << " = " << b[i] << "\n";
        }
    }

    void mostrarFormulas() {
        cout << "\n                                                                FORMULAS DE JACOBI \n\n\n";
        int anchoTerminal = 125;
        int nColumnas = n;
        int colAncho = 1;
        int totalAncho = nColumnas * colAncho + 1;
        int margenIzq = max(0, (anchoTerminal - totalAncho)/2);

        for(int i=0;i<n;i++){
            cout << string(margenIzq, ' ');
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
    vector<vector<double>> tabla;           // filas: iteraciones, columnas: variables
    vector<vector<double>> errores;         // errores relativos entre iteraciones consecutivas (misma dimensión que tabla)
    vector<vector<double>> errorFinal;      // error relativo respecto a la última iteración (misma dimensión que tabla)
    vector<double> resultadoFinal;
    bool ejecutado = false;
    bool metodoConvergente = false;

    JacobiSolver(SistemaLineal *s){S=s;}

    void resolver() {
        int n = S->n;
        vector<double> x(n,0.0);
        vector<double> x_new(n,0.0);
        metodoConvergente = false;

        // Verificar diagonal no cero
        for(int i=0;i<n;i++){
            if(S->A[i][i]==0){
                cout << "\n                                                      Error: hay cero en diagonal. Jacobi no puede continuar.\n";
                return;
            }
        }

        // Advertencia si no es diagonal dominante
        if(!S->esDiagonalDominante(S->A)){
            cout << "\n                                                      ADVERTENCIA: La matriz no es estrictamente diagonal dominante.\n"
                 << "\n                                                      Jacobi puede no converger o generar valores infinitos.\n";
        }

        metodoConvergente = true;

        // Inicializar tablas
        tabla.clear();
        errores.clear();
        errorFinal.clear();

        tabla.push_back(x);                         // iteración 0 (valores iniciales)
        errores.push_back(vector<double>(n,0.0));  // errores iter 0 = 0

        int iter=0;
        double tolerancia = pow(10,-S->fixDecimales);
        double error = tolerancia+1;
        const int MAX_ITER = 200000; // límite seguro

        while(error>tolerancia && iter < MAX_ITER){
            iter++;
            for(int i=0;i<n;i++){
                double suma=0;
                for(int j=0;j<n;j++)
                    if(j!=i) suma += S->A[i][j]*x[j];
                x_new[i]=(S->b[i]-suma)/S->A[i][i];
            }

            // calcular error absoluto máximo entre iteraciones
            error = 0;
            for(int i=0;i<n;i++){
                double dif = abs(x_new[i]-x[i]);
                if(dif>error) error=dif;
            }

            // calcular errores relativos por componente usando la fórmula |(x_new - x_old)/x_new|
            vector<double> errIter(n, 0.0);
            for(int i=0;i<n;i++){
                if(abs(x_new[i]) > 1e-15)
                    errIter[i] = fabs((x_new[i] - x[i]) / x_new[i]);
                else
                    errIter[i] = 0.0;
            }
            errores.push_back(errIter);

            // actualizar x y guardar iteración
            x = x_new;
            tabla.push_back(x);
        }

        resultadoFinal = x;
        ejecutado = true;

        // calcular errorFinal: para cada iteración k y componente i:
        // errorFinal[k][i] = |(resultadoFinal[i] - tabla[k][i]) / resultadoFinal[i]|
        errorFinal.resize(tabla.size(), vector<double>(n,0.0));
        for(size_t k=0;k<tabla.size();k++){
            for(int i=0;i<n;i++){
                if(abs(resultadoFinal[i]) > 1e-15)
                    errorFinal[k][i] = fabs((resultadoFinal[i] - tabla[k][i]) / resultadoFinal[i]);
                else
                    errorFinal[k][i] = 0.0;
            }
        }

        // ---- Imprimir resultados centrados ----
        int anchoTerminal = 125;
        int nColumnas = n;
        int colAncho = 1;
        int totalAncho = nColumnas * colAncho + 1;
        int margenIzq = max(0, (anchoTerminal - totalAncho)/2);

        cout << "\n" << string(margenIzq, ' ') << "  RESULTADO FINAL \n\n\n";
        cout << string(margenIzq, ' ') << "Iteraciones necesarias: " << iter << "\n\n";

        cout << fixed << setprecision(S->fixDecimales);
        for(int i=0;i<n;i++){
            cout << string(margenIzq, ' ');
            cout << "x" << i+1 << " = " << x[i] << "\n";
        }
    }

    void mostrarTablaCompleta() {
        if(S->n == 0){
            cout << "\n                                                      Debe ingresar un sistema primero.\n";
            return;
        }

        // Si la tabla aun no tiene datos, calculamos automáticamente
        if(tabla.empty()){
            cout << "\n                                                      ADVERTENCIA: Jacobi no se ha ejecutado. Calculando tabla automáticamente...\n";
            resolver();
            // NOTA: resolver() ya llenó tabla, errores y errorFinal
            // dejamos 'ejecutado' tal como quedó (true)
        }

        int n = S->n;
        int anchoTerminal = 140;
        int colAncho = 18; // ancho por columna para que quepan decimales
        int totalAncho = (n * 3 + 1) * (colAncho/2); // estimación
        int margenIzq = 15;

        cout << "\n" << string(margenIzq, ' ') << "                                                 TABLA DE ITERACIONES \n";
        // Cabecera: Iteracion | x1 x2 ... | Err x1 Err x2 ... | ErrF x1 ErrF x2 ...
        cout << string(margenIzq, ' ') << setw(10) << "Iter";
        for(int i=0;i<n;i++) cout << setw(colAncho) << ("x" + to_string(i+1));
        for(int i=0;i<n;i++) cout << setw(colAncho) << ("Err x" + to_string(i+1));
        for(int i=0;i<n;i++) cout << setw(colAncho) << ("ErrF x" + to_string(i+1));
        cout << "\n";

        // Filas
        cout << fixed << setprecision(S->fixDecimales);
        for(size_t k=0;k<tabla.size();k++){
            cout << string(margenIzq, ' ') << setw(10) << k;
            // valores x
            for(int i=0;i<n;i++) cout << setw(colAncho) << tabla[k][i];
            // errores relativos entre iteraciones (errores): note errores[0] corresponds to iter 0 (zeros)
            for(int i=0;i<n;i++){
                // errores has same size as tabla (we pushed at start)
                double val = (k < errores.size()) ? errores[k][i] : 0.0;
                cout << setw(colAncho) << val;
            }
            // errores respecto a la ultima (errorFinal)
            for(int i=0;i<n;i++){
                double valF = (k < errorFinal.size()) ? errorFinal[k][i] : 0.0;
                cout << setw(colAncho) << valF;
            }
            cout << "\n";
        }
    }

    void mostrarTablaRango(int ini,int fin){
        if(S->n == 0){
            cout << "                                                      Debe ingresar un sistema primero.\n";
            return;
        }

        if(tabla.empty()){
            cout << "\n                                                      ADVERTENCIA: Jacobi no se ha ejecutado. Calculando tabla automaticamente...\n";
            resolver();
        }

        if(ini<0 || fin>= (int)tabla.size() || ini>fin){
            cout << "                                                      ADVERTENCIA: Rango invalido. Iteraciones disponibles: 0 a " << tabla.size()-1 << "\n";
            return;
        }

        int n = S->n;
        int anchoTerminal = 140;
        int colAncho = 18;
        int totalAncho = (n * 3 + 1) * (colAncho/2);
        int margenIzq = 15;

        cout << "\n" << string(margenIzq, ' ') << "                           TABLA DE ITERACIONES EN RANGO\n";
        cout << string(margenIzq, ' ') << setw(10) << "Iter";
        for(int i=0;i<n;i++) cout << setw(colAncho) << ("x" + to_string(i+1));
        for(int i=0;i<n;i++) cout << setw(colAncho) << ("Err x" + to_string(i+1));
        for(int i=0;i<n;i++) cout << setw(colAncho) << ("ErrF x" + to_string(i+1));
        cout << "\n";

        cout << fixed << setprecision(S->fixDecimales);
        for(int k=ini;k<=fin;k++){
            cout << string(margenIzq, ' ') << setw(10) << k;
            for(int i=0;i<n;i++) cout << setw(colAncho) << tabla[k][i];
            for(int i=0;i<n;i++){
                double val = (k < errores.size()) ? errores[k][i] : 0.0;
                cout << setw(colAncho) << val;
            }
            for(int i=0;i<n;i++){
                double valF = (k < errorFinal.size()) ? errorFinal[k][i] : 0.0;
                cout << setw(colAncho) << valF;
            }
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
            for(size_t k=0;k<tabla.size();k++){
                archivo << "Iteracion " << k << ": ";
                for(int i=0;i<S->n;i++)
                    archivo << "x" << i+1 << "=" << tabla[k][i] << "  ";
                // errores relativos entre iteraciones
                for(int i=0;i<S->n;i++)
                    archivo << "Err x" << i+1 << "=" << ((k < errores.size()) ? errores[k][i] : 0.0) << "  ";
                // errores respecto a la ultima
                for(int i=0;i<S->n;i++)
                    archivo << "ErrF x" << i+1 << "=" << ((k < errorFinal.size()) ? errorFinal[k][i] : 0.0) << "  ";
                archivo << "\n";
            }
        } else {
            archivo << "\n                                                      No es posible generar resultados ni tabla dado a que el método no puede converger.\n";
        }
        archivo.close();
        cout << "\n                                                      Archivo 'resultado.txt' creado correctamente.\n";
    }
};

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





