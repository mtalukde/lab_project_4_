#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// Function to compute x*cos(3x) using series expansion
double series_cos3x(double x, double epsilon, int& iterations) {
    double sum = 0.0;
    double term = x;  // Z_0 = x
    int n = 0;

    do {
        sum += term;
        n++;

        // Compute next term using recurrence relation:
        // Z_n = Z_{n-1} * (-9x^2 / ((2n)(2n-1)))
        double denominator = (2.0 * n) * (2.0 * n - 1.0);
        term = term * (-9.0 * x * x / denominator);

    } while (fabs(term) >= epsilon && n < 1000);

    iterations = n;
    return sum;
}

// Standard library function
double library_cos3x(double x) {
    return x * cos(3.0 * x);
}

// Function to compute discrepancy δ
double calculate_delta(double series_value, double library_value) {
    return sqrt(fabs(series_value * series_value - library_value * library_value));
}

int main() {
    double epsilon, x_start, x_end, delta_x, x_ideal;

    // Input parameters
    cout << "Enter precision ε: ";
    cin >> epsilon;
    cout << "Enter interval start x_start: ";
    cin >> x_start;
    cout << "Enter interval end x_end: ";
    cin >> x_end;
    cout << "Enter step Δx: ";
    cin >> delta_x;
    // Ask for ideal x value
    cout << "\nEnter ideal x value (x_ideal) for detailed analysis: ";
    cin >> x_ideal;

    // Table 1: computations for different x values
    cout << "\nTable 1 - Computation for different x values (ε = " << epsilon << ")" << endl;
    cout << "=============================================" << endl;
    cout << setw(10) << "x" << setw(15) << "f(x)" << setw(15) << "F(x)"
         << setw(12) << "δ" << setw(12) << "Iterations" << endl;
    cout << "---------------------------------------------" << endl;

    for (double x = x_start; x <= x_end; x += delta_x) {
        int iterations;
        double series_val = series_cos3x(x, epsilon, iterations);
        double library_val = library_cos3x(x);
        double delta = calculate_delta(series_val, library_val);

        cout << fixed << setprecision(6);
        cout << setw(10) << x << setw(15) << series_val << setw(15) << library_val
             << setw(12) << delta << setw(12) << iterations << endl;
    }



    // Table 2: computations for different precision values with ideal x
    cout << "\nTable 2 - Computation for different precision values (x_ideal = " << x_ideal << ")" << endl;
    cout << "==================================================================" << endl;
    cout << setw(12) << "ε" << setw(15) << "f(x_ideal)" << setw(15) << "F(x_ideal)"
         << setw(12) << "δ" << setw(12) << "Iterations" << endl;
    cout << "------------------------------------------------------------------" << endl;

    double epsilons[] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};

    for (double eps : epsilons) {
        int iterations;
        double series_val = series_cos3x(x_ideal, eps, iterations);
        double library_val = library_cos3x(x_ideal);
        double delta = calculate_delta(series_val, library_val);

        cout << scientific << setprecision(0);
        cout << setw(12) << eps;
        cout << fixed << setprecision(8);
        cout << setw(15) << series_val << setw(15) << library_val
             << setw(12) << delta << setw(12) << iterations << endl;
    }

    return 0;
}
