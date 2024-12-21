#include <iostream>
#include <time.h>

using namespace std;
void showMatrix(double** matrixA, double* matrixB, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << matrixA[i][j] << ' ';
        }
        cout << "| " << matrixB[i] << endl;
    }
    cout << endl;
}

void showMatrix(double* matrix, int n)
{
    for (int i = 0; i < n; i++)
    {
        cout << matrix[i] << endl;
    }
    cout << endl;
}

double* Gauss(double** A, double* B, int n)
{
    double* X = new double[n];

    for (int k = 0; k < n; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            if (abs(A[k][k]) < abs(A[i][k]))
            {
                swap(A[i], A[k]);
                swap(B[i], B[k]);
            }
        }

        double el = A[k][k];

        if (el == 0)
        {
            cerr << "the leading element can't be equal to 0";
            exit(0);
        }

        for (int i = k; i < n; i++)
        {
            A[k][i] /= el;
        }
        B[k] /= el;

        for (int i = k + 1; i < n; i++)
        {
            double num = A[i][k];
            for (int j = k; j < n; j++)
            {
                A[i][j] -= num * A[k][j];
            }
            B[i] -= num * B[k];
        }
    }

    for (int k = n - 1; k >= 0; k--)
    {
        X[k] = B[k];
        for (int i = n - 1; i > k; i--)
        {
            X[k] -= A[k][i] * X[i];
        }
    }

    return X;
}

double* residualVector(double** A, double* B, double* X, int n)
{
    double* F = new double[n];

    for (int i = 0; i < n; i++)
    {
        double sum = 0;

        for (int j = 0; j < n; j++)
        {
            sum += A[i][j] * X[j];
        }

        F[i] = sum - B[i];
    }

    return F;
}

double norm(double* vec, int n)
{
    double result = vec[0];

    for (int i = 0; i < n; i++)
    {
        result = max(vec[i], result);
    }

    return result;
}

double* auxiliaryB(double** A, double* X, int n)
{
    double* B2 = new double[n] {};

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            B2[i] += A[i][j] * X[j];
        }
    }

    return B2;
}

double accuracy(double** A, double* X1, int n)
{
    double acc = 0;

    double* diffX = new double[n] {};
    double* B2 = auxiliaryB(A, X1, n);
    double* X2 = Gauss(A, B2, n);

    for (int i = 0; i < n; i++)
    {
        diffX[i] = X2[i] - X1[i];
    }

    acc = norm(diffX, n) / norm(X1, n);

    return acc;
}

int main()
{
    int n = 3;

    double** A = new double* [n];

    for (int i = 0; i < n; i++)
    {
        A[i] = new double [n] {};
    }

    double* B = new double [n] {};

    A[0][0] = 2.30;
    A[0][1] = 5.70;
    A[0][2] = -0.80;
    A[1][0] = 3.50;
    A[1][1] = -2.70;
    A[1][2] = 5.30;
    A[2][0] = 1.70;
    A[2][1] = 2.30;
    A[2][2] = -1.80;

    B[0] = -6.49;
    B[1] = 19.20;
    B[2] = -5.09;

    double** copyOfA = new double* [n];

    for (int i = 0; i < n; i++)
    {
        copyOfA[i] = new double [n] {};

        for (int j = 0; j < n; j++)
        {
            copyOfA[i][j] = A[i][j];
        }
    }

    double* copyOfB = new double[n];
    for (int i = 0; i < n; i++)
    {
        copyOfB[i] = B[i];
    }

    cout << "\tA | B:" << endl;
    showMatrix(A, B, n);

    double* X = Gauss(A, B, n);

    cout << "\tX:" << endl;
    showMatrix(X, n);

    double* F = residualVector(copyOfA, copyOfB, X, n);
    double normOfF = norm(F, n);

    cout << "norm of F = " << normOfF << endl;

    cout << "error = " << accuracy(copyOfA, X, n) << endl;

    for (int i = 0; i < n; i++)
    {
        delete[] A[i];
    }

    delete[] A;
    delete[] B;
    delete[] X;
    delete[] F;
}