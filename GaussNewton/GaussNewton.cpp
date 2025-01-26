#define _USE_MATH_DEFINES
#include <iostream>
#include <Windows.h>
#include <iomanip>

using namespace std;

struct Matrix {
private:
    unsigned int width;//n - number of perem
    unsigned int height;//m - number of func

public:
    double** matr;

    Matrix() {
        width = height = 0;
        matr = nullptr;
    }

    Matrix(unsigned int width, unsigned int height) {
        this->width = width;
        this->height = height;
        matr = new double* [width];
        for (unsigned int i = 0; i < width; i++) {
            matr[i] = new double[height];
            for (unsigned int j = 0; j < height; j++) {
                matr[i][j] = 0;
            }
        }
    }

    Matrix(const Matrix& A) {
        width = A.width;
        height = A.height;
        matr = new double* [width];
        for (unsigned int i = 0; i < width; i++) {
            matr[i] = new double[height];
            for (unsigned int j = 0; j < height; j++) {
                matr[i][j] = A.matr[i][j];
            }
        }
    }

    ~Matrix() {
        for (unsigned int i = 0; i < width; i++) {
            delete matr[i];
        }
        delete matr;
    }

    double GetWidth() {
        return width;
    }

    double GetHeight() {
        return height;
    }

    void Print() {
        for (unsigned int i = 0; i < height; i++) {
            for (unsigned int j = 0; j < width; j++) {
                cout << right << fixed << setprecision(3) << setw(7) << matr[j][i] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    void PrintFull() {
        for (unsigned int i = 0; i < height; i++) {
            for (unsigned int j = 0; j < width; j++) {
                cout << setprecision(10) << matr[j][i] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    //friend ofstream& operator<< (ofstream& out, const Matrix& outMatr);
    //friend ifstream& operator>> (ifstream& in, const Matrix& inMatr);

    void ReadFromConsole() {
        for (unsigned int i = 0; i < height; i++) {
            for (unsigned int j = 0; j < width; j++) {
                cin >> matr[j][i];
            }
        }
    }

    Matrix& operator=(const Matrix& newMatr) {
        if (this == &newMatr) {
            return *this;
        }
        width = newMatr.width;
        height = newMatr.height;
        matr = new double* [width];
        for (unsigned int i = 0; i < width; i++) {
            matr[i] = new double[height];
            for (unsigned int j = 0; j < height; j++) {
                matr[i][j] = newMatr.matr[i][j];
            }
        }
        return *this;
    }

    friend Matrix operator+(const Matrix& left, const Matrix& right);
    friend Matrix operator-(const Matrix& left, const Matrix& right);
    friend Matrix operator*(const Matrix& left, const Matrix& right);

    Matrix& operator+=(const Matrix& right) {
        for (unsigned int i = 0; i < width; i++) {
            for (unsigned int j = 0; j < height; j++) {
                matr[i][j] += right.matr[i][j];
            }
        }
        return *this;
    }

    Matrix& operator-=(const Matrix& right) {
        for (unsigned int i = 0; i < width; i++) {
            for (unsigned int j = 0; j < height; j++) {
                matr[i][j] -= right.matr[i][j];
            }
        }
        return *this;
    }

    Matrix Transpose() {
        Matrix res(height, width);
        for (unsigned int i = 0; i < width; i++) {
            for (unsigned int j = 0; j < height; j++) {
                res.matr[j][i] = matr[i][j];
            }
        }
        return res;
    }

    Matrix FindInverse() {
        Matrix res(width, height);
        Matrix tmp(*this);
        for (unsigned int i = 0; i < height; i++) {
            res.matr[i][i] = 1;
        }
        double div;
        for (unsigned int i = 0; i < height-1; i++) {
            for (unsigned int k = i+1; k < height; k++) {
                div = tmp.matr[i][k] / tmp.matr[i][i];
                //cout << " div = " << div << endl;
                for (unsigned int j = i; j < width; j++) {
                    tmp.matr[j][k] -= tmp.matr[j][i] * div;
                }
                for (unsigned int j = 0; j < width; j++) {
                    res.matr[j][k] -= res.matr[j][i] * div;
                }
            }
        }
        for (unsigned int i = height-1; i > 0; i--) {
            //cout << " i = " << i << endl;
            for (unsigned int k = 0; k < i; k++) {
                div = tmp.matr[i][k] / tmp.matr[i][i];
                //cout << " div = " << div << endl;
                //cout << " k = " << k << endl;
                for (unsigned int j = i; j < width; j++) {
                    tmp.matr[j][k] -= tmp.matr[j][i] * div;
                }
                for (unsigned int j = 0; j < width; j++) {
                    //cout << " j " << j << endl;
                    res.matr[j][k] -= res.matr[j][i] * div;
                }
            }
        }
        for (unsigned int i = 0; i < height; i++) {
            for (unsigned int j = 0; j < width; j++) {
                res.matr[j][i] /= tmp.matr[i][i];
            }
        }
        return res;
    }
};

Matrix operator+(const Matrix& left, const Matrix& right) {
    Matrix res(left.width, left.height);
    for (unsigned int i = 0; i < res.width; i++) {
        for (unsigned int j = 0; j < res.height; j++) {
            res.matr[i][j] = left.matr[i][j] + right.matr[i][j];
        }
    }
    return res;
}

Matrix operator-(const Matrix& left, const Matrix& right) {
    Matrix res(left.width, left.height);
    for (unsigned int i = 0; i < res.width; i++) {
        for (unsigned int j = 0; j < res.height; j++) {
            res.matr[i][j] = left.matr[i][j] - right.matr[i][j];
        }
    }
    return res;
}

Matrix operator*(const Matrix& left, const Matrix& right) {
    Matrix res(right.width, left.height);
    for (unsigned int i = 0; i < res.width; i++) {
        for (unsigned int j = 0; j < res.height; j++) {
            for (unsigned int k = 0; k < right.height; k++) {
                res.matr[i][j] += left.matr[k][j]*right.matr[i][k];
            }
        }
    }
    return res;
}

struct TestTable {
public:
    unsigned int numberOfCases;
    double* x;
    double* y;
    unsigned int numberOfParam;
    double (*f)(double, Matrix&);//type of solution

    TestTable(unsigned int numOfCases, double* X, double* Y, unsigned int numOfParams, double (*f)(double, Matrix&)) {
        numberOfCases = numOfCases;
        x = X;
        y = Y;
        numberOfParam = numOfParams;
        this->f = f;
    }
};

Matrix CalculateJakobiMatrix(const TestTable& testTable, Matrix& curParams, double h) {
    Matrix res(testTable.numberOfParam, testTable.numberOfCases);
    for (unsigned int i = 0; i < res.GetHeight(); i++) {
        for (unsigned int j = 0; j < res.GetWidth(); j++) {
            curParams.matr[j][0] -= h;
            res.matr[j][i] = testTable.y[i] - testTable.f(testTable.x[i], curParams);//f(x,..,beta[i]-h,..)
            curParams.matr[j][0] += 2 * h;
            res.matr[j][i] -= testTable.y[i] - testTable.f(testTable.x[i], curParams);//f(x,..,beta[i]+h,..) for derivative calc 
            curParams.matr[j][0] -= h;
            res.matr[j][i] /= (2 * h);
        }
    }
    return res;
}

double TestF1(double x, Matrix& param) {
    return x * param.matr[0][0] / (param.matr[1][0] + x);
}

int main()
{
    int counter = 0;
    double* Test1X = new double[7];
    Test1X[0] = 0.038; Test1X[1] = 0.194; Test1X[2] = 0.425; Test1X[3] = 0.626; Test1X[4] = 1.253; Test1X[5] = 2.500; Test1X[6] = 3.740;
    double* Test1Y = new double[7]; 
    Test1Y[0] = 0.050; Test1Y[1] = 0.127; Test1Y[2] = 0.094; Test1Y[3] = 0.2122; Test1Y[4] = 0.2729; Test1Y[5] = 0.2665; Test1Y[6] = 0.3317;
    TestTable Test1(7, Test1X, Test1Y, 2, TestF1);

    double eps = 0.001;
    double h = 0.01;
    Matrix curParams(Test1.numberOfParam,1);
        curParams.matr[0][0] = 0.9;
        curParams.matr[1][0] = 0.2;
    Matrix prevParams;
    Matrix Jak;
    Matrix Delta;
    Matrix tmp;
    Matrix errors(1,Test1.numberOfCases);
    double MaxParamDifference;
    cout << " Start Params " << endl;
    curParams.Print();
    do{
       MaxParamDifference = DBL_MIN;
       prevParams = curParams;
       Jak = CalculateJakobiMatrix(Test1, curParams, h);
      // Jak.Print();
       tmp = Jak.Transpose();
       //tmp.Print();
       Delta = tmp * Jak;
       //Delta.Print();
       Delta = Delta.FindInverse();
       //Delta.Print();// Jak.Transpose();
       Delta = Delta * tmp;
       //Delta.Print();
       for (unsigned int i = 0; i < Test1.numberOfCases; i++) {
           errors.matr[0][i] = Test1.y[i] - Test1.f(Test1.x[i], curParams);
       }
       Delta = Delta * errors;
       //Delta.Print();
       curParams += Delta.Transpose();
       cout << " Iter i = " << counter << endl;
       cout << " Errors " << endl;
       errors.Print();
       cout << " Current params " << endl;
       curParams.Print();
       for (unsigned int i = 0; i < Test1.numberOfParam; i++) {
           if (MaxParamDifference < abs(curParams.matr[i][0] - prevParams.matr[i][0])) {
               MaxParamDifference = abs(curParams.matr[i][0] - prevParams.matr[i][0]);
           }
       }
       counter++;
    } while (MaxParamDifference > eps && counter < 1000);
    
    if (counter < 1000) {
        cout << " Result " << endl;
        curParams.PrintFull();
    }

    return 0;
}
