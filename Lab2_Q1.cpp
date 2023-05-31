#include <sys/time.h>
#include <iostream>
using namespace std;


double get_time() {
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}

void matmul(int SIZE, int **MatrixA, int **MatrixB, int **MatrixC){
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE; j++){
            MatrixC[i][j] = 0;
            for (int k = 0; k < SIZE; k++){
                MatrixC[i][j] = MatrixC[i][j]+MatrixA[i][k] * MatrixB[k][j];
            }
        }
    }
}
void PrintMatrix(int **MatrixA, int N){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            cout << MatrixA[i][j] << " ";
        }
        cout << endl;
    }
}

void Matrix_Sum(int N, int** MatrixA, int** MatrixB, int** Sum_Matrix){

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Sum_Matrix[i][j] = MatrixA[i][j] + MatrixB[i][j];
        }
    }
}

void Matrix_Sub(int N, int** MatrixA, int** MatrixB, int** Sub_Matrix){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Sub_Matrix[i][j] = MatrixA[i][j] -MatrixB[i][j];
        }
    }
}

void Matrix_Mul(int N, int** MatrixA, int** MatrixB, int** Mul_Matrix){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Mul_Matrix[i][j] = 0;
            for (int k = 0; k < N; k++){
                Mul_Matrix[i][j] = Mul_Matrix[i][j] + MatrixA[i][k] * MatrixB[k][j];
            }
        }
    }
}

void Strassen(int N, int** MatrixA, int** MatrixB, int** MatrixC){

    int** MatrixA11;
    int** MatrixA12;
    int** MatrixA21;
    int** MatrixA22;

    int** MatrixB11;
    int** MatrixB12;
    int** MatrixB21;
    int** MatrixB22;

    int** MatrixC11;
    int** MatrixC12;
    int** MatrixC21;
    int** MatrixC22;

    MatrixA11 = new int*[N/2];
    MatrixA12 = new int*[N/2];
    MatrixA21 = new int*[N/2];
    MatrixA22 = new int*[N/2];

    MatrixB11 = new int*[N/2];
    MatrixB12 = new int*[N/2];
    MatrixB21 = new int*[N/2];
    MatrixB22 = new int*[N/2];

    MatrixC11 = new int*[N/2];
    MatrixC12 = new int*[N/2];
    MatrixC21 = new int*[N/2];
    MatrixC22 = new int*[N/2];
    for (int i = 0; i < N/2; i++)
    {
        MatrixA11[i] = new int[N/2];
        MatrixA12[i] = new int[N / 2];
        MatrixA21[i] = new int[N / 2];
        MatrixA22[i] = new int[N / 2];

        MatrixB11[i] = new int[N / 2];
        MatrixB12[i] = new int[N / 2];
        MatrixB21[i] = new int[N / 2];
        MatrixB22[i] = new int[N / 2];

        MatrixC11[i] = new int[N / 2];
        MatrixC12[i] = new int[N / 2];
        MatrixC21[i] = new int[N / 2];
        MatrixC22[i] = new int[N / 2];
    }

    for (int i = 0; i < N / 2; i++){
        for (int j = 0; j < N / 2; j++){
            MatrixA11[i][j] = MatrixA[i][j];
            MatrixA12[i][j] = MatrixA[i][j + N / 2];
            MatrixA21[i][j] = MatrixA[i + N / 2][j];
            MatrixA22[i][j] = MatrixA[i + N / 2][j + N / 2];

            MatrixB11[i][j] = MatrixB[i][j];
            MatrixB12[i][j] = MatrixB[i][j + N / 2];
            MatrixB21[i][j] = MatrixB[i + N / 2][j];
            MatrixB22[i][j] = MatrixB[i + N / 2][j + N / 2];
        }
    }

    int** MatrixS1=new int*[N/2];
    int** MatrixS2 = new int*[N/2];
    int** MatrixS3 = new int*[N/2];
    int** MatrixS4 = new int*[N / 2];
    int** MatrixS5 = new int*[N / 2];
    int** MatrixS6 = new int*[N / 2];
    int** MatrixS7 = new int*[N / 2];
    int** MatrixS8 = new int*[N / 2];
    int** MatrixS9 = new int*[N / 2];
    int** MatrixS10 = new int*[N / 2];

    for (int i = 0; i < N / 2; i++)
    {
        MatrixS1[i] = new int[N / 2];
        MatrixS2[i] = new int[N / 2];
        MatrixS3[i] = new int[N / 2];
        MatrixS4[i] = new int[N / 2];
        MatrixS5[i] = new int[N / 2];
        MatrixS6[i] = new int[N / 2];
        MatrixS7[i] = new int[N / 2];
        MatrixS8[i] = new int[N / 2];
        MatrixS9[i] = new int[N / 2];
        MatrixS10[i] = new int[N / 2];
    }

    Matrix_Sub(N/2, MatrixB12, MatrixB22, MatrixS1);//S1 = B12 - B22
    Matrix_Sum(N / 2, MatrixA11, MatrixA12, MatrixS2);//S2 = A11 + A12
    Matrix_Sum(N / 2, MatrixA21, MatrixA22, MatrixS3);//S3 = A21 + A22
    Matrix_Sub(N / 2, MatrixB21, MatrixB11, MatrixS4);//S4 = B21 - B11
    Matrix_Sum(N / 2, MatrixA11, MatrixA22, MatrixS5);//S5 = A11 + A22
    Matrix_Sum(N / 2, MatrixB11, MatrixB22, MatrixS6);//S6 = B11 + B22
    Matrix_Sub(N / 2, MatrixA12, MatrixA22, MatrixS7);//S7 = A12 - A22
    Matrix_Sum(N / 2, MatrixB21, MatrixB22, MatrixS8);//S8 = B21 + B22
    Matrix_Sub(N / 2, MatrixA11, MatrixA21, MatrixS9);//S9 = A11 - A21
    Matrix_Sum(N / 2, MatrixB11, MatrixB12, MatrixS10);//S10 = B11 + B12

    int** MatrixP1 = new int*[N / 2];
    int** MatrixP2 = new int*[N / 2];
    int** MatrixP3 = new int*[N / 2];
    int** MatrixP4 = new int*[N / 2];
    int** MatrixP5 = new int*[N / 2];
    int** MatrixP6 = new int*[N / 2];
    int** MatrixP7 = new int*[N / 2];

    for (int i = 0; i < N / 2; i++)
    {
        MatrixP1[i] = new int[N / 2];
        MatrixP2[i] = new int[N / 2];
        MatrixP3[i] = new int[N / 2];
        MatrixP4[i] = new int[N / 2];
        MatrixP5[i] = new int[N / 2];
        MatrixP6[i] = new int[N / 2];
        MatrixP7[i] = new int[N / 2];
    }
    Matrix_Mul(N / 2, MatrixA11, MatrixS1, MatrixP1);//P1 = A11 • S1
    Matrix_Mul(N / 2, MatrixS2, MatrixB22, MatrixP2);//P2 = S2 • B22
    Matrix_Mul(N / 2, MatrixS3, MatrixB11, MatrixP3);//P3 = S3 • B11
    Matrix_Mul(N / 2, MatrixA22, MatrixS4, MatrixP4);//P4 = A22 • S4
    Matrix_Mul(N / 2, MatrixS5, MatrixS6, MatrixP5);//P5 = S5 • S6
    Matrix_Mul(N / 2, MatrixS7, MatrixS8, MatrixP6);//P6 = S7 • S8
    Matrix_Mul(N / 2, MatrixS9, MatrixS10, MatrixP7);//P7 = S9 • S10

    //根据以上7个结果计算C矩阵
    Matrix_Sum(N / 2, MatrixP5, MatrixP4, MatrixC11); //C11 = P5 + P4 - P2 + P6
    Matrix_Sub(N / 2, MatrixC11, MatrixP2, MatrixC11);
    Matrix_Sum(N / 2, MatrixC11, MatrixP6, MatrixC11);
    Matrix_Sum(N / 2, MatrixP1, MatrixP2, MatrixC12);//C12 = P1 + P2
    Matrix_Sum(N / 2, MatrixP3, MatrixP4, MatrixC21);	//C21 = P3 + P4
    Matrix_Sum(N / 2, MatrixP5, MatrixP1, MatrixC22);	//C22 = P5 + P1 - P3 - P7
    Matrix_Sub(N / 2, MatrixC22, MatrixP3, MatrixC22);
    Matrix_Sub(N / 2, MatrixC22, MatrixP7, MatrixC22);

    for (int i = 0; i < N / 2; i++){
        for (int j = 0; j < N / 2; j++){
            MatrixC[i][j] = MatrixC11[i][j];
            MatrixC[i][j+N/2] = MatrixC12[i][j];
            MatrixC[i+N/2][j] = MatrixC21[i][j];
            MatrixC[i+N/2][j+N/2] = MatrixC22[i][j];
        }
    }
}

void NormalMul_Matrix(int N, int **MatrixA, int **MatrixB, int **MatrixC){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            MatrixC[i][j] = 0;
            for (int k = 0; k < N; k++){
                MatrixC[i][j] = MatrixC[i][j]+MatrixA[i][k] * MatrixB[k][j];
            }
        }
    }
}
void Init_Matrix(int N,int** MatrixA, int** MatrixB){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            MatrixA[i][j] = rand() % 10 + 1;
            MatrixB[i][j] = rand() % 10 + 1;
        }
    }

}

void Test_Matrix(int N, int** MatrixA, int** MatrixB){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            MatrixA[i][j] = 1;
            MatrixB[i][j] = 2;
        }
    }
}

int main(){
    int N;
    cout << "Please enter the matrix size: ";
    cin >> N;
    float avg_time_strassen = 0.0f;
    float avg_time_matmul = 0.0f;
    for (int m = 0; m < 32; m++) {
        int **MatrixA = new int *[N];
        int **MatrixB = new int *[N];
        int **MatrixC = new int *[N];
        int** MatrixCm = new int*[N];

        for (int i = 0; i < N; i++) {
            MatrixA[i] = new int[N];
            MatrixB[i] = new int[N];
            MatrixC[i] = new int[N];
            MatrixCm[i] = new int[N];

        }

        Init_Matrix(N, MatrixA, MatrixB);
        //Test_Matrix(N, MatrixA, MatrixB);
        cout << "A Matrix: (For saving your time, I commented the //PrintMatrix(MatrixA, N); in the following part)" << endl;
        //PrintMatrix(MatrixA, N);

        cout << "B Matrix: (If you want to see the details of the matrix, you can go to the Q1.cpp to uncomment the PrintMatrix() in int main().) " << endl;
        //PrintMatrix(MatrixB, N);

        auto t = get_time();
        cout << "Strassen:" << endl;
        t = get_time();
        Strassen(N, MatrixA, MatrixB, MatrixC);
        avg_time_strassen+=get_time() - t;
        //PrintMatrix(MatrixC, N);

        cout << "matmul:" << endl;
        matmul(N, MatrixA, MatrixB, MatrixCm);
        avg_time_matmul+=get_time() - t;
        //PrintMatrix(MatrixCm, N);


    }
    printf("Input matrix size is : %d, Avg Time for matmul: %f\n",N,avg_time_matmul / 32);

    printf("Input matrix size is : %d, Avg Time for Strassen Alg.: %f\n",N,avg_time_strassen / 32);
}
