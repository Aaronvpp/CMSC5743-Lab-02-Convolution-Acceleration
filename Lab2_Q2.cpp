#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cassert>
#include <sys/time.h>
using namespace std;

float G[12] = {1, 0, 0,
               0.5, 0.5, 0.5,
               0.5, -0.5, 0.5,
               0, 0, 1};
float GT[12] = {1, 0.5, 0.5, 0,
                0, 0.5, -0.5, 0,
                0, 0.5, 0.5, 1};

double get_time() {
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}
void dot(float* A, int row_A, int col_A, float* B, int row_B, int col_B, float* C) {
    assert(col_A == row_B);
    for (int i = 0; i < row_A; i++)
        for (int j = 0; j < col_B; j++)
            //  C[i*col_A + j] = 0;
            for (int p = 0; p < col_A; p++)
                C[i * col_B + j] += A[i * col_A + p] * B[p * col_B + j];
}

void multi(float* A, int row_A, int col_A, float* B, int row_B, int col_B, float* C) {
    assert(row_A == row_B && col_A == col_B);
    for (int i = 0; i < row_A; i++)
        for (int j = 0; j < col_A; j++)
            C[col_A * i + j] = A[col_A * i + j] * B[col_A * i + j];
}

void winograd_transforme_g(float* g, float* transformed_g) {
    float Gg[12] = {0};

    dot(G, 4, 3, g, 3, 3, Gg); //U=GgGT
    dot(Gg, 4, 3, GT, 3, 4, transformed_g);
}

void winograd(float* U, float* d, float* result) {
    float BTd[16] = {0};
    float V[16] = {0};
    float UV[16] = {0};
    float ATUV[8] = {0};

    // start dot(BT, 4, 4, d, 4, 4, BTd);
    BTd[0] = d[0] - d[8];
    BTd[1] = d[1] - d[9];
    BTd[2] = d[2] - d[10];
    BTd[3] = d[3] - d[11];

    BTd[4] = d[4] + d[8];
    BTd[5] = d[5] + d[9];
    BTd[6] = d[6] + d[10];
    BTd[7] = d[7] + d[11];

    BTd[8] = -d[4] + d[8];
    BTd[9] = -d[5] + d[9];
    BTd[10] = -d[6] + d[10];
    BTd[11] = -d[7] + d[11];

    BTd[12] = d[4] - d[12];
    BTd[13] = d[5] - d[13];
    BTd[14] = d[6] - d[14];
    BTd[15] = d[7] - d[15];
    // end dot(BT, 4, 4, d, 4, 4, BTd);

    // start dot(BTd, 4, 4, B, 4, 4, V);
    V[0] = BTd[0] - BTd[2];
    V[4] = BTd[4] - BTd[6];
    V[8] = BTd[8] - BTd[10];
    V[12] = BTd[12] - BTd[14];

    V[1] = BTd[1] + BTd[2];
    V[5] = BTd[5] + BTd[6];
    V[9] = BTd[9] + BTd[10];
    V[13] = BTd[13] + BTd[14];

    V[2] = -BTd[1] + BTd[2];
    V[6] = -BTd[5] + BTd[6];
    V[10] = -BTd[9] + BTd[10];
    V[14] = -BTd[13] + BTd[14];

    V[3] = BTd[1] - BTd[3];
    V[7] = BTd[5] - BTd[7];
    V[11] = BTd[9] - BTd[11];
    V[15] = BTd[13] - BTd[15];
    // end dot(BTd, 4, 4, B, 4, 4, V);

    multi(U, 4, 4, V, 4, 4, UV);

    // start dot(AT, 2, 4, UV, 4, 4, ATUV);
    ATUV[0] = UV[0] + UV[4] + UV[8];
    ATUV[1] = UV[1] + UV[5] + UV[9];
    ATUV[2] = UV[2] + UV[6] + UV[10];
    ATUV[3] = UV[3] + UV[7] + UV[11];

    ATUV[4] = UV[4] - UV[8] - UV[12];
    ATUV[5] = UV[5] - UV[9] - UV[13];
    ATUV[6] = UV[6] - UV[10] - UV[14];
    ATUV[7] = UV[7] - UV[11] - UV[15];
    // end dot(AT, 2, 4, UV, 4, 4, ATUV);

    // start dot(ATUV, 2, 4, A, 4, 2, result);
    result[0] += (ATUV[0] + ATUV[1] + ATUV[2]);
    result[2] += (ATUV[4] + ATUV[5] + ATUV[6]);
    result[1] += (ATUV[1] - ATUV[2] - ATUV[3]);
    result[3] += (ATUV[5] - ATUV[6] - ATUV[7]);
    // end dot(ATUV, 2, 4, A, 4, 2, result);
}
//float get_data(float* data_im, int c, int im_w, int im_h, int row, int col, int ph, int pw) {
//
//    row = row - ph;
//    col = col - pw;
//    if(row < 0 || row >= im_h || col < 0 || col >= im_w) {
//        return 0;
//    }
//
//    return data_im[c * im_w * im_h + row * im_w + col];
//
//}
//void im2col( float* data_im,
//            const int im_c,
//            const int im_w,
//            const int im_h,
//            const int kw,
//            const int kh,
//            const int ph,
//            const int pw,
//            const int sh,
//            const int sw,
//            float* data_col,
//            const int col_w,
//            const int col_h) {
//
//    // win_w and win_h are the stop times of the kernel in the image.
//    int win_w = (im_w + 2 * pw - kw + 1) / sw;
//    int win_h = (im_h + 2 * ph - kh + 1) / sh;
//    int x;
//    int y;
//
//    for (int i = 0; i< col_h; i++) {
//
//        x = i % win_w;
//        y = i / win_h;
//        for(int j = 0; j < col_w; j++) {
//            int c = j / (kw * kh);
//            int kj = j % kw;
//            int ki = j / kw;
//
//            int row = y * sh + ki;
//            int col = x * sw + kj;
//
//            data_col[i * col_w + j] = get_data(data_im, c, im_w, im_h, row, col, ph, pw);
//        }
//    }
//    cout << data_col[0] <<endl;
//
//}

int main(int argc, char const* argv[]) {
    float g[] = {1, 2, 3,
                 4, 5, 6,
                 7, 8, 9};//Kernel

    int SIZE = 56;
    //Build the feature map;
    float MatrixA[SIZE*SIZE];
    for (int i = 0; i < (SIZE*SIZE); i++){
            MatrixA[i] = rand() % 10 + 1;
    }

    int N=4, M = 3*3*N;
    float *test_col = new float [M];

//    cout << size(d) << endl;
//    im2col(d,0,4,4,3,3,0,0,1,1,test_col,9,4);
//    cout << size(d) << endl;
//    for (int i = 0; i < M; i++){
//        cout << test_col[i] << endl;
//    }
    //Splitting the matrix:
    float Split_Matrix[26*26][16];
    float avg_time_winograd = 0.0f;
    for (int i = 0; i < (26*26); i++) {
        Split_Matrix[i][0] = MatrixA[i];Split_Matrix[i][1] = MatrixA[i+1];Split_Matrix[i][2] = MatrixA[i+2];Split_Matrix[i][3] = MatrixA[i+3];
        Split_Matrix[i][4] = MatrixA[i+56];Split_Matrix[i][5] = MatrixA[i+57];Split_Matrix[i][6] = MatrixA[i+58];Split_Matrix[i][7] = MatrixA[i+59];
        Split_Matrix[i][8] = MatrixA[i+112];Split_Matrix[i][9] = MatrixA[i+113];Split_Matrix[i][10] = MatrixA[i+114];Split_Matrix[i][11] = MatrixA[i+115];
        Split_Matrix[i][12] = MatrixA[i+168];Split_Matrix[i][13] = MatrixA[i+169];Split_Matrix[i][14] = MatrixA[i+170];Split_Matrix[i][15] = MatrixA[i+171];
    }

    //Assign each splited matrix TO d[]:4*4 for Winograd
    cout << "This is the result for 3 channel input feature map: 56*56, kernal: 3*3 and 64 output channel Winograd" << endl;
    auto t = get_time();
    t = get_time();
    for (int k = 0; k < 64; k++) {
        for (int c = 0; c < 3; c++) {//3 input channels
            for (int i = 0; i < (26 * 26); i++) {
                float d[16] = {0};
                d[0] = Split_Matrix[i][0];
                d[1] = Split_Matrix[i][1];
                d[2] = Split_Matrix[i][2];
                d[3] = Split_Matrix[i][3];
                d[4] = Split_Matrix[i][4];
                d[5] = Split_Matrix[i][5];
                d[6] = Split_Matrix[i][6];
                d[7] = Split_Matrix[i][7];
                d[8] = Split_Matrix[i][8];
                d[9] = Split_Matrix[i][9];
                d[10] = Split_Matrix[i][10];
                d[11] = Split_Matrix[i][11];
                d[12] = Split_Matrix[i][12];
                d[13] = Split_Matrix[i][13];
                d[14] = Split_Matrix[i][14];
                d[15] = Split_Matrix[i][15];
                float result[4] = {0};
                float *transformed_g = (float *) calloc(16, sizeof(float));
                winograd_transforme_g(g, transformed_g);
                winograd(transformed_g, d, result);
                //adding 3 input channels' result together, not the final size, just for time calculation
                float channel_3_result[3][4];
                if (c >= 1){
                    channel_3_result[c][0] = channel_3_result[c-1][0] + result[0];
                    channel_3_result[c][1] = channel_3_result[c-1][1] + result[1];
                    channel_3_result[c][2] = channel_3_result[c-1][2] + result[2];
                    channel_3_result[c][3] = channel_3_result[c-1][3] + result[3];
                }
                else {
                    channel_3_result[c][0] = result[0];
                    channel_3_result[c][1] = result[1];
                    channel_3_result[c][2] = result[2];
                    channel_3_result[c][3] = result[3];
                }

                //printing all the final output matrix
//                for (int j = 0; j < 4; j++) {
//                    cout << channel_3_result[c][j] << endl;
//                }
            }
        }
    }
    avg_time_winograd+=get_time() - t;
    cout << "This is the result for input channel: 3, feature map: 56*56, kernal: 3*3 and  output channel: 64 Winograd\n" << endl;
    cout << "(If you want to see the output matrix, please go to the Lab2_Q2 to uncomment 'printing all the final output matrix' around line 245)\n" << endl;
    printf("The running time for Winograd is %f",avg_time_winograd);
    //end of Winograd


    return 0;
}
