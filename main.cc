#include <iostream>

struct Matrix{
    size_t sz;
    double* arr; 

    Matrix(size_t sz) {
        this->sz = sz;
        this->arr = new double[sz * sz];
    }

    inline double& index(size_t x, size_t y) {
        return this->arr[sz*x + y];
    }

    ~Matrix() {
        delete arr;
    }

    void print() {
        printf("[\n");
        for(int i = 0; i < sz; i++) {
            printf("[");
            for(int j = 0; j < sz; j++) {
                printf("%.2f,", this->index(j,i));
            }
            printf("]\n");
        }
        printf("]\n");
    }
};

// Multiply a and b, put result in res
void mmult(Matrix* a, Matrix* b, Matrix* res) {
    if (a->sz != b->sz || a->sz != res->sz) {
        printf("Wrong sizes\n");
        exit(1);
    }

    size_t n = a->sz;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                res->index(i, j) += a->index(i,k) * b->index(k, j);
            }
        }
    }


}

int main() {
    Matrix* a = new Matrix(3);
    a->index(0,0) = 1;
    a->index(0,2) = 1;
    a->index(1,2) = 1;
    a->print();

    Matrix* b = new Matrix(3);
    b->index(0,0) = 2;
    b->index(1,1) = 2;
    b->index(2,2) = 2;


    Matrix* res = new Matrix(3);

    mmult(a, b, res);

    res->print();
}