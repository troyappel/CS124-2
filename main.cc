#include <iostream>
#include <assert.h>

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
            printf("   [");
            for(int j = 0; j < sz; j++) {
                printf("%.2f,", this->index(j,i));
            }
            printf("]\n");
        }
        printf("]\n");
    }
};

struct MatPak {
    Matrix* m;
    size_t s_x;
    size_t s_y;

    size_t e_x;
    size_t e_y;
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

void madd_s(MatPak a, MatPak b, MatPak res) {
    assert(a.m->sz >= a.e_x);
    assert(a.m->sz >= a.e_y);

    assert(b.m->sz >= b.e_x);
    assert(b.m->sz >= b.e_y);

    assert(a.e_y - a.s_y == b.e_y - b.s_y);
    assert(a.e_x - a.s_x == b.e_x - b.s_x);

    for(int i = 0; i < a.e_x - a.s_x; i++) {
        for (int j = 0; j < a.e_y - a.s_y; j++) {
            res.m->index(i + res.s_x, j + res.s_y) 
                = a.m->index(i + a.s_x, j + a.s_y)
                + b.m->index(i + b.s_x, j + b.s_y);
        }
    }

}

// This is OK in place
void madd(Matrix* a, Matrix* b, Matrix* res) {
    madd_s(MatPak{a, 0, 0, a->sz, a->sz}, MatPak{b, 0, 0, b->sz, b->sz}, MatPak{res, 0, 0, res->sz, res->sz});
}



void mmult_strassen() {

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
    b->index(1,2) = 3.14;
    b->print();


    Matrix* res = new Matrix(3);

    madd_s(MatPak{a, 1, 1, 3, 3}, MatPak{b, 1, 0, 3, 2}, MatPak{res, 0, 0, 2, 2});

    res->print();
}