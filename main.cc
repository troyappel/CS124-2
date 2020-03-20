#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <string.h> 

size_t CUTOFF = 10;

struct Matrix{
    size_t sz;
    double* arr; 

    Matrix(size_t sz) {
        this->sz = sz;
        this->arr = new double[sz * sz];
        this->clear();
    }

    inline double& index(size_t x, size_t y) {
        return this->arr[sz*x + y];
    }

    inline void clear() {
        memset(arr, 0, sz*sz*sizeof(double));
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

struct MatPak;

MatPak make(Matrix* m,size_t s_x,size_t s_y, size_t e_x, size_t e_y);



struct MatPak {
    Matrix* m;
    size_t s_x;
    size_t s_y;

    size_t e_x;
    size_t e_y;

    double nothing;

    bool padded;

    inline double& index(size_t x, size_t y) {
        if(padded && (x == e_x - s_x - 1 || y == e_y - s_y - 1)) {
            nothing = 0;
            return nothing;
        }
        return m->index(s_x + x, s_y + y);
    };

    inline MatPak subPak(size_t s_x,size_t s_y, size_t e_x, size_t e_y) {
        return make(m, this->s_x + s_x, this->s_y + s_y, this->s_x + e_x, this->s_y + e_y);
    };

};

// Always returns an even-size MatPak, padded with 0 if bounds given are odd 
MatPak make(Matrix* m,size_t s_x,size_t s_y, size_t e_x, size_t e_y) {
    bool padded;
    size_t n = e_x - e_y;
    if (n & 1) {
        e_x++;
        e_y++;
        padded = true;
    } else {
        padded = false;
    }

    return MatPak{m, s_x, s_y, e_x, e_y, padded};

}


void mmult_s(MatPak a, MatPak b, MatPak res) {
    size_t n = a.e_x - a.s_x;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                res.index(i, j) += a.index(i,k) * b.index(k, j);
            }
        }
    }
}


// Multiply a and b, put result in res
void mmult(Matrix* a, Matrix* b, Matrix* res) {
    res->clear();
    mmult_s(MatPak{a, 0, 0, a->sz, a->sz}, MatPak{b, 0, 0, b->sz, b->sz}, MatPak{res, 0, 0, res->sz, res->sz});
}



void madd_s(MatPak a, MatPak b, MatPak res) {
    // assert(a.m->sz >= a.e_x);
    // assert(a.m->sz >= a.e_y);

    // assert(b.m->sz >= b.e_x);
    // assert(b.m->sz >= b.e_y);

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

void msub_s(MatPak a, MatPak b, MatPak res) {
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
                - b.m->index(i + b.s_x, j + b.s_y);
        }
    }

}

// This is OK in place
void madd(Matrix* a, Matrix* b, Matrix* res) {
    madd_s(MatPak{a, 0, 0, a->sz, a->sz}, MatPak{b, 0, 0, b->sz, b->sz}, MatPak{res, 0, 0, res->sz, res->sz});
}

void cmult(Matrix* a, double c) {
    size_t n = a->sz;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a->index(i,j) *= c;
        }
    }
}

void mmult_strassen(Matrix* a, Matrix* b, Matrix* res);
void mmult_strassen_s(MatPak a, MatPak b, MatPak res);

void mmult_strassen(Matrix* a, Matrix* b, Matrix* res) {
    res->clear();
    mmult_strassen_s(MatPak{a, 0, 0, a->sz, a->sz}, MatPak{b, 0, 0, b->sz, b->sz}, MatPak{res, 0, 0, res->sz, res->sz});
}

void mmult_strassen_s(MatPak a, MatPak b, MatPak res) {
    size_t n2 = (a.e_x - a.s_x) / 2;
    size_t n = a.e_x - a.s_x;

    if (n <= 1) {
        res.index(0,0) = a.index(0,0) * b.index(0,0);
        return;
    }

    if (n <= CUTOFF) {
        mmult_s(a,b,res);
        return;
    }

    Matrix* M1 = new Matrix(n2);
    Matrix* M2 = new Matrix(n2);
    Matrix* M3 = new Matrix(n2);
    Matrix* M4 = new Matrix(n2);
    Matrix* M5 = new Matrix(n2);
    Matrix* M6 = new Matrix(n2);
    Matrix* M7 = new Matrix(n2);

    Matrix* scratch1 = new Matrix(n2);
    Matrix* scratch2 = new Matrix(n2);

    madd_s(a.subPak( 0, 0, n2, n2), a.subPak( n2, n2, n, n), MatPak{scratch1, 0, 0, n2, n2});
    madd_s(b.subPak( 0, 0, n2, n2), b.subPak( n2, n2, n, n), MatPak{scratch2, 0, 0, n2, n2});
    mmult_strassen(scratch1, scratch2, M1);


    madd_s(a.subPak( 0, n2, n2, n), a.subPak( n2, n2, n, n), MatPak{scratch1, 0, 0, n2, n2});
    mmult_strassen_s(MatPak{scratch1, 0,0,n2,n2}, b.subPak(0,0,n2,n2), MatPak{M2, 0, 0, M2->sz, M2->sz});

    msub_s(b.subPak( n2, 0, n, n2), b.subPak( n2, n2, n, n), MatPak{scratch1, 0, 0, n2, n2});
    mmult_strassen_s(MatPak{scratch1, 0,0,n2,n2}, a.subPak(0,0,n2,n2), MatPak{M3, 0, 0, n2, n2});

    msub_s(b.subPak( 0, n2, n2, n), b.subPak( 0, 0, n2, n2), MatPak{scratch1, 0, 0, n2, n2});
    mmult_strassen_s(MatPak{scratch1, 0,0,n2,n2}, a.subPak(n2,n2,n,n), MatPak{M4, 0, 0, n2, n2});

    madd_s(a.subPak( 0, 0, n2, n2), a.subPak( n2, 0, n, n2), MatPak{scratch1, 0, 0, n2, n2});
    mmult_strassen_s(MatPak{scratch1, 0,0,n2,n2}, b.subPak(n2,n2,n,n), MatPak{M5, 0, 0, M5->sz, M5->sz});

    msub_s(a.subPak( 0, n2, n2, n), a.subPak( 0,0,n2,n2), MatPak{scratch1, 0, 0, n2, n2});
    madd_s(b.subPak( 0,0,n2,n2), b.subPak( n2,0,n,n2), MatPak{scratch2, 0, 0, n2, n2});
    mmult_strassen(scratch1, scratch2, M6);

    msub_s(a.subPak( n2, 0, n, n2), a.subPak( n2,n2,n,n), MatPak{scratch1, 0, 0, n2, n2});
    madd_s(b.subPak(  0, n2, n2, n), b.subPak( n2,n2,n,n), MatPak{scratch2, 0, 0, n2, n2});
    mmult_strassen(scratch1, scratch2, M7);


    // Put together result

    madd_s(MatPak{M1, 0,0,n2,n2}, MatPak{M4, 0,0,n2,n2}, MatPak{res.m, 0, 0, n2, n2});
    msub_s(MatPak{res.m, 0, 0, n2, n2}, MatPak{M5, 0,0,n2,n2}, MatPak{res.m, 0, 0, n2, n2});
    madd_s(MatPak{res.m, 0, 0, n2, n2}, MatPak{M7, 0,0,n2,n2}, MatPak{res.m, 0, 0, n2, n2});

    madd_s(MatPak{M3, 0,0,n2,n2}, MatPak{M5, 0,0,n2,n2}, MatPak{res.m, n2, 0, n, n2});

    madd_s(MatPak{M2, 0,0,n2,n2}, MatPak{M4, 0,0,n2,n2}, MatPak{res.m, 0,n2, n2, n});


    msub_s(MatPak{M1, 0,0,n2,n2}, MatPak{M2, 0,0,n2,n2}, MatPak{res.m, n2, n2, n, n});
    madd_s(MatPak{res.m, n2, n2, n, n}, MatPak{M3, 0,0,n2,n2}, MatPak{res.m, n2, n2, n, n});
    madd_s(MatPak{res.m, n2, n2, n, n}, MatPak{M6, 0,0,n2,n2}, MatPak{res.m, n2, n2, n, n});

    // M1->print();
    // M2->print();
    // M3->print();
    // M4->print();
    // M5->print();
    // M6->print();
    // M7->print();


    delete M1;
    delete M2;
    delete M3;
    delete M4;
    delete M5;
    delete M6;
    delete M7;

    delete scratch1;
    delete scratch2;

}

int main() {
    Matrix* a = new Matrix(3);
    a->index(0,1) = 1;
    a->index(1,0) = 1;
    a->index(2,2) = 1;


    a->print();

    Matrix* b = new Matrix(3);
    b->index(0,0) = 2;
    b->index(1,1) = 2;
    b->index(2,2) = 2;
    b->print();


    Matrix* res = new Matrix(3);

    mmult_strassen(a,b,res);
    res->print();

    mmult(a,b,res);
    res ->print();
}