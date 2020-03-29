#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <string.h> 
#include <fstream>
#include <assert.h>
#include <chrono>

size_t CUTOFF = 3;

struct Matrix{
    size_t sz;
    long* arr; 

    Matrix(size_t sz) {
        this->sz = sz;
        this->arr = new long[sz * sz];
        this->clear();
    }

    inline long& index(size_t x, size_t y) {
        return this->arr[sz*x + y];
    }

    inline void clear() {
        memset(arr, 0, sz*sz*sizeof(long));
    }

    ~Matrix() {
        delete arr;
    }

    void print() {
        printf("[\n");
        for(int i = 0; i < sz; i++) {
            printf("   [");
            for(int j = 0; j < sz; j++) {
                printf("%ld,", this->index(i,j));
            }
            printf("]\n");
        }
        printf("]\n");
    }
};

struct MatPak;

MatPak make(Matrix* m,size_t s_x,size_t s_y, size_t e_x, size_t e_y);

long nothing = 0;

struct MatPak {
    Matrix* m;
    size_t s_x;
    size_t s_y;

    size_t e_x;
    size_t e_y;

    bool padded_r;
    bool padded_d;

    inline long& index(size_t x, size_t y) {
        return m->index(s_x + x, s_y + y);
    };

    inline MatPak subPak(size_t s_x,size_t s_y, size_t e_x, size_t e_y) {
        return make(m, this->s_x + s_x, this->s_y + s_y, this->s_x + e_x, this->s_y + e_y);
    };

    inline MatPak subSet(size_t i, size_t j) {
        size_t n = this->e_x - this->s_x;
        size_t s_x0 = (i == 1) ? 0 : (n+1)/2;
        size_t s_y0 = (j == 1) ? 0 : (n+1)/2;

        size_t e_x0 = s_x0 + (n+1)/2;
        size_t e_y0 = s_y0 + (n+1)/2;


        printf("(%ld, %ld), (%ld, %ld)\n", s_x0, s_y0, e_x0, e_y0);
        return subPak(s_x0, s_y0, e_x0, e_y0);
    }

};

// Always returns an even-size MatPak, padded with 0 if bounds given extend past edge
MatPak make(Matrix* m,size_t s_x,size_t s_y, size_t e_x, size_t e_y) {
    bool padded_r = false;
    bool padded_d = false;
    size_t n_x = e_x - s_x;
    if (e_x > m->sz) {
        padded_r = true;
    }

    size_t n_y = e_y - s_y;
    if (e_y > m->sz) {
        padded_d = true;
    }


    return MatPak{m, s_x, s_y, e_x, e_y, padded_r, padded_d};

}

// MatPak submat(MatPak)


void mmult_s(MatPak a, MatPak b, MatPak res) {
    size_t n = a.e_x - a.s_x;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {

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
    assert(a.e_y - a.s_y == b.e_y - b.s_y);
    assert(a.e_x - a.s_x == b.e_x - b.s_x);

    assert(a.e_y - a.s_y == res.e_y - res.s_y);
    assert(a.e_x - a.s_x == res.e_x - res.s_x);


    size_t n = a.e_x - a.s_x;

    for(int i = 0; i < a.e_x - a.s_x; i++) {
        for (int j = 0; j < a.e_y - a.s_y; j++) {
            res.index(i, j) 
                = a.index(i, j)
                + b.index(i, j);
        }
    }

}

void msub_s(MatPak a, MatPak b, MatPak res) {
    // assert(a.m->sz >= a.e_x);
    // assert(a.m->sz >= a.e_y);

    // assert(b.m->sz >= b.e_x);
    // assert(b.m->sz >= b.e_y);

    assert(a.e_y - a.s_y == b.e_y - b.s_y);
    assert(a.e_x - a.s_x == b.e_x - b.s_x);

    for(int i = 0; i < a.e_x - a.s_x; i++) {
        for (int j = 0; j < a.e_y - a.s_y; j++) {
            res.index(i, j) 
                = a.index(i, j)
                - b.index(i, j);
        }
    }

}

// This is OK in place
void madd(Matrix* a, Matrix* b, Matrix* res) {
    madd_s(MatPak{a, 0, 0, a->sz, a->sz}, MatPak{b, 0, 0, b->sz, b->sz}, MatPak{res, 0, 0, res->sz, res->sz});
}

void cmult(Matrix* a, long c) {
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
    size_t n2 = (a.e_x - a.s_x + 1) / 2;
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

    MatPak scratch1_p = make(scratch1, 0, 0, n2, n2);
    MatPak scratch2_p = make(scratch2, 0, 0, n2, n2);

    // std::cout << n << "hither" << n2 << "\n" << std::flush;

    madd_s(a.subSet(1,1), a.subSet(2,2), scratch1_p);
    madd_s(b.subSet(1,1), b.subSet(2,2), scratch2_p);
    mmult_strassen(scratch1, scratch2, M1);


    madd_s(a.subSet(2,1), a.subSet(2,2), scratch1_p);
    mmult_strassen_s(scratch1_p, b.subSet(1,1), make(M2, 0, 0, M2->sz, M2->sz));

    msub_s(b.subSet(1,2), b.subSet(2,2), scratch1_p);
    mmult_strassen_s(a.subSet(1,1), scratch1_p, make(M3, 0, 0, n2, n2));

    msub_s(b.subSet(2, 1), b.subSet(1, 1), scratch1_p);
    mmult_strassen_s(a.subSet(2,2), scratch1_p, make(M4, 0, 0, n2, n2));

    madd_s(a.subSet(1,1),a.subSet(1,2), scratch1_p);
    mmult_strassen_s(scratch1_p, b.subSet(2,2), make(M5, 0, 0, n2, n2));

    msub_s(a.subSet(2,1), a.subSet(1,1), scratch1_p);
    madd_s(b.subSet(1,1), b.subSet(1,2), scratch2_p);
    mmult_strassen(scratch1, scratch2, M6);

    msub_s(a.subSet(1,2), a.subSet(2,2), scratch1_p);
    madd_s(b.subSet(2,1), b.subSet(2,2), scratch2_p);
    mmult_strassen(scratch1, scratch2, M7);

    M7->print();

    M1->print();
    M2->print();
    M3->print();
    M4->print();
    M5->print();
    M6->print();
    M7->print();


    // Put together result
    // Have to be more careful here, to keep bounds right!

    madd_s(make(M1, 0,0,n2,n2), make(M4, 0,0,n2,n2), res.subSet(1,1));
    msub_s(res.subSet(1,1), make(M5, 0,0,n2,n2), res.subSet(1,1));
    madd_s(res.subSet(1,1), make(M7, 0,0,n2,n2), res.subSet(1,1));

    madd_s(make(M3, 0,0,n2,n2), make(M5, 0,0,n2,n2), res.subSet(1,2));

    madd_s(make(M2, 0,0,n2,n2), make(M4, 0,0,n2,n2), res.subSet(2,1));


    msub_s(make(M1, 0,0,n2,n2), make(M2, 0,0,n2,n2), res.subSet(2,2));
    madd_s(res.subSet(2,2), make(M3, 0,0,n2,n2), res.subSet(2,2));
    madd_s(res.subSet(2,2), make(M6, 0,0,n2,n2), res.subSet(2,2));

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

bool are_equal(Matrix* a, Matrix* b) {
    for (size_t i = 0; i < a->sz*a->sz; i++) {
        if(a->arr[i] != b->arr[i]) {
            return false;
        }
    }
    return true;
}

int main(int argc, char** argv) {

    size_t cutoff_a = strtoul(argv[2],nullptr, 10);
    if(cutoff_a != 0) {
        CUTOFF = cutoff_a;
    }

    // if(strtoul(argv[2],nullptr, 10))

    size_t dim = strtoul(argv[2],nullptr, 10);

    std::ifstream infile; 

    infile.open(argv[3]);

    if (!infile) {
        printf("File N/A\n");
        exit(1);
    }

    Matrix* a = new Matrix(dim);
    Matrix* b = new Matrix(dim);

    for(size_t i = 0; i < dim * dim; i++) {
        infile >> a->arr[i];
    }

    for(size_t i = 0; i < dim * dim; i++) {
        infile >> b->arr[i];
    }

    infile.close();

    // a->print();
    // b->print();



    printf("Multiplying strassen...\n");

    unsigned long total_i = 0;
    unsigned long total_s = 0;

    for(int i = 0; i < 1; i++) {
        // {
        //     Matrix* res = new Matrix(dim);
        //     auto begin = std::chrono::high_resolution_clock::now();
        //     mmult(a,b,res);
        //     auto end = std::chrono::high_resolution_clock::now();
        //     total_i += std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
        //     delete res;
        // }
        {
            Matrix* res = new Matrix(dim);
            auto begin = std::chrono::high_resolution_clock::now();
            mmult_strassen(a,b,res);
            auto end = std::chrono::high_resolution_clock::now();
            total_s += std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
            delete res;
        }
    }

    // printf("Iterative: %lu μs\n", total_i/10);

    printf("Strassen: %lu μs\n", total_s/1);


    // assert(are_equal(res, res_s));
}