//
// Created by janos on 12/14/20.
//

#pragma once

#include <cholmod.h>

class Triplets;

class SparseMatrix {
public:
    SparseMatrix(size_t m = 0, size_t n = 0, size_t nnz = 0);

    explicit SparseMatrix(Triplets& T);

    SparseMatrix(SparseMatrix&& B);

    SparseMatrix& operator=(SparseMatrix&& B);

    ~SparseMatrix();

    static SparseMatrix identity(size_t m, size_t n);

    SparseMatrix transpose() const;

    size_t nRows() const;

    size_t nCols() const;

    size_t nnz() const;

    // returns norm. 0: Infinity, 1: 1-norm
    double norm(int norm) const;

    cholmod_sparse* toCholmod();

    // math
    friend SparseMatrix operator*(const SparseMatrix& A, double s);
    friend SparseMatrix operator+(const SparseMatrix& A, const SparseMatrix& B);
    friend SparseMatrix operator-(const SparseMatrix& A, const SparseMatrix& B);
    friend SparseMatrix operator*(const SparseMatrix& A, const SparseMatrix& B);

    friend SparseMatrix& operator*=(SparseMatrix& A, double s);
    friend SparseMatrix& operator+=(SparseMatrix& A, const SparseMatrix& B);
    friend SparseMatrix& operator-=(SparseMatrix& A, const SparseMatrix& B);

protected:
    // clear
    void clear();

    // member
    cholmod_sparse *data;
};

class Triplets {
public:
    // constructor
    Triplets(size_t m, size_t n);

    // destructor
    ~Triplets();

    // add entry
    void add(size_t i, size_t j, double x);

    // returns choldmod representation
    cholmod_triplet* toCholmod();

protected:
    // increases capacity
    void increaseCapacity();

    // clear
    void clear();

    // member
    cholmod_triplet *data;
    size_t m, n, capacity;
};