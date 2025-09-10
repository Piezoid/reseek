#pragma once

#include "myutils.h"

template<class T> class span
    {
private:
    T *m_data;
    uint m_size;

public:
    span(T *data, uint size) : m_data(data), m_size(size) {}
    T *data() const { return m_data; }
    uint size() const { return m_size; }
    T &operator[](uint index) const { return m_data[index]; }
    };


template<class T> class Matrix {
    private:
    T *m_data;
    uint m_rows;
    uint m_cols;
    bool m_OwnData = false;

public:
    // Default constructor - creates empty matrix
    Matrix() : m_data(nullptr), m_rows(0), m_cols(0), m_OwnData(false) {}
    
    Matrix(T *data, uint rows, uint cols, bool ownData) : m_data(data), m_rows(rows), m_cols(cols), m_OwnData(ownData) {}
    
    Matrix(const Matrix& other) = delete;

    //Matrix(const Matrix& other) : m_data(other.m_data), m_rows(other.m_rows), m_cols(other.m_cols), m_OwnData(false) {}
    
    Matrix(Matrix&& other) noexcept : m_data(other.m_data), m_rows(other.m_rows), m_cols(other.m_cols), m_OwnData(other.m_OwnData) 
        {
        other.m_data = nullptr;
        other.m_rows = 0;
        other.m_cols = 0;
        other.m_OwnData = false;
        }
    
    Matrix& operator=(const Matrix& other) = delete;
    // Matrix& operator=(const Matrix& other) 
    //     {
    //     if (this != &other) 
    //         {
    //         if (m_OwnData) myfree(m_data);
    //         m_data = other.m_data;
    //         m_rows = other.m_rows;
    //         m_cols = other.m_cols;
    //         m_OwnData = false;
    //         }
    //     return *this;
    //     }
    
    Matrix& operator=(Matrix&& other) noexcept 
        {
        if (this != &other) 
            {
            if (m_OwnData) myfree(m_data);
            m_data = other.m_data;
            m_rows = other.m_rows;
            m_cols = other.m_cols;
            m_OwnData = other.m_OwnData;
            other.m_data = nullptr;
            other.m_rows = 0;
            other.m_cols = 0;
            other.m_OwnData = false;
            }
        return *this;
        }
    
    static Matrix<T> Allocate(uint rows, uint cols) { return Matrix<T>(myalloc(T, rows*cols), rows, cols, false); }
    static Matrix<T> FromOwnedData(T *data, uint rows, uint cols) { return Matrix<T>(data, rows, cols, true); }

    T *data() const { return m_data; }
    uint Rows() const { return m_rows; }
    uint Cols() const { return m_cols; }

    span<T> operator[](uint row) { return span<T>(m_data + row*m_cols, m_cols); }
    span<const T> operator[](uint row) const { return span<const T>(m_data + row*m_cols, m_cols); }

    ~Matrix() { if (m_OwnData) myfree(m_data); }
    };

// Conversion function from vector<vector<byte>> to Matrix<byte>
template<class T>
Matrix<T> VectorToMatrix(const vector<vector<T>> &vec) {
    if (vec.empty()) {
        return Matrix<T>();
    }
    
    uint rows = vec.size();
    uint cols = vec[0].size();
    T *data = myalloc(T, rows * cols);
    
    for (uint i = 0; i < rows; ++i) {
        for (uint j = 0; j < cols; ++j) {
            data[i * cols + j] = vec[i][j];
        }
    }
    
    return Matrix<T>::FromOwnedData(data, rows, cols);
}

// Conversion function from Matrix<byte> to vector<vector<byte>>
template<class T>
vector<vector<T>> MatrixToVector(const Matrix<T> &matrix) {
    vector<vector<T>> result;
    if (matrix.Rows() == 0 || matrix.Cols() == 0) {
        return result;
    }
    
    result.resize(matrix.Rows());
    for (uint i = 0; i < matrix.Rows(); ++i) {
        result[i].resize(matrix.Cols());
        for (uint j = 0; j < matrix.Cols(); ++j) {
            result[i][j] = matrix[i][j];
        }
    }
    
    return result;
}


