#pragma once

#include "myutils.h"

template<class T> class span
    {
private:
    T *m_data;
    size_t m_size;

public:
    span() : m_data(nullptr), m_size(0) {}
    span(T *data, size_t size) : m_data(data), m_size(size) {}
    T *data() const { return m_data; }
    size_t size() const { return m_size; }
    T &operator[](size_t index) const { return m_data[(size_t) index]; }
    };


template<class T> class Matrix {
    private:
    T *m_data;
    size_t m_rows;
    size_t m_cols;
    bool m_OwnData = false;

public:
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
    
    static Matrix<T> Allocate(uint rows, uint cols) { return Matrix<T>(myalloc64(T, rows*cols), rows, cols, true); }
    static Matrix<T> FromOwnedData(T *data, uint rows, uint cols) { return Matrix<T>(data, rows, cols, true); }
    static Matrix<T> FromSharedData(T *data, uint rows, uint cols) { return Matrix<T>(data, rows, cols, false); }

    T *data() const { return m_data; }
    uint Rows() const { return m_rows; }
    uint Cols() const { return m_cols; }

    span<T> operator[](uint row) { return span<T>(m_data + (size_t) row * (size_t) m_cols, m_cols); }
    span<const T> operator[](uint row) const { return span<const T>(m_data + (size_t) row * (size_t) m_cols, m_cols); }

    ~Matrix() { if (m_OwnData) myfree(m_data); }
    };

template<class T> class SquareMatrix {
        private:
        T *m_data;
        size_t m_size;
        bool m_OwnData = false;

    public:
        SquareMatrix() : m_data(nullptr), m_size(0), m_OwnData(false) {}

        SquareMatrix(T *data, uint size, bool ownData) : m_data(data), m_size(size), m_OwnData(ownData) {}

        SquareMatrix(const SquareMatrix& other) = delete;

        //Matrix(const Matrix& other) : m_data(other.m_data), m_rows(other.m_rows), m_cols(other.m_cols), m_OwnData(false) {}
        SquareMatrix(SquareMatrix&& other) noexcept : m_data(other.m_data), m_size(other.m_size), m_OwnData(other.m_OwnData)
            {
            other.m_data = nullptr;
            other.m_size = 0;
            other.m_OwnData = false;
            }

        SquareMatrix& operator=(const SquareMatrix& other) = delete;
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

        SquareMatrix& operator=(SquareMatrix&& other) noexcept
            {
            if (this != &other)
                {
                if (m_OwnData) myfree(m_data);
                m_data = other.m_data;
                m_size = other.m_size;
                m_OwnData = other.m_OwnData;
                other.m_data = nullptr;
                other.m_size = 0;
                other.m_OwnData = false;
                }
            return *this;
            }

        static SquareMatrix<T> Allocate(uint size) { return SquareMatrix<T>(myalloc64(T, size*size), size, true); }
        static SquareMatrix<T> FromOwnedData(T *data, uint size) { return SquareMatrix<T>(data, size, true); }
        static SquareMatrix<T> FromSharedData(T *data, uint size) { return SquareMatrix<T>(data, size, false); }

        T *data() const { return m_data; }
        uint Rows() const { return m_size; }
        uint Cols() const { return m_size; }

        span<T> operator[](uint row) { return span<T>(m_data + (size_t) row * (size_t) m_size, m_size); }
        span<const T> operator[](uint row) const { return span<const T>(m_data + (size_t) row * (size_t) m_size, m_size); }

        ~SquareMatrix() { if (m_OwnData) myfree((void*) m_data); }
    };