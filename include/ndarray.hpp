#ifndef NDARRAY_H_INCLUDED
#define NDARRAY_H_INCLUDED

#include <cstddef>
#include <stdexcept>


/* ****************************************************************************

    Templated version of a 2/3 D array implementation

    Typical use:

    NDarray<double> array(3,3);

    which creates a 3x3 array of doubles.

   ****************************************************************************/

template <typename Type>
class NDarray
{
public:
    NDarray();
    NDarray(size_t nx);
    NDarray(size_t ny, size_t nx);
    NDarray(size_t nz, size_t ny, size_t nx);

    size_t size() const;
    size_t size(size_t dimension) const;
    size_t get_dimensions() const;
    void resize(size_t nx);
    void resize(size_t ny, size_t nx);
    void resize(size_t nz, size_t ny, size_t nx);

    Type& operator() (size_t i);        // Subscript operators sometimes comes as singleton
    Type  operator() (size_t i) const;  // Subscript operators sometimes comes as singleton
    Type& operator() (size_t j, size_t i);        // Subscript operators often come in pairs
    Type  operator() (size_t j, size_t i) const;  // Subscript operators often come in pairs
    Type& operator() (size_t k, size_t j, size_t i);        // or in triplets
    Type  operator() (size_t k, size_t j, size_t i) const;  // or in triplets

    ~NDarray();                              // Destructor
    NDarray(const NDarray& m);               // Copy constructor
    NDarray& operator= (const NDarray& m);   // Assignment operator
    NDarray& operator= (const Type* pvalues);// Assignment operator for pointer support.

private:
    size_t nz_, ny_, nx_, dimensions_;
    Type* data_;
};


// Constructors
template <typename Type>
NDarray<Type>::NDarray() : nz_{0}, ny_ {0}, nx_ {0}, dimensions_ {0}, data_ {new Type[1]}
{
    //data_ = new Type[1];
}

template <typename Type>
NDarray<Type>::NDarray(size_t nx) : nz_{0}, ny_{0}, nx_ (nx), dimensions_ {1}, data_ {new Type[nx]}
{
    if (nx == 0)
        throw std::runtime_error("NDarray constructor has 0 size");
    //data_ = new Type[nx];
}

template <typename Type>
NDarray<Type>::NDarray(size_t ny, size_t nx) : nz_{0}, ny_ (ny), nx_ (nx), dimensions_ (2), data_ {new Type[ny * nx]}
{
    if (ny == 0 || nx == 0)
        throw std::runtime_error("NDarray constructor has 0 size");
    //data_ = new Type[ny * nx];
}

template <typename Type>
NDarray<Type>::NDarray(size_t nz, size_t ny, size_t nx) : nz_ (nz), ny_ (ny), nx_ (nx), dimensions_ (3), data_ {new Type[nz*ny*nx]}
{
    if (nz ==0 || ny == 0 || nx == 0)
        throw std::runtime_error("NDarray constructor has 0 size");
    //data_ = new Type[nz * ny * nx];
}


template <typename Type>
size_t NDarray<Type>::get_dimensions() const
{
    return dimensions_;
}

template <typename Type>
size_t NDarray<Type>::size(size_t dimension) const
{
    if (dimension>dimensions_)
        throw std::runtime_error("Object has not sufficient dimensions.");
    size_t s;
    if (dimension == 1)
        s = nx_;
    else if (dimension == 2)
        s = ny_;
    else if (dimension == 3)
        s = nz_;
    else
        throw std::runtime_error("Only 1, 2 and 3 dimension arrays are supported.");
    return s;
}

template <typename Type>
size_t NDarray<Type>::size() const
{
    size_t s;
    if (dimensions_ == 1)
        s = nx_;
    else if (dimensions_ == 2)
        s = ny_ * nx_;
    else if (dimensions_ == 3)
        s = nz_ * ny_ * nx_;
    else
        throw std::runtime_error("Only 1, 2 and 3 dimension arrays are supported.");
    return s;
}


template <typename Type>
void NDarray<Type>::resize(size_t nx)
{
    dimensions_ = 1;
    nx_ = nx;
    delete [] data_;
    data_ = new Type[nx_];
}


template <typename Type>
void NDarray<Type>::resize(size_t ny, size_t nx)
{
    dimensions_ = 2;
    ny_ = ny;
    nx_ = nx;
    delete [] data_;
    data_ = new Type[ny_ * nx_];
}

template <typename Type>
void NDarray<Type>::resize(size_t nz, size_t ny, size_t nx)
{
    dimensions_ = 3;
    nz_ = nz;
    ny_ = ny;
    nx_ = nx;
    delete [] data_;
    data_ = new Type[nz_ * ny_ * nx_];
}


template <typename Type>
NDarray<Type>::NDarray(const NDarray<Type>& m) : nz_ {}, ny_ {}, nx_ {}, dimensions_ {}, data_ {}// Copy constructor
{
    size_t n;
    nx_ = m.nx_;
    ny_ = m.ny_;
    nz_ = m.nz_;
    dimensions_ = m.dimensions_;
    if (dimensions_==2)
        n = ny_ * nx_;
    else
        n = nz_ * ny_ * nx_;
    data_ = new Type[n];
    for(size_t i=0; i<n; ++i)
    {
        data_[i] = m.data_[i];
    }
}

template <typename Type>
NDarray<Type>& NDarray<Type>::operator= (const NDarray<Type>& m)   // Assignment operator
{
    size_t n;
    nx_ = m.nx_;
    ny_ = m.ny_;
    nz_ = m.nz_;
    dimensions_ = m.dimensions_;
    if (dimensions_==2)
        n = ny_ * nx_;
    else
        n = nz_ * ny_ * nx_;
    delete [] data_; // delete any existing data.
    data_ = new Type[n];
    for(size_t i=0; i<n; ++i)
    {
        data_[i] = m.data_[i];
    }
    return *this;
}

template <typename Type>
NDarray<Type>& NDarray<Type>::operator= (const Type* pvalues)   // Assignment operator for pointer support
{
    size_t n;
    if (dimensions_==2)
        n = ny_ * nx_;
    else
        n = nz_ * ny_ * nx_;
    delete [] data_; // delete any existing data.
    data_ = new Type[n];
    for(size_t i=0; i<n; ++i)
    {
        data_[i] = *pvalues;
        pvalues++;
    }
    return *this;
}

// Destructor

template <typename Type>
NDarray<Type>::~NDarray()
{
      delete[] data_;
}

template <typename Type>
Type& NDarray<Type>::operator() (size_t i)
{
    if (dimensions_ == 1)
    {
        if (i >= nx_)
            throw std::runtime_error("NDarray subscript out of bounds");
    }
    else if (dimensions_ == 2)
    {
        if (i >= nx_*ny_)
            throw std::runtime_error("NDarray subscript out of bounds");
    }
    else if (dimensions_ == 3)
    {
        if (i >= nx_*ny_*nz_)
            throw std::runtime_error("NDarray subscript out of bounds");
    }
    return data_[i];
}


template <typename Type>
Type NDarray<Type>::operator() (size_t i) const
{
    if (dimensions_ == 1)
    {
        if (i >= nx_)
            throw std::runtime_error("NDarray subscript out of bounds");
    }
    else if (dimensions_ == 2)
    {
        if (i >= nx_*ny_)
            throw std::runtime_error("NDarray subscript out of bounds");
    }
    else if (dimensions_ == 3)
    {
        if (i >= nx_*ny_*nz_)
            throw std::runtime_error("NDarray subscript out of bounds");
    }
    return data_[i];
}

template <typename Type>
Type& NDarray<Type>::operator() (size_t j, size_t i)
{
    if (dimensions_ != 2)
        throw std::runtime_error("Dimension error. This instance requires exactly two arguments.");
    if (j >= ny_ || i >= nx_)
        throw std::runtime_error("NDarray subscript out of bounds");
    return data_[nx_*j + i];
}


template <typename Type>
Type NDarray<Type>::operator() (size_t j, size_t i) const
{
    if (dimensions_ != 2)
        throw std::runtime_error("Dimension error. This instance requires exactly two arguments.");
    if (j >= ny_ || i >= nx_)
        throw std::runtime_error("const NDarray subscript out of bounds");
    return data_[nx_*j + i];
}


template <typename Type>
Type& NDarray<Type>::operator() (size_t k, size_t j, size_t i)
{
    if (dimensions_ != 3)
        throw std::runtime_error("Dimension error. This instance requires exactly three arguments.");
    if (k>=nz_ || j >= ny_ || i >= nx_)
        throw std::runtime_error("NDarray subscript out of bounds");
    return data_[(ny_*nx_)*k + nx_*j + i];
}


template <typename Type>
Type NDarray<Type>::operator() (size_t k, size_t j, size_t i) const
{
    if (dimensions_ != 3)
        throw std::runtime_error("Dimension error. This instance requires exactly three arguments.");
    if (k >= nz_ || j >= ny_ || i >= nx_)
        throw std::runtime_error("NDarray subscript out of bounds");
    return data_[(ny_*nx_)*k + nx_*j + i];
}


#endif // NDARRAY_H_INCLUDED
