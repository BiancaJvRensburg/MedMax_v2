#ifndef POINT3_H__
#define POINT3_H__

#include <cassert>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>


template< typename T >
class point3
{
public:
    typedef T               type_t;

    point3< T >( T x_ , T y_ , T z_) { v[0] = x_; v[1] = y_; v[2] = z_; }

    template< class point_t > point3< T >( point_t const & p )
    {
        v[0] = (T)(p[0]);
        v[1] = (T)(p[1]);
        v[2] = (T)(p[2]);
    }

    point3< T >(){v[0] = 0; v[1] = 0; v[2] = 0;}

    inline  T x() const {return v[0];}
    inline  T y() const {return v[1];}
    inline  T z() const {return v[2];}

    inline  T operator [] (unsigned int c) const
    {
        return v[c];
    }
    inline  T & operator [] (unsigned int c)
    {
        return v[c];
    }

    static point3<T> Zero() { return point3<T>(0,0,0); }

    void setZero()
    {
        v[0] = 0;
        v[1] = 0;
        v[2] = 0;
    }

    // You cannot template over anything here, but maybe you could template over typename T2 for operator += (const point3< T2 > & other)
    void operator += (const point3< T > & other)
    {
        v[0] += other.x();
        v[1] += other.y();
        v[2] += other.z();
    }
    void operator -= (const point3< T > & other)
    {
        v[0] -= other.x();
        v[1] -= other.y();
        v[2] -= other.z();
    }

    // This is going to create problems if the compiler needs to resolve umbiguous casts, but it's the cleaner way to do it
    void operator *= (int s)
    {
        v[0] *= s;
        v[1] *= s;
        v[2] *= s;
    }
    void operator *= (unsigned int s)
    {
        v[0] *= s;
        v[1] *= s;
        v[2] *= s;
    }
    void operator *= (float s)
    {
        v[0] *= s;
        v[1] *= s;
        v[2] *= s;
    }
    void operator *= (double s)
    {
        v[0] *= s;
        v[1] *= s;
        v[2] *= s;
    }
    void operator /= (int s)
    {
        v[0] /= s;
        v[1] /= s;
        v[2] /= s;
    }
    void operator /= (unsigned int s)
    {
        v[0] /= s;
        v[1] /= s;
        v[2] /= s;
    }
    void operator /= (float s)
    {
        v[0] /= s;
        v[1] /= s;
        v[2] /= s;
    }
    void operator /= (double s)
    {
        v[0] /= s;
        v[1] /= s;
        v[2] /= s;
    }

    point3<T> getOrthogonal() const
    {
        if( v[0] == 0 )
        {
            return point3<T>( 0 , v[2] , -v[1] );
        }
        else if( v[1] == 0 )
        {
            return point3<T>( v[2] , 0 , -v[0] );
        }

        return point3<T>( v[1] , -v[0] , 0 );
    }

    T norm() const
    {
        return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    }

    void normalize()
    {
        T _n = norm();
        v[0] /= _n;
        v[1] /= _n;
        v[2] /= _n;
    }
    void setNorm( T newNorm )
    {
        T _factor = newNorm / norm();
        v[0] *= _factor;
        v[1] *= _factor;
        v[2] *= _factor;
    }

    inline static
    T dot( const point3< T > & p1 , const point3< T > & p2 )
    {
        return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
    }

    inline static
    point3< T > cross( const point3< T > & p1 , const point3< T > & p2 )
    {
        return point3< T >(
                    p1.y() * p2.z() - p1.z() * p2.y(),
                    p1.z() * p2.x() - p1.x() * p2.z(),
                    p1.x() * p2.y() - p1.y() * p2.x()
                    );
    }

    point3<T> direction() const {
        T n = norm();
        return (*this/n);
    }

    static
    point3< T > rotateVectorSimilarly( point3< T > const & vectorToRotate , point3<T> const & srcVector , point3<T> const & srcVectorRotated ) {
        point3<T> const & axis = point3<T>::cross(srcVector , srcVectorRotated);
        if( sqrt(axis.norm()) < 0.000000001 ) return vectorToRotate; // if there is no rotation. careful with these epsilons...
        point3<T> const & v1 = srcVector.direction();
        point3<T> const & v2 = srcVectorRotated.direction();
        point3<T> const & a = axis.direction();
        point3<T> const & b1 = point3<T>::cross(a , v1);
        point3<T> const & b2 = point3<T>::cross(a , v2);

        return (point3<T>::dot( v1 , vectorToRotate ))*v2 + (point3<T>::dot( a , vectorToRotate ))*a + (point3<T>::dot( b1 , vectorToRotate ))*b2 ;
    }


    bool isnan() const {
        return std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2]);
    }


private:
    T v[3];
};






template< typename T > inline
point3< T > operator + (const point3< T > & p1 , const point3< T > & p2 )
{
    return point3< T >( p1.x() + p2.x() , p1.y() + p2.y() , p1.z() + p2.z() );
}
template< typename T > inline
point3< T > operator - (const point3< T > & p1 , const point3< T > & p2 )
{
    return point3< T >( p1.x() - p2.x() , p1.y() - p2.y() , p1.z() - p2.z() );
}


template< typename T > inline
point3< T > operator - (const point3< T > & p2 )
{
    return point3< T >( - p2.x() , - p2.y() , - p2.z() );
}

template< typename T > inline
point3< T > operator * (const point3< T > & p , int s)
{
    return point3< T >( s*p.x() , s*p.y()  , s*p.z() );
}
template< typename T > inline
point3< T > operator * (const point3< T > & p , float s)
{
    return point3< T >( s*p.x() , s*p.y() , s*p.z()  );
}
template< typename T > inline
point3< T > operator * (const point3< T > & p , double s)
{
    return point3< T >( s*p.x() , s*p.y()  , s*p.z() );
}
template< typename T > inline
point3< T > operator * ( int s , const point3< T > & p )
{
    return point3< T >( s*p.x() , s*p.y()  , s*p.z() );
}
template< typename T > inline
point3< T > operator * ( float s , const point3< T > & p )
{
    return point3< T >( s*p.x() , s*p.y() , s*p.z()  );
}
template< typename T > inline
point3< T > operator * ( double s , const point3< T > & p )
{
    return point3< T >( s*p.x() , s*p.y()  , s*p.z() );
}


template< typename T > inline
point3< T > operator / (const point3< T > & p , int s)
{
    return point3< T >( p.x()/s , p.y()/s  , p.z()/s );
}
template< typename T > inline
point3< T > operator / (const point3< T > & p , float s)
{
    return point3< T >( p.x()/s , p.y()/s  , p.z()/s );
}
template< typename T > inline
point3< T > operator / (const point3< T > & p , double s)
{
    return point3< T >( p.x()/s , p.y()/s  , p.z()/s );
}


template< typename T > inline
T operator * (const point3< T > & p1 , const point3< T > & p2)
{
    return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}


template< typename T > inline
point3<T> operator % (const point3<T> & p1 , const point3<T> & p2)
{
    return point3<T>::cross(p1,p2);
}




template< typename T > inline std::ostream & operator << (std::ostream & s , point3< T > const & p)
{
    s << p[0] << " \t" << p[1] << " \t" << p[2];
    return s;
}





typedef point3< float >    point3f;
typedef point3< double >   point3d;


#endif
