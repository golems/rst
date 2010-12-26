// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2006-2008 Benoit Jacob <jacob.benoit.1@gmail.com>
// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#ifndef EIGEN_PRODUCT_H
#define EIGEN_PRODUCT_H

/** \class GeneralProduct
  * \ingroup Core_Module
  *
  * \brief Expression of the product of two general matrices or vectors
  *
  * \param LhsNested the type used to store the left-hand side
  * \param RhsNested the type used to store the right-hand side
  * \param ProductMode the type of the product
  *
  * This class represents an expression of the product of two general matrices.
  * We call a general matrix, a dense matrix with full storage. For instance,
  * This excludes triangular, selfadjoint, and sparse matrices.
  * It is the return type of the operator* between general matrices. Its template
  * arguments are determined automatically by ProductReturnType. Therefore,
  * GeneralProduct should never be used direclty. To determine the result type of a
  * function which involves a matrix product, use ProductReturnType::Type.
  *
  * \sa ProductReturnType, MatrixBase::operator*(const MatrixBase<OtherDerived>&)
  */
template<typename Lhs, typename Rhs, int ProductType = internal::product_type<Lhs,Rhs>::value>
class GeneralProduct;

enum {
  Large = 2,
  Small = 3
};

namespace internal {

template<int Rows, int Cols, int Depth> struct product_type_selector;

template<int Size, int MaxSize> struct product_size_category
{
  enum { is_large = MaxSize == Dynamic ||
                    Size >= EIGEN_CACHEFRIENDLY_PRODUCT_THRESHOLD,
         value = is_large  ? Large
               : Size == 1 ? 1
                           : Small
  };
};

template<typename Lhs, typename Rhs> struct product_type
{
  typedef typename remove_all<Lhs>::type _Lhs;
  typedef typename remove_all<Rhs>::type _Rhs;
  enum {
    MaxRows  = _Lhs::MaxRowsAtCompileTime,
    Rows  = _Lhs::RowsAtCompileTime,
    MaxCols  = _Rhs::MaxColsAtCompileTime,
    Cols  = _Rhs::ColsAtCompileTime,
    MaxDepth = EIGEN_SIZE_MIN_PREFER_FIXED(_Lhs::MaxColsAtCompileTime,
                                           _Rhs::MaxRowsAtCompileTime),
    Depth = EIGEN_SIZE_MIN_PREFER_FIXED(_Lhs::ColsAtCompileTime,
                                        _Rhs::RowsAtCompileTime),
    LargeThreshold = EIGEN_CACHEFRIENDLY_PRODUCT_THRESHOLD
  };

  // the splitting into different lines of code here, introducing the _select enums and the typedef below,
  // is to work around an internal compiler error with gcc 4.1 and 4.2.
private:
  enum {
    rows_select = product_size_category<Rows,MaxRows>::value,
    cols_select = product_size_category<Cols,MaxCols>::value,
    depth_select = product_size_category<Depth,MaxDepth>::value
  };
  typedef product_type_selector<rows_select, cols_select, depth_select> selector;

public:
  enum {
    value = selector::ret
  };
#ifdef EIGEN_DEBUG_PRODUCT
  static void debug()
  {
      EIGEN_DEBUG_VAR(Rows);
      EIGEN_DEBUG_VAR(Cols);
      EIGEN_DEBUG_VAR(Depth);
      EIGEN_DEBUG_VAR(rows_select);
      EIGEN_DEBUG_VAR(cols_select);
      EIGEN_DEBUG_VAR(depth_select);
      EIGEN_DEBUG_VAR(value);
  }
#endif
};


/* The following allows to select the kind of product at compile time
 * based on the three dimensions of the product.
 * This is a compile time mapping from {1,Small,Large}^3 -> {product types} */
// FIXME I'm not sure the current mapping is the ideal one.
template<int M, int N>  struct product_type_selector<M,N,1>              { enum { ret = OuterProduct }; };
template<int Depth>     struct product_type_selector<1,    1,    Depth>  { enum { ret = InnerProduct }; };
template<>              struct product_type_selector<1,    1,    1>      { enum { ret = InnerProduct }; };
template<>              struct product_type_selector<Small,1,    Small>  { enum { ret = CoeffBasedProductMode }; };
template<>              struct product_type_selector<1,    Small,Small>  { enum { ret = CoeffBasedProductMode }; };
template<>              struct product_type_selector<Small,Small,Small>  { enum { ret = CoeffBasedProductMode }; };
template<>              struct product_type_selector<Small, Small, 1>    { enum { ret = LazyCoeffBasedProductMode }; };
template<>              struct product_type_selector<Small, Large, 1>    { enum { ret = LazyCoeffBasedProductMode }; };
template<>              struct product_type_selector<Large, Small, 1>    { enum { ret = LazyCoeffBasedProductMode }; };
template<>              struct product_type_selector<1,    Large,Small>  { enum { ret = CoeffBasedProductMode }; };
template<>              struct product_type_selector<1,    Large,Large>  { enum { ret = GemvProduct }; };
template<>              struct product_type_selector<1,    Small,Large>  { enum { ret = CoeffBasedProductMode }; };
template<>              struct product_type_selector<Large,1,    Small>  { enum { ret = CoeffBasedProductMode }; };
template<>              struct product_type_selector<Large,1,    Large>  { enum { ret = GemvProduct }; };
template<>              struct product_type_selector<Small,1,    Large>  { enum { ret = CoeffBasedProductMode }; };
template<>              struct product_type_selector<Small,Small,Large>  { enum { ret = GemmProduct }; };
template<>              struct product_type_selector<Large,Small,Large>  { enum { ret = GemmProduct }; };
template<>              struct product_type_selector<Small,Large,Large>  { enum { ret = GemmProduct }; };
template<>              struct product_type_selector<Large,Large,Large>  { enum { ret = GemmProduct }; };
template<>              struct product_type_selector<Large,Small,Small>  { enum { ret = GemmProduct }; };
template<>              struct product_type_selector<Small,Large,Small>  { enum { ret = GemmProduct }; };
template<>              struct product_type_selector<Large,Large,Small>  { enum { ret = GemmProduct }; };

} // end namespace internal

/** \class ProductReturnType
  * \ingroup Core_Module
  *
  * \brief Helper class to get the correct and optimized returned type of operator*
  *
  * \param Lhs the type of the left-hand side
  * \param Rhs the type of the right-hand side
  * \param ProductMode the type of the product (determined automatically by internal::product_mode)
  *
  * This class defines the typename Type representing the optimized product expression
  * between two matrix expressions. In practice, using ProductReturnType<Lhs,Rhs>::Type
  * is the recommended way to define the result type of a function returning an expression
  * which involve a matrix product. The class Product should never be
  * used directly.
  *
  * \sa class Product, MatrixBase::operator*(const MatrixBase<OtherDerived>&)
  */
template<typename Lhs, typename Rhs, int ProductType>
struct ProductReturnType
{
  // TODO use the nested type to reduce instanciations ????
//   typedef typename internal::nested<Lhs,Rhs::ColsAtCompileTime>::type LhsNested;
//   typedef typename internal::nested<Rhs,Lhs::RowsAtCompileTime>::type RhsNested;

  typedef GeneralProduct<Lhs/*Nested*/, Rhs/*Nested*/, ProductType> Type;
};

template<typename Lhs, typename Rhs>
struct ProductReturnType<Lhs,Rhs,CoeffBasedProductMode>
{
  typedef typename internal::nested<Lhs, Rhs::ColsAtCompileTime, typename internal::plain_matrix_type<Lhs>::type >::type LhsNested;
  typedef typename internal::nested<Rhs, Lhs::RowsAtCompileTime, typename internal::plain_matrix_type<Rhs>::type >::type RhsNested;
  typedef CoeffBasedProduct<LhsNested, RhsNested, EvalBeforeAssigningBit | EvalBeforeNestingBit> Type;
};

template<typename Lhs, typename Rhs>
struct ProductReturnType<Lhs,Rhs,LazyCoeffBasedProductMode>
{
  typedef typename internal::nested<Lhs, Rhs::ColsAtCompileTime, typename internal::plain_matrix_type<Lhs>::type >::type LhsNested;
  typedef typename internal::nested<Rhs, Lhs::RowsAtCompileTime, typename internal::plain_matrix_type<Rhs>::type >::type RhsNested;
  typedef CoeffBasedProduct<LhsNested, RhsNested, NestByRefBit> Type;
};

// this is a workaround for sun CC
template<typename Lhs, typename Rhs>
struct LazyProductReturnType : public ProductReturnType<Lhs,Rhs,LazyCoeffBasedProductMode>
{};

/***********************************************************************
*  Implementation of Inner Vector Vector Product
***********************************************************************/

// FIXME : maybe the "inner product" could return a Scalar
// instead of a 1x1 matrix ??
// Pro: more natural for the user
// Cons: this could be a problem if in a meta unrolled algorithm a matrix-matrix
// product ends up to a row-vector times col-vector product... To tackle this use
// case, we could have a specialization for Block<MatrixType,1,1> with: operator=(Scalar x);

namespace internal {

template<typename Lhs, typename Rhs>
struct traits<GeneralProduct<Lhs,Rhs,InnerProduct> >
 : traits<Matrix<typename scalar_product_traits<typename Lhs::Scalar, typename Rhs::Scalar>::ReturnType,1,1> >
{};

}

template<typename Lhs, typename Rhs>
class GeneralProduct<Lhs, Rhs, InnerProduct>
  : internal::no_assignment_operator,
    public Matrix<typename internal::scalar_product_traits<typename Lhs::Scalar, typename Rhs::Scalar>::ReturnType,1,1>
{
    typedef Matrix<typename internal::scalar_product_traits<typename Lhs::Scalar, typename Rhs::Scalar>::ReturnType,1,1> Base;
  public:
    GeneralProduct(const Lhs& lhs, const Rhs& rhs)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<typename Lhs::RealScalar, typename Rhs::RealScalar>::value),
        YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)

      Base::coeffRef(0,0) = (lhs.transpose().cwiseProduct(rhs)).sum();
    }

    /** Convertion to scalar */
    operator const typename Base::Scalar() const {
      return Base::coeff(0,0);
    }
};

/***********************************************************************
*  Implementation of Outer Vector Vector Product
***********************************************************************/

namespace internal {
template<int StorageOrder> struct outer_product_selector;

template<typename Lhs, typename Rhs>
struct traits<GeneralProduct<Lhs,Rhs,OuterProduct> >
 : traits<ProductBase<GeneralProduct<Lhs,Rhs,OuterProduct>, Lhs, Rhs> >
{};

}

template<typename Lhs, typename Rhs>
class GeneralProduct<Lhs, Rhs, OuterProduct>
  : public ProductBase<GeneralProduct<Lhs,Rhs,OuterProduct>, Lhs, Rhs>
{
  public:
    EIGEN_PRODUCT_PUBLIC_INTERFACE(GeneralProduct)

    GeneralProduct(const Lhs& lhs, const Rhs& rhs) : Base(lhs,rhs)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<typename Lhs::RealScalar, typename Rhs::RealScalar>::value),
        YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)
    }

    template<typename Dest> void scaleAndAddTo(Dest& dest, Scalar alpha) const
    {
      internal::outer_product_selector<(int(Dest::Flags)&RowMajorBit) ? RowMajor : ColMajor>::run(*this, dest, alpha);
    }
};

namespace internal {

template<> struct outer_product_selector<ColMajor> {
  template<typename ProductType, typename Dest>
  static EIGEN_DONT_INLINE void run(const ProductType& prod, Dest& dest, typename ProductType::Scalar alpha) {
    typedef typename Dest::Index Index;
    // FIXME make sure lhs is sequentially stored
    // FIXME not very good if rhs is real and lhs complex while alpha is real too
    const Index cols = dest.cols();
    for (Index j=0; j<cols; ++j)
      dest.col(j) += (alpha * prod.rhs().coeff(j)) * prod.lhs();
  }
};

template<> struct outer_product_selector<RowMajor> {
  template<typename ProductType, typename Dest>
  static EIGEN_DONT_INLINE void run(const ProductType& prod, Dest& dest, typename ProductType::Scalar alpha) {
    typedef typename Dest::Index Index;
    // FIXME make sure rhs is sequentially stored
    // FIXME not very good if lhs is real and rhs complex while alpha is real too
    const Index rows = dest.rows();
    for (Index i=0; i<rows; ++i)
      dest.row(i) += (alpha * prod.lhs().coeff(i)) * prod.rhs();
  }
};

} // end namespace internal

/***********************************************************************
*  Implementation of General Matrix Vector Product
***********************************************************************/

/*  According to the shape/flags of the matrix we have to distinghish 3 different cases:
 *   1 - the matrix is col-major, BLAS compatible and M is large => call fast BLAS-like colmajor routine
 *   2 - the matrix is row-major, BLAS compatible and N is large => call fast BLAS-like rowmajor routine
 *   3 - all other cases are handled using a simple loop along the outer-storage direction.
 *  Therefore we need a lower level meta selector.
 *  Furthermore, if the matrix is the rhs, then the product has to be transposed.
 */
namespace internal {

template<typename Lhs, typename Rhs>
struct traits<GeneralProduct<Lhs,Rhs,GemvProduct> >
 : traits<ProductBase<GeneralProduct<Lhs,Rhs,GemvProduct>, Lhs, Rhs> >
{};

template<int Side, int StorageOrder, bool BlasCompatible>
struct gemv_selector;

} // end namespace internal

template<typename Lhs, typename Rhs>
class GeneralProduct<Lhs, Rhs, GemvProduct>
  : public ProductBase<GeneralProduct<Lhs,Rhs,GemvProduct>, Lhs, Rhs>
{
  public:
    EIGEN_PRODUCT_PUBLIC_INTERFACE(GeneralProduct)

    typedef typename Lhs::Scalar LhsScalar;
    typedef typename Rhs::Scalar RhsScalar;

    GeneralProduct(const Lhs& lhs, const Rhs& rhs) : Base(lhs,rhs)
    {
//       EIGEN_STATIC_ASSERT((internal::is_same<typename Lhs::Scalar, typename Rhs::Scalar>::value),
//         YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY)
    }

    enum { Side = Lhs::IsVectorAtCompileTime ? OnTheLeft : OnTheRight };
    typedef typename internal::conditional<int(Side)==OnTheRight,_LhsNested,_RhsNested>::type MatrixType;

    template<typename Dest> void scaleAndAddTo(Dest& dst, Scalar alpha) const
    {
      eigen_assert(m_lhs.rows() == dst.rows() && m_rhs.cols() == dst.cols());
      internal::gemv_selector<Side,(int(MatrixType::Flags)&RowMajorBit) ? RowMajor : ColMajor,
                       bool(internal::blas_traits<MatrixType>::HasUsableDirectAccess)>::run(*this, dst, alpha);
    }
};

namespace internal {

// The vector is on the left => transposition
template<int StorageOrder, bool BlasCompatible>
struct gemv_selector<OnTheLeft,StorageOrder,BlasCompatible>
{
  template<typename ProductType, typename Dest>
  static void run(const ProductType& prod, Dest& dest, typename ProductType::Scalar alpha)
  {
    Transpose<Dest> destT(dest);
    enum { OtherStorageOrder = StorageOrder == RowMajor ? ColMajor : RowMajor };
    gemv_selector<OnTheRight,OtherStorageOrder,BlasCompatible>
      ::run(GeneralProduct<Transpose<typename ProductType::_RhsNested>,Transpose<typename ProductType::_LhsNested>, GemvProduct>
        (prod.rhs().transpose(), prod.lhs().transpose()), destT, alpha);
  }
};

template<> struct gemv_selector<OnTheRight,ColMajor,true>
{
  template<typename ProductType, typename Dest>
  static void run(const ProductType& prod, Dest& dest, typename ProductType::Scalar alpha)
  {
    typedef typename ProductType::Index Index;
    typedef typename ProductType::LhsScalar   LhsScalar;
    typedef typename ProductType::RhsScalar   RhsScalar;
    typedef typename ProductType::Scalar      ResScalar;
    typedef typename ProductType::RealScalar  RealScalar;
    typedef typename ProductType::ActualLhsType ActualLhsType;
    typedef typename ProductType::ActualRhsType ActualRhsType;
    typedef typename ProductType::LhsBlasTraits LhsBlasTraits;
    typedef typename ProductType::RhsBlasTraits RhsBlasTraits;
    typedef Map<Matrix<ResScalar,Dynamic,1>, Aligned> MappedDest;

    ActualLhsType actualLhs = LhsBlasTraits::extract(prod.lhs());
    ActualRhsType actualRhs = RhsBlasTraits::extract(prod.rhs());

    ResScalar actualAlpha = alpha * LhsBlasTraits::extractScalarFactor(prod.lhs())
                                  * RhsBlasTraits::extractScalarFactor(prod.rhs());

    enum {
      // FIXME find a way to allow an inner stride on the result if packet_traits<Scalar>::size==1
      EvalToDestAtCompileTime = Dest::InnerStrideAtCompileTime==1,
      ComplexByReal = (NumTraits<LhsScalar>::IsComplex) && (!NumTraits<RhsScalar>::IsComplex)
    };

    bool alphaIsCompatible = (!ComplexByReal) || (imag(actualAlpha)==RealScalar(0));
    bool evalToDest = EvalToDestAtCompileTime && alphaIsCompatible;
    
    RhsScalar compatibleAlpha = get_factor<ResScalar,RhsScalar>::run(actualAlpha);

    ResScalar* actualDest;
    if (evalToDest)
    {
      actualDest = &dest.coeffRef(0);
    }
    else
    {
      actualDest = ei_aligned_stack_new(ResScalar,dest.size());
      if(!alphaIsCompatible)
      {
        MappedDest(actualDest, dest.size()).setZero();
        compatibleAlpha = RhsScalar(1);
      }
      else
        MappedDest(actualDest, dest.size()) = dest;
    }

    general_matrix_vector_product
      <Index,LhsScalar,ColMajor,LhsBlasTraits::NeedToConjugate,RhsScalar,RhsBlasTraits::NeedToConjugate>::run(
        actualLhs.rows(), actualLhs.cols(),
        &actualLhs.const_cast_derived().coeffRef(0,0), actualLhs.outerStride(),
        actualRhs.data(), actualRhs.innerStride(),
        actualDest, 1,
        compatibleAlpha);

    if (!evalToDest)
    {
      if(!alphaIsCompatible)
        dest += actualAlpha * MappedDest(actualDest, dest.size());
      else
        dest = MappedDest(actualDest, dest.size());
      ei_aligned_stack_delete(ResScalar, actualDest, dest.size());
    }
  }
};

template<> struct gemv_selector<OnTheRight,RowMajor,true>
{
  template<typename ProductType, typename Dest>
  static void run(const ProductType& prod, Dest& dest, typename ProductType::Scalar alpha)
  {
    typedef typename ProductType::LhsScalar LhsScalar;
    typedef typename ProductType::RhsScalar RhsScalar;
    typedef typename ProductType::Scalar    ResScalar;
    typedef typename ProductType::Index Index;
    typedef typename ProductType::ActualLhsType ActualLhsType;
    typedef typename ProductType::ActualRhsType ActualRhsType;
    typedef typename ProductType::_ActualRhsType _ActualRhsType;
    typedef typename ProductType::LhsBlasTraits LhsBlasTraits;
    typedef typename ProductType::RhsBlasTraits RhsBlasTraits;

    ActualLhsType actualLhs = LhsBlasTraits::extract(prod.lhs());
    ActualRhsType actualRhs = RhsBlasTraits::extract(prod.rhs());

    ResScalar actualAlpha = alpha * LhsBlasTraits::extractScalarFactor(prod.lhs())
                                  * RhsBlasTraits::extractScalarFactor(prod.rhs());

    enum {
      // FIXME I think here we really have to check for packet_traits<Scalar>::size==1
      // because in this case it is fine to have an inner stride
      DirectlyUseRhs = ((packet_traits<RhsScalar>::size==1) || (_ActualRhsType::Flags&ActualPacketAccessBit))
                     && (!(_ActualRhsType::Flags & RowMajorBit))
    };

    RhsScalar* rhs_data;
    if (DirectlyUseRhs)
       rhs_data = &actualRhs.const_cast_derived().coeffRef(0);
    else
    {
      rhs_data = ei_aligned_stack_new(RhsScalar, actualRhs.size());
      Map<typename _ActualRhsType::PlainObject>(rhs_data, actualRhs.size()) = actualRhs;
    }

    general_matrix_vector_product
      <Index,LhsScalar,RowMajor,LhsBlasTraits::NeedToConjugate,RhsScalar,RhsBlasTraits::NeedToConjugate>::run(
        actualLhs.rows(), actualLhs.cols(),
        &actualLhs.const_cast_derived().coeffRef(0,0), actualLhs.outerStride(),
        rhs_data, 1,
        &dest.coeffRef(0,0), dest.innerStride(),
        actualAlpha);

    if (!DirectlyUseRhs) ei_aligned_stack_delete(RhsScalar, rhs_data, prod.rhs().size());
  }
};

template<> struct gemv_selector<OnTheRight,ColMajor,false>
{
  template<typename ProductType, typename Dest>
  static void run(const ProductType& prod, Dest& dest, typename ProductType::Scalar alpha)
  {
    typedef typename Dest::Index Index;
    // TODO makes sure dest is sequentially stored in memory, otherwise use a temp
    const Index size = prod.rhs().rows();
    for(Index k=0; k<size; ++k)
      dest += (alpha*prod.rhs().coeff(k)) * prod.lhs().col(k);
  }
};

template<> struct gemv_selector<OnTheRight,RowMajor,false>
{
  template<typename ProductType, typename Dest>
  static void run(const ProductType& prod, Dest& dest, typename ProductType::Scalar alpha)
  {
    typedef typename Dest::Index Index;
    // TODO makes sure rhs is sequentially stored in memory, otherwise use a temp
    const Index rows = prod.rows();
    for(Index i=0; i<rows; ++i)
      dest.coeffRef(i) += alpha * (prod.lhs().row(i).cwiseProduct(prod.rhs().transpose())).sum();
  }
};

} // end namespace internal

/***************************************************************************
* Implementation of matrix base methods
***************************************************************************/

/** \returns the matrix product of \c *this and \a other.
  *
  * \note If instead of the matrix product you want the coefficient-wise product, see Cwise::operator*().
  *
  * \sa lazyProduct(), operator*=(const MatrixBase&), Cwise::operator*()
  */
template<typename Derived>
template<typename OtherDerived>
inline const typename ProductReturnType<Derived,OtherDerived>::Type
MatrixBase<Derived>::operator*(const MatrixBase<OtherDerived> &other) const
{
  // A note regarding the function declaration: In MSVC, this function will sometimes
  // not be inlined since DenseStorage is an unwindable object for dynamic
  // matrices and product types are holding a member to store the result.
  // Thus it does not help tagging this function with EIGEN_STRONG_INLINE.
  enum {
    ProductIsValid =  Derived::ColsAtCompileTime==Dynamic
                   || OtherDerived::RowsAtCompileTime==Dynamic
                   || int(Derived::ColsAtCompileTime)==int(OtherDerived::RowsAtCompileTime),
    AreVectors = Derived::IsVectorAtCompileTime && OtherDerived::IsVectorAtCompileTime,
    SameSizes = EIGEN_PREDICATE_SAME_MATRIX_SIZE(Derived,OtherDerived)
  };
  // note to the lost user:
  //    * for a dot product use: v1.dot(v2)
  //    * for a coeff-wise product use: v1.cwiseProduct(v2)
  EIGEN_STATIC_ASSERT(ProductIsValid || !(AreVectors && SameSizes),
    INVALID_VECTOR_VECTOR_PRODUCT__IF_YOU_WANTED_A_DOT_OR_COEFF_WISE_PRODUCT_YOU_MUST_USE_THE_EXPLICIT_FUNCTIONS)
  EIGEN_STATIC_ASSERT(ProductIsValid || !(SameSizes && !AreVectors),
    INVALID_MATRIX_PRODUCT__IF_YOU_WANTED_A_COEFF_WISE_PRODUCT_YOU_MUST_USE_THE_EXPLICIT_FUNCTION)
  EIGEN_STATIC_ASSERT(ProductIsValid || SameSizes, INVALID_MATRIX_PRODUCT)
#ifdef EIGEN_DEBUG_PRODUCT
  internal::product_type<Derived,OtherDerived>::debug();
#endif
  return typename ProductReturnType<Derived,OtherDerived>::Type(derived(), other.derived());
}

/** \returns an expression of the matrix product of \c *this and \a other without implicit evaluation.
  *
  * The returned product will behave like any other expressions: the coefficients of the product will be
  * computed once at a time as requested. This might be useful in some extremely rare cases when only
  * a small and no coherent fraction of the result's coefficients have to be computed.
  *
  * \warning This version of the matrix product can be much much slower. So use it only if you know
  * what you are doing and that you measured a true speed improvement.
  *
  * \sa operator*(const MatrixBase&)
  */
template<typename Derived>
template<typename OtherDerived>
const typename LazyProductReturnType<Derived,OtherDerived>::Type
MatrixBase<Derived>::lazyProduct(const MatrixBase<OtherDerived> &other) const
{
  enum {
    ProductIsValid =  Derived::ColsAtCompileTime==Dynamic
                   || OtherDerived::RowsAtCompileTime==Dynamic
                   || int(Derived::ColsAtCompileTime)==int(OtherDerived::RowsAtCompileTime),
    AreVectors = Derived::IsVectorAtCompileTime && OtherDerived::IsVectorAtCompileTime,
    SameSizes = EIGEN_PREDICATE_SAME_MATRIX_SIZE(Derived,OtherDerived)
  };
  // note to the lost user:
  //    * for a dot product use: v1.dot(v2)
  //    * for a coeff-wise product use: v1.cwiseProduct(v2)
  EIGEN_STATIC_ASSERT(ProductIsValid || !(AreVectors && SameSizes),
    INVALID_VECTOR_VECTOR_PRODUCT__IF_YOU_WANTED_A_DOT_OR_COEFF_WISE_PRODUCT_YOU_MUST_USE_THE_EXPLICIT_FUNCTIONS)
  EIGEN_STATIC_ASSERT(ProductIsValid || !(SameSizes && !AreVectors),
    INVALID_MATRIX_PRODUCT__IF_YOU_WANTED_A_COEFF_WISE_PRODUCT_YOU_MUST_USE_THE_EXPLICIT_FUNCTION)
  EIGEN_STATIC_ASSERT(ProductIsValid || SameSizes, INVALID_MATRIX_PRODUCT)

  return typename LazyProductReturnType<Derived,OtherDerived>::Type(derived(), other.derived());
}

#endif // EIGEN_PRODUCT_H