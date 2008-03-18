// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDGEOMETRY_HH
#define DUNE_IDENTITYGRIDGEOMETRY_HH

/** \file
 * \brief The IdentityGridGeometry class and its specializations
 */

#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>

namespace Dune {



  // forward declaration
  template<int mydim, int coorddim, class GridImp>
  class IdentityGridGeometry;




  template<int mydim, int coorddim, class GridImp>
  class IdentityGridMakeableGeometry :
    public Geometry<mydim, coorddim, GridImp, IdentityGridGeometry>
  {
  public:

    // The codimension of this entitypointer wrt the host grid
    enum {CodimInHostGrid = GridImp::HostGridType::dimension - mydim};

    enum {DimensionWorld = GridImp::HostGridType::dimensionworld};

    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridGeometryType;

    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridLocalGeometryType;

    // select appropiate hostgrid geometry via typeswitch
    typedef typename SelectType<coorddim==DimensionWorld, HostGridGeometryType, HostGridLocalGeometryType>::Type HostGridGeometry;


    //! \todo Please doc me !
    IdentityGridMakeableGeometry(const HostGridGeometry& hostGeometry)
      : Geometry<mydim, coorddim, GridImp, IdentityGridGeometry>(IdentityGridGeometry<mydim, coorddim, GridImp>(hostGeometry))
    {};

  };




  template<int mydim, int coorddim, class GridImp>
  class IdentityGridGeometry :
    public GeometryDefaultImplementation <mydim, coorddim, GridImp, IdentityGridGeometry>
  {
  private:

    typedef typename GridImp::ctype ctype;


  public:

    // The codimension of this entitypointer wrt the host grid
    enum {CodimInHostGrid = GridImp::HostGridType::dimension - mydim};
    enum {DimensionWorld = GridImp::HostGridType::dimensionworld};

    // select appropiate hostgrid geometry via typeswitch
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridGeometryType;
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::Geometry HostGridLocalGeometryType;

    typedef typename SelectType<coorddim==DimensionWorld, HostGridGeometryType, HostGridLocalGeometryType>::Type HostGridGeometry;


    /** Default constructor.
     */
    IdentityGridGeometry(const HostGridGeometry& hostGeometry)
      : hostGeometry_(hostGeometry)
    {}


    /** \brief Return the element type identifier
     */
    GeometryType type () const {
      return hostGeometry_.type();
    }


    //! return the number of corners of this element. Corners are numbered 0...n-1
    int corners () const {
      return hostGeometry_.corners();
    }


    //! access to coordinates of corners. Index is the number of the corner
    const FieldVector<ctype, coorddim>& operator[] (int i) const {
      return hostGeometry_[i];
    }


    /** \brief Maps a local coordinate within reference element to
     * global coordinate in element  */
    FieldVector<ctype, coorddim> global (const FieldVector<ctype, mydim>& local) const {
      return hostGeometry_.global(local);
    }


    /** \brief Maps a global coordinate within the element to a
     * local coordinate in its reference element */
    FieldVector<ctype, mydim> local (const FieldVector<ctype, coorddim>& global) const {
      return hostGeometry_.local(global);
    }


    //! Returns true if the point is in the current element
    bool checkInside(const FieldVector<ctype, mydim> &local) const {
      return hostGeometry_.checkInside(local);
    }


    /**
     */
    ctype integrationElement (const FieldVector<ctype, mydim>& local) const {
      return hostGeometry_.integrationElement(local);
    }


    //! The Jacobian matrix of the mapping from the reference element to this element
    const FieldMatrix<ctype, mydim,mydim>& jacobianInverseTransposed (const FieldVector<ctype, mydim>& local) const {
      return hostGeometry_.jacobianInverseTransposed(local);
    }


    const HostGridGeometry& hostGeometry_;

  };


}  // namespace Dune

#endif
