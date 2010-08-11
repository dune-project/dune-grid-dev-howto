// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDHOSTGRIDACCESS_HH
#define DUNE_IDENTITYGRIDHOSTGRIDACCESS_HH

#include <string>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class >
  class IdentityGrid;



  // HostGridAccess
  // --------------

  template< class Grid >
  struct HostGridAccess;



  /** \class HostGridAccess
   *  \brief provides access to host grid objects
   *
   *  \tparam  Grid  meta grid, whose host grid shall be accessed
   *
   *  \nosubgrouping
   */
  template< class HG >
  struct HostGridAccess< IdentityGrid< HG > >
  {
    /** \name Exported Types
     * \{ */

    typedef IdentityGrid< HG > Grid;

    //! type of HostGrid
    typedef HG HostGrid;

    /** \} */

    /** \brief A Traits struct that collects return types of class member methods.
     *
     *  \tparam codim codimension
     */
    template< int codim >
    struct Codim
    {
      //! type of the IdGrid entity
      typedef typename Grid::template Codim< codim >::Entity Entity;
      //! type of the IdGrid entity pointer
      typedef typename Grid::template Codim< codim >::EntityPointer EntityPointer;

      //! type of the host entity
      typedef typename HostGrid::template Codim< codim >::Entity HostEntity;
      //! type of the host entity pointer
      typedef typename HostGrid::template Codim< codim >::EntityPointer HostEntityPointer;
    };

    //! type of the IdGrid leaf intersection
    typedef typename Grid::Traits::LeafIntersection LeafIntersection;
    //! type of the IdGrid level intersection
    typedef typename Grid::Traits::LevelIntersection LevelIntersection;

    //! type of the host leaf intersection
    typedef typename HostGrid::Traits::LeafIntersection HostLeafIntersection;
    //! type of the host level intersection
    typedef typename HostGrid::Traits::LevelIntersection HostLevelIntersection;

    /** \brief Get underlying HostGrid.
     *  \param[in]  grid  grid, whose host grid shall be returned
     *  \returns HostGrid
     */
    static const HostGrid &hostGrid ( const Grid &grid )
    {
      return grid.getHostGrid();
    }

    template< class Entity >
    static const typename Codim< Entity::codimension >::HostEntity &
    hostEntity ( const Entity &entity )
    {
      return hostEntity< Entity::codimension >( entity );
    }

    template< int codim >
    static const typename Codim< codim >::HostEntity &
    hostEntity ( const typename Codim< codim >::Entity &entity )
    {
      return *Grid::getRealImplementation( entity ).hostEntity_;
    }

    static const HostLeafIntersection &
    hostIntersection ( const LeafIntersection &intersection )
    {
      return *Grid::getRealImplementation( intersection ).hostIterator_;
    }

    static const HostLevelIntersection &
    hostIntersection ( const LevelIntersection &intersection )
    {
      return *Grid::getRealImplementation( intersection ).hostIterator_;
    }
  };

}

#endif // #ifndef DUNE_IDENTITYGRIDHOSTGRIDACCESS_HH
