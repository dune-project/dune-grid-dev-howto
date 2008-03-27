// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRID_INDEXSETS_HH
#define DUNE_IDENTITYGRID_INDEXSETS_HH

/** \file
 * \brief The index and id sets for the IdentityGrid class
 */

#include <vector>

namespace Dune {


  template <class GridImp>
  struct IdentityGridLevelIndexSetTypes
  {
    //! The types
    template<int cd>
    struct Codim
    {
      template<PartitionIteratorType pitype>
      struct Partition
      {
        typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LevelIterator Iterator;
      };
    };
  };




  /** \todo Take the index types from the host grid */
  template<class GridImp>
  class IdentityGridLevelIndexSet :
    public IndexSetDefaultImplementation<GridImp,IdentityGridLevelIndexSet<GridImp>,IdentityGridLevelIndexSetTypes<GridImp> >
  {
  public:

    typedef typename remove_const<GridImp>::type::HostGridType HostGrid;

    enum {dim = GridImp::dimension};

    typedef IndexSet<GridImp,IdentityGridLevelIndexSet<GridImp>,IdentityGridLevelIndexSetTypes<GridImp> > Base;

    //! get index of an entity
    template<int codim>
    int index (const typename GridImp::Traits::template Codim<codim>::Entity& e) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).template index<codim>(*grid_->template getHostEntity<codim>(e));
    }


    //! get index of subEntity of a codim 0 entity
    template<int codim>
    int subIndex (const typename GridImp::Traits::template Codim<0>::Entity& e, int i) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).template subIndex<codim>(*grid_->template getHostEntity<0>(e), i);
    }


    //! get number of entities of given codim, type and on this level
    int size (int codim) const {
      return grid_->hostgrid_->levelIndexSet(level_).size(codim);
    }


    //! get number of entities of given codim, type and on this level
    int size (GeometryType type) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).size(type);
    }


    /** \brief Deliver all geometry types used in this grid */
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return grid_->hostgrid_->levelIndexSet(level_).geomTypes(codim);
    }


    //! one past the end on this level
    template<int cd, PartitionIteratorType pitype>
    typename Base::template Codim<cd>::template Partition<pitype>::Iterator begin () const
    {
      return grid_->template lbegin<cd,pitype>(level_);
    }


    //! Iterator to one past the last entity of given codim on level for partition type
    template<int cd, PartitionIteratorType pitype>
    typename Base::template Codim<cd>::template Partition<pitype>::Iterator end () const
    {
      return grid_->template lend<cd,pitype>(level_);
    }


    /** \brief Set up the index set */
    void update(const GridImp& grid, int level)
    {
      grid_ = &grid;
      level_ = level;
    }


    GridImp* grid_;

    int level_;
  };




  template <class GridImp>
  struct IdentityGridLeafIndexSetTypes
  {
    //! The types
    template<int cd>
    struct Codim
    {
      template<PartitionIteratorType pitype>
      struct Partition
      {
        typedef typename GridImp::template Codim<cd>::template Partition<pitype>::LeafIterator Iterator;
      };
    };
  };




  template<class GridImp>
  class IdentityGridLeafIndexSet :
    public IndexSetDefaultImplementation<GridImp,IdentityGridLeafIndexSet<GridImp>,IdentityGridLeafIndexSetTypes<GridImp> >
  {
    typedef typename remove_const<GridImp>::type::HostGridType HostGrid;

  public:


    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instatiated yet.
     */
    enum {dim = remove_const<GridImp>::type::dimension};

    typedef IndexSet<GridImp,IdentityGridLeafIndexSet<GridImp>,IdentityGridLeafIndexSetTypes<GridImp> > Base;


    //! constructor stores reference to a grid and level
    IdentityGridLeafIndexSet (const GridImp& grid)
      : grid_(&grid)
    {}


    //! get index of an entity
    /*
        We use the RemoveConst to extract the Type from the mutable class,
        because the const class is not instatiated yet.
     */
    template<int codim>
    int index (const typename remove_const<GridImp>::type::template Codim<codim>::Entity& e) const
    {
      return grid_->hostgrid_->leafIndexSet().template index<codim>(*grid_->template getHostEntity<codim>(e));
    }


    //! get index of subEntity of a codim 0 entity
    /*
        We use the RemoveConst to extract the Type from the mutable class,
        because the const class is not instatiated yet.
     */
    template<int codim>
    int subIndex (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i) const
    {
      return grid_->hostgrid_->leafIndexSet().template subIndex<codim>(*grid_->template getHostEntity<0>(e),i);
    }


    //! get number of entities of given type
    int size (GeometryType type) const
    {
      return grid_->hostgrid_->leafIndexSet().size(type);
    }


    //! get number of entities of given codim
    int size (int codim) const
    {
      return grid_->hostgrid_->leafIndexSet().size(codim);
    }


    /** \brief Deliver all geometry types used in this grid */
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return grid_->hostgrid_->leafIndexSet().geomTypes(codim);
    }


    //! one past the end on this level
    template<int codim, PartitionIteratorType pitype>
    typename Base::template Codim<codim>::template Partition<pitype>::Iterator begin () const
    {
      return grid_->template leafbegin<codim,pitype>();
    }


    //! Iterator to one past the last entity of given codim on level for partition type
    template<int codim, PartitionIteratorType pitype>
    typename Base::template Codim<codim>::template Partition<pitype>::Iterator end () const
    {
      return grid_->template leafend<codim,pitype>();
    }


    /** \todo Currently we support only vertex and element indices */
    void update(const GridImp& grid)
    {
      grid_ = &grid;
    }


    GridImp* grid_;
  };




  template <class GridImp>
  class IdentityGridGlobalIdSet :
    public IdSet<GridImp,IdentityGridGlobalIdSet<GridImp>,
        typename remove_const<GridImp>::type::HostGridType::Traits::GlobalIdSet::IdType>
  {

    typedef typename remove_const<GridImp>::type::HostGridType HostGrid;


  public:
    //! constructor stores reference to a grid
    IdentityGridGlobalIdSet (const GridImp& g) : grid_(&g) {}

    //! define the type used for persistent indices
    typedef typename HostGrid::Traits::GlobalIdSet::IdType GlobalIdType;


    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instatiated yet.
     */
    template<int cd>
    GlobalIdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      // Return id of the host entity
      return grid_->hostgrid_->globalIdSet().id(*grid_->getRealImplementation(e).hostEntity_);
    }


    //! get id of subEntity
    /*
        We use the remove_const to extract the Type from the mutable class,
        because the const class is not instatiated yet.
     */
    template<int cc>
    GlobalIdType subId (const typename remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i) const
    {
      // Return sub id of the host entity
      return grid_->hostgrid_->globalIdSet().template subId<cc>(*grid_->getRealImplementation(e).hostEntity_,i);
    }


    /** \todo Should be private */
    void update() {}


    const GridImp* grid_;
  };




  template<class GridImp>
  class IdentityGridLocalIdSet :
    public IdSet<GridImp,IdentityGridLocalIdSet<GridImp>,
        typename remove_const<GridImp>::type::HostGridType::Traits::LocalIdSet::IdType>
  {
  private:

    typedef typename remove_const<GridImp>::type::HostGridType HostGrid;


  public:
    //! define the type used for persistent local ids
    typedef typename HostGrid::Traits::LocalIdSet::IdType LocalIdType;


    //! constructor stores reference to a grid
    IdentityGridLocalIdSet (const GridImp& g) : grid_(&g) {}


    //! get id of an entity
    /*
        We use the remove_const to extract the Type from the mutable class,
        because the const class is not instatiated yet.
     */
    template<int cd>
    LocalIdType id (const typename remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      // Return id of the host entity
      return grid_->hostgrid_->localIdSet().id(*grid_->getRealImplementation(e).hostEntity_);
    }


    //! get id of subEntity
    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instatiated yet.
     */
    template<int cc>
    LocalIdType subId (const typename remove_const<GridImp>::type::template Codim<0>::Entity& e, int i) const
    {
      // Return sub id of the host entity
      return grid_->hostgrid_->localIdSet().template subId<cc>(*grid_->getRealImplementation(e).hostEntity_,i);
    }


    /** \todo Should be private */
    void update() {}


    const GridImp* grid_;
  };


}  // namespace Dune


#endif
