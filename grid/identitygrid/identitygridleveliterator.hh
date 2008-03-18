// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDLEVELITERATOR_HH
#define DUNE_IDENTITYGRIDLEVELITERATOR_HH

/** \file
 * \brief The IdentityGridLevelIterator class and its specializations
 */

namespace Dune {




  //**********************************************************************
  //
  // --IdentityGridLevelIterator
  /** \brief Iterator over all entities of a given codimension and level of a grid.
   * \ingroup IdentityGrid
   */
  template<int codim, PartitionIteratorType pitype, class GridImp>
  class IdentityGridLevelIterator :
    public Dune::IdentityGridEntityPointer <codim,GridImp>,
    public LevelIteratorDefaultImplementation <codim,pitype,GridImp,IdentityGridLevelIterator>
  {
  private:

    enum {dim = GridImp::dimension};


  public:

    //! Constructor
    explicit IdentityGridLevelIterator(const GridImp* subGrid, int level)
      : IdentityGridEntityPointer<codim,GridImp>(subGrid, subGrid->hostgrid_->template lbegin<codim>(level)),
        hostGridLevelIterator_(subGrid->hostgrid_->template lbegin<codim>(level)),
        hostGridLevelEndIterator_(subGrid->hostgrid_->template lend<codim>(level))
    {
      this->virtualEntity_.setToTarget(hostGridLevelIterator_);
    }


    /** \brief Constructor which create the end iterator
        \param endDummy Here only to distinguish it from the other constructor
     */
    explicit IdentityGridLevelIterator(const GridImp* subGrid, int level, bool endDummy)
      :
        IdentityGridEntityPointer<codim,GridImp>(subGrid, subGrid->hostgrid_->template lend<codim>(level)),
        hostGridLevelIterator_(subGrid->hostgrid_->template lend<codim>(level)),
        hostGridLevelEndIterator_(subGrid->hostgrid_->template lend<codim>(level))
    {}


    //! prefix increment
    void increment() {
      ++hostGridLevelIterator_;
      this->virtualEntity_.setToTarget(hostGridLevelIterator_);
    }


  private:

    // LevelIterator to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::Traits::template Codim<codim>::LevelIterator HostGridLevelIterator;

    //! \todo Please doc me !
    HostGridLevelIterator hostGridLevelIterator_;

    //! \todo Please doc me !
    HostGridLevelIterator hostGridLevelEndIterator_;

  };


}  // namespace Dune

#endif
