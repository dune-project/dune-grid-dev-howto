// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRID_HH
#define DUNE_IDENTITYGRID_HH

/** \file
 * \brief The IdentityGrid class
 */

#include <string>
#include <map>

#include <dune/common/collectivecommunication.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/common/bitfield.hh>
#include <dune/grid/common/grid.hh>
#include <dune/common/timer.hh>

// The components of the Identitygrid interface
#include "identitygrid/identitygridgeometry.hh"
#include "identitygrid/identitygridentity.hh"
#include "identitygrid/identitygridentitypointer.hh"
#include "identitygrid/identitygridintersectioniterator.hh"
#include "identitygrid/identitygridleveliterator.hh"
#include "identitygrid/identitygridleafiterator.hh"
#include "identitygrid/identitygridhierarchiciterator.hh"
#include "identitygrid/identitygridindexsets.hh"

namespace Dune {

  // Forward declaration
  template <class HostGrid>
  class IdentityGrid;




  template<int dim, class HostGrid>
  struct IdentityGridFamily
  {
    typedef GridTraits<
        dim,
        HostGrid::dimensionworld,
        Dune::IdentityGrid<HostGrid>,
        IdentityGridGeometry,
        IdentityGridEntity,
        IdentityGridEntityPointer,
        IdentityGridLevelIterator,
        IdentityGridLeafIntersectionIterator,
        IdentityGridLevelIntersectionIterator,
        IdentityGridHierarchicIterator,
        IdentityGridLeafIterator,
        IdentityGridLevelIndexSet< const IdentityGrid<HostGrid> >,
        IdentityGridLevelIndexSetTypes< const IdentityGrid<HostGrid> >,
        IdentityGridLeafIndexSet< const IdentityGrid<HostGrid> >,
        IdentityGridLeafIndexSetTypes< const IdentityGrid<HostGrid> >,
        IdentityGridGlobalIdSet< const IdentityGrid<HostGrid> >,
        typename HostGrid::Traits::GlobalIdSet::IdType,
        IdentityGridLocalIdSet< const IdentityGrid<HostGrid> >,
        typename HostGrid::Traits::LocalIdSet::IdType,
        CollectiveCommunication<IdentityGrid<HostGrid> >
        > Traits;
  };




  //**********************************************************************
  //
  // --IdentityGrid
  //
  //**********************************************************************

  /** \brief [<em> provides \ref Dune::Grid </em>]
   *
   */
  template <class HostGrid>
  class IdentityGrid :
    public GridDefaultImplementation  <HostGrid::dimension, HostGrid::dimensionworld, double, IdentityGridFamily<HostGrid::dimension,HostGrid> >
  {

    friend class IdentityGridLevelIndexSet<const IdentityGrid<HostGrid> >;
    friend class IdentityGridLeafIndexSet<const IdentityGrid<HostGrid> >;
    friend class IdentityGridGlobalIdSet<const IdentityGrid<HostGrid> >;
    friend class IdentityGridLocalIdSet<const IdentityGrid<HostGrid> >;
    friend class IdentityGridHierarchicIterator<const IdentityGrid<HostGrid> >;
    friend class IdentityGridLevelIntersectionIterator<const IdentityGrid<HostGrid> >;
    friend class IdentityGridLeafIntersectionIterator<const IdentityGrid<HostGrid> >;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class IdentityGridLevelIterator;

    template<int codim, PartitionIteratorType pitype, class GridImp_>
    friend class IdentityGridLeafIterator;


    template<int codim_, int dim_, class GridImp_>
    friend class IdentityGridEntity;

  public:

    /** \todo Should not be public */
    typedef HostGrid HostGridType;

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    //! type of the used GridFamily for this grid
    typedef IdentityGridFamily<HostGrid::dimension,HostGrid>  GridFamily;

    //! the Traits
    typedef typename IdentityGridFamily<HostGrid::dimension,HostGrid>::Traits Traits;

    //! The type used to store coordinates, inherited from the HostGrid
    typedef typename HostGrid::ctype ctype;


    /** \brief Constructor
     */
    explicit IdentityGrid(HostGrid& hostgrid) :
      hostgrid_(&hostgrid),
      leafIndexSet_(*this),
      globalIdSet_(*this),
      localIdSet_(*this)
    {
      setIndices();
    }


    //! Desctructor
    ~IdentityGrid()
    {
      // Delete level index sets
      for (size_t i=0; i<levelIndexSets_.size(); i++)
        if (levelIndexSets_[i])
          delete (levelIndexSets_[i]);
    }


    //! return grid name
    std::string name() const
    {
      return "IdentityGrid";
    }


    //! Return maximum level defined in this grid. Levels are numbered
    //! 0 ... maxlevel with 0 the coarsest level.
    int maxLevel() const {
      return hostgrid_->maxLevel();
    }


    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
      return IdentityGridLevelIterator<codim,All_Partition, const IdentityGrid<HostGrid> >(this, level);
    }


    //! one past the end on this level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
      return IdentityGridLevelIterator<codim,All_Partition, const IdentityGrid<HostGrid> >(this, level, true);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
      return IdentityGridLevelIterator<codim,PiType, const IdentityGrid<HostGrid> >(this, level);
    }


    //! one past the end on this level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
      return IdentityGridLevelIterator<codim,PiType, const IdentityGrid<HostGrid> >(this, level, true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
      return IdentityGridLeafIterator<codim,All_Partition, const IdentityGrid<HostGrid> >(this);
    }


    //! one past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend() const {
      return IdentityGridLeafIterator<codim,All_Partition, const IdentityGrid<HostGrid> >(this, true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
      return IdentityGridLeafIterator<codim,PiType, const IdentityGrid<HostGrid> >(this);
    }


    //! one past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
      return IdentityGridLeafIterator<codim,PiType, const IdentityGrid<HostGrid> >(this, true);
    }


    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const {
      return hostgrid_->size(level,codim);
    }


    //! number of leaf entities per codim in this process
    int size (int codim) const {
      return leafIndexSet().size(codim);
    }


    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const {
      return levelIndexSets_[level]->size(type);
    }


    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const
    {
      return leafIndexSet().size(type);
    }


    /** \brief Access to the GlobalIdSet */
    const typename Traits::GlobalIdSet& globalIdSet() const {
      return globalIdSet_;
    }


    /** \brief Access to the LocalIdSet */
    const typename Traits::LocalIdSet& localIdSet() const {
      return localIdSet_;
    }


    /** \brief Access to the LevelIndexSets */
    const typename Traits::LevelIndexSet& levelIndexSet(int level) const
    {
      if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
      return *levelIndexSets_[level];
    }


    /** \brief Access to the LeafIndexSet */
    const typename Traits::LeafIndexSet& leafIndexSet() const
    {
      return leafIndexSet_;
    }


    /** @name Grid Refinement Methods */
    /*@{*/


#warning Adaptation stuff is not compiled!
#if 0
    /** global refinement
     * \todo optimize implementation
     */
    void globalRefine (int refCount)
    {
      typedef typename Traits::template Codim<0>::LeafIterator ElementLeafIterator;
      ElementLeafIterator it = leafbegin<0>();
      ElementLeafIterator end = leafend<0>();
      for(; it!=end; ++it)
        mark(1, it);
      preAdapt();
      adapt();
    }

    /** \brief Mark entity for refinement
     *
     * This only works for entities of codim 0.
     * The parameter is currently ignored
     *
     * \return <ul>
     * <li> true, if marking was succesfull </li>
     * <li> false, if marking was not possible </li>
     * </ul>
     */
    bool mark(int refCount, const typename Traits::template Codim<0>::EntityPointer & e)
    {
      return hostgrid_->mark(refCount, getHostEntityPointer(e));
    }


    /** \brief Return refinement mark for entity
     *
     * \return refinement mark (1,0,-1)
     */
    int getMark(const typename Traits::template Codim<0>::EntityPointer & e) const
    {
      if ((adaptationStep_==preAdaptDone)or (adaptationStep_==adaptDone))
        DUNE_THROW(InvalidStateException, "You can not use getMark() after preAdapt() or adapt() !");

      if (not (e->isLeaf()))
        return 0;

      int level = e->level();
      int index = levelIndexSet(level).index(*e);

      if (refinementMark_[level][index])
        return 1;
      if (coarseningMark_[level][index])
        return -1;

      return 0;
    }


    //! \todo Please doc me !
    bool preAdapt(){
      if (adaptationStep_==preAdaptDone)
      {
        std::cout << "You already called call preAdapt() ! Aborting preAdapt()." << std::endl;
        return false;
      }
      if (adaptationStep_==adaptDone)
      {
        std::cout << "You did not call postAdapt() after adapt() ! Calling postAdapt() automatically." << std::endl;
        postAdapt();
      }
      adaptationStep_ = nothingDone;
      hostAdaptationStep_ = nothingDone;

            #ifdef SUBGRID_VERBOSE
      std::cout << "preadapt 1 (start)" << std::endl;
            #endif

      typedef typename Traits::template Codim<0>::LevelIterator ElementLevelIterator;
      typedef typename Traits::template Codim<0>::LeafIterator ElementLeafIterator;
      typedef typename Traits::template Codim<0>::Entity::HierarchicIterator HierarchicIterator;
      typedef typename Traits::LevelIndexSet LevelIndexSet;
      typedef typename Traits::template Codim<0>::Entity::EntityPointer ElementPointer;

      std::vector<int> maxElementLevel(size(dim), 0);

      ElementLeafIterator it = leafbegin<0>();
      ElementLeafIterator end = leafend<0>();
      for (; it!=end; ++it)
      {
        int index = levelIndexSet(it->level()).index(*it);
        ElementPointer element = it;
        int newLevel = it->level();
        if (refinementMark_[it->level()][index])
          ++newLevel;
        if (coarseningMark_[it->level()][index])
          --newLevel;

        for (int i = 0; i<it->template count<dim>(); ++i)
        {
          int nodeIndex = leafIndexSet().template subIndex<dim>(*it, i);
          if (maxElementLevel[nodeIndex] < newLevel)
            maxElementLevel[nodeIndex] = newLevel;
        }
      }

      // do the following twice to recognize all cases including coarsening
      for (int i=0; i<2; ++i)
      {
        for (int level=maxLevel(); level>=0; --level)
        {
          ElementLevelIterator it = lbegin<0>(level);
          ElementLevelIterator end = lend<0>(level);
          for (; it!=end; ++it)
          {
            if(it->isLeaf())
            {
              int index = levelIndexSet(level).index(*it);

              int maxNeighborLevel = 0;
              for (int i = 0; i<it->template count<dim>(); ++i)
              {
                int nodeIndex = leafIndexSet().template subIndex<dim>(*it, i);
                if (maxNeighborLevel < maxElementLevel[nodeIndex])
                  maxNeighborLevel = maxElementLevel[nodeIndex];
              }
              int minLevel = maxNeighborLevel-maxLevelDiff_;
              for (int i = 0; i<it->template count<dim>(); ++i)
              {
                int nodeIndex = leafIndexSet().template subIndex<dim>(*it, i);
                if (maxElementLevel[nodeIndex] < minLevel)
                  maxElementLevel[nodeIndex] = minLevel;
              }

              if ((coarseningMark_[level][index])and (level-1<minLevel))
                coarseningMark_[level][index] = false;

              if (not (coarseningMark_[level][index])and (level<minLevel))
                refinementMark_[level][index] = true;
            }
          }
        }
      }

            #ifdef SUBGRID_VERBOSE
      std::cout << "preadapt 2 (refinement/coarsening marks for neighbours set)" << std::endl;
            #endif

      bool mightBeCoarsened = false;
      bool callHostgridAdapt = false;

      it = leafbegin<0>();
      end = leafend<0>();
      for (; it!=end; ++it)
      {
        int level = it->level();
        int index = levelIndexSet(level).index(*it);

        // if element is marked for coarsening
        if (coarseningMark_[level][index])
        {
          int fatherLevel = it->father()->level();
          int fatherIndex = levelIndexSet(fatherLevel).index(*(it->father()));

          // if father is not processed
          //    then look for brothers
          if (not (refinementMark_[fatherLevel][fatherIndex])and not (coarseningMark_[fatherLevel][fatherIndex]))
          {
            // iterate over all brothers
            HierarchicIterator hIt = it->father()->hbegin(level);
            HierarchicIterator hEnd = it->father()->hend(level);
            for (; hIt!=hEnd; ++hIt)
            {
              // if brother is not marked for coarsening
              //    then mark father for not coarsening and stop
              int brotherIndex = levelIndexSet(level).index(*hIt);
              if (not (coarseningMark_[level][brotherIndex]))
              {
                refinementMark_[fatherLevel][fatherIndex] = true;
                break;
              }
            }
            // if father was not marked for not coarsening
            //    then mark it for coarsening
            if (not (refinementMark_[fatherLevel][fatherIndex]))
              coarseningMark_[fatherLevel][fatherIndex] = true;
          }

          // if father is marked for not coarsening
          //    then unset element mark for coarsening
          if (refinementMark_[fatherLevel][fatherIndex])
            coarseningMark_[level][index] = false;
        }
        else
        {
          if (refinementMark_[level][index])
          {
            if (getRealImplementation(*it).hostEntity_->isLeaf())
            {
              hostgrid_->mark(1, getRealImplementation(*it).hostEntity_);
              callHostgridAdapt = true;
            }
          }
        }
      }

            #ifdef SUBGRID_VERBOSE
      std::cout << "preadapt 3 (checked coarsening marks, marked host grid)" << std::endl;
            #endif

      if (callHostgridAdapt)
      {
        hostgrid_->preAdapt();

        hostAdaptationStep_ = preAdaptDone;

                #ifdef SUBGRID_VERBOSE
        std::cout << "preadapt 4 (host grid preadapt called)" << std::endl;
                #endif
      }

      adaptationStep_ = preAdaptDone;

      return mightBeCoarsened;
    }


    //! Triggers the grid refinement process
    bool adapt()
    {
      if ((adaptationStep_==nothingDone)or (adaptationStep_==postAdaptDone))
      {
        std::cout << "You did not call preAdapt() before adapt() ! Calling preAdapt() automatically." << std::endl;
        preAdapt();
      }
      if (adaptationStep_==adaptDone)
      {
        std::cout << "You already called call adapt() ! Aborting adapt()." << std::endl;
        return false;
      }

      typedef typename Traits::template Codim<0>::LevelIterator ElementLevelIterator;
      typedef typename HostGrid::Traits::template Codim<0>::LevelIterator HostElementLevelIterator;
      typedef typename HostGrid::template Codim<0>::Entity::HierarchicIterator HostHierarchicIterator;

      typedef typename Traits::GlobalIdSet GlobalIdSet;
      typedef typename GlobalIdSet::IdType GlobalIdType;

            #ifdef SUBGRID_VERBOSE
      std::cout << "adapt 1 (start)" << std::endl;
            #endif

      std::map<GlobalIdType, bool> elements;
      for(int level=0; level <= maxLevel(); ++level)
      {
        ElementLevelIterator it = lbegin<0>(level);
        ElementLevelIterator end = lend<0>(level);
        for (; it!=end; ++it)
        {
          int level = it->level();
          int index = levelIndexSet(level).index(*it);

          // if element is leaf
          //    then if not marked for coarsening store it and its refinement mark
          // if element is not leaf
          //    then store it and ignore refinement mark
          if (it->isLeaf())
          {

            if (not (coarseningMark_[level][index]))
              elements[globalIdSet().id(*it)] = refinementMark_[level][index];
          }
          else
            elements[globalIdSet().id(*it)] = false;
        }
      }

            #ifdef SUBGRID_VERBOSE
      std::cout << "adapt 2 (element ids stored)" << std::endl;
            #endif

      if (hostAdaptationStep_==preAdaptDone)
      {
        hostgrid_->adapt();

        hostAdaptationStep_ = adaptDone;

                #ifdef SUBGRID_VERBOSE
        std::cout << "adapt 3 (host grid adapt called)" << std::endl;
                #endif
      }

      createBegin();
      // recreate entity mark vectors

      // mark entities
      for(int level=0; level <= hostgrid_->maxLevel(); ++level)
      {
        HostElementLevelIterator it = hostgrid_->lbegin<0>(level);
        HostElementLevelIterator end = hostgrid_->lend<0>(level);
        for (; it!=end; ++it)
        {
          typename std::map<GlobalIdType, bool>::iterator contained = elements.find(hostgrid_->globalIdSet().id(*it));
          if (contained != elements.end())
          {
            add(it);
            if (contained->second)
            {
              HostHierarchicIterator hIt = it->hbegin(it->level()+1);
              HostHierarchicIterator hEnd = it->hend(it->level()+1);
              for (; hIt!=hEnd; ++hIt)
                add(hIt);
            }
          }
        }
      }

            #ifdef SUBGRID_VERBOSE
      std::cout << "adapt 4 (subgrid elements added)" << std::endl;
            #endif

      createEnd();

            #ifdef SUBGRID_VERBOSE
      std::cout << "adapt 5 (subgrid created)" << std::endl;
            #endif

      for(int level=0; level < maxLevel(); ++level)
      {
        ElementLevelIterator it = lbegin<0>(level);
        ElementLevelIterator end = lend<0>(level);
        for (; it!=end; ++it)
        {
          typename std::map<GlobalIdType, bool>::const_iterator e = elements.find(globalIdSet().id(*it));
          if (e != elements.end())
          {
            if (e->second == true)
              refinementMark_[level][levelIndexSet(level).index(*it)] = true;
          }
        }
      }

            #ifdef SUBGRID_VERBOSE
      std::cout << "adapt 6 (wasRefined marks set)" << std::endl;
            #endif

      adaptationStep_ = adaptDone;

      return true;
    }


    /** \brief Clean up refinement markers */
    void postAdapt(){
      if (adaptationStep_==nothingDone)
      {
        std::cout << "You did not call preAdapt() before adapt() ! Calling preAdapt() automatically." << std::endl;
        preAdapt();
      }
      if (adaptationStep_==preAdaptDone)
      {
        std::cout << "You did not call adapt() before postAdaptt() ! Calling adapt() automatically." << std::endl;
        adapt();
      }
      if (adaptationStep_==postAdaptDone)
      {
        std::cout << "You already called call postAdapt() ! Aborting postAdapt()." << std::endl;
        return;
      }

            #ifdef SUBGRID_VERBOSE
      std::cout << "postadapt 1 (start)" << std::endl;
            #endif

      for (int level=0; level<=maxLevel(); ++level)
        refinementMark_[level].unsetAll();

      if (hostAdaptationStep_ == adaptDone)
      {
        hostgrid_->postAdapt();

        hostAdaptationStep_ = postAdaptDone;

                #ifdef SUBGRID_VERBOSE
        std::cout << "postadapt 2 (host grid postadapt called)" << std::endl;
                #endif
      }

      adaptationStep_ = postAdaptDone;

      return;
    }

    /*@}*/
#endif

    /** \brief Size of the overlap on the leaf level */
    unsigned int overlapSize(int codim) const {
      return hostgrid_->overlapSize(codim);
    }


    /** \brief Size of the ghost cell layer on the leaf level */
    unsigned int ghostSize(int codim) const {
      return hostgrid_->ghostSize(codim);
    }


    /** \brief Size of the overlap on a given level */
    unsigned int overlapSize(int level, int codim) const {
      return hostgrid_->overlapSize(level,codim);
    }


    /** \brief Size of the ghost cell layer on a given level */
    unsigned int ghostSize(int level, int codim) const {
      return hostgrid_->ghostSize(level,codim);
    }


#if 0
    /** \brief Distributes this grid over the available nodes in a distributed machine
     *
     * \param minlevel The coarsest grid level that gets distributed
     * \param maxlevel does currently get ignored
     */
    void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
      DUNE_THROW(NotImplemented, "IdentityGrid::loadBalance()");
    }

    /** \brief The communication interface
     *  @param T: array class holding data associated with the entities
     *  @param P: type used to gather/scatter data in and out of the message buffer
     *  @param codim: communicate entites of given codim
     *  @param if: one of the predifined interface types, throws error if it is not implemented
     *  @param level: communicate for entities on the given level
     *
     *  Implements a generic communication function sending an object of type P for each entity
     *  in the intersection of two processors. P has two methods gather and scatter that implement
     *  the protocol. Therefore P is called the "protocol class".
     */
    template<class T, template<class> class P, int codim>
    void communicate (T& t, InterfaceType iftype, CommunicationDirection dir, int level);

    /*! The new communication interface

       communicate objects for all codims on a given level
     */
    template<class DataHandle>
    void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
    {}

    template<class DataHandle>
    void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
    {}
#endif


    /** dummy collective communication */
    const CollectiveCommunication<IdentityGrid>& comm () const
    {
      return ccobj;
    }


    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    //! Returns the hostgrid this IdentityGrid lives in
    HostGridType& getHostGrid() const
    {
      return *hostgrid_;
    }


    //! Returns the hostgrid entity encapsulated in given subgrid entity
    template <int codim>
    typename HostGrid::Traits::template Codim<codim>::EntityPointer getHostEntity(const typename Traits::template Codim<codim>::Entity& e) const
    {
      return getRealImplementation(e).hostEntity_;
    }


    /** \brief Track hostgrid adaptation
     *
     * Returns true if *adapt() method of hostgrid was called
     * during last run of *adapt() in subgrid
     */
    bool hostGridAdapted ()
    {
      return (hostAdaptationStep_!=nothingDone);
    }

  protected:

    //! The host grid which contains the actual grid hierarchy structure
    HostGrid* hostgrid_;

  private:

    //! compute the grid indices and ids
    void setIndices()
    {
      localIdSet_.update();

      globalIdSet_.update();

      // //////////////////////////////////////////
      //   Create the index sets
      // //////////////////////////////////////////
      for (int i=levelIndexSets_.size(); i<=maxLevel(); i++) {
        IdentityGridLevelIndexSet<const IdentityGrid<HostGrid> >* p
          = new IdentityGridLevelIndexSet<const IdentityGrid<HostGrid> >();
        levelIndexSets_.push_back(p);
      }

      for (int i=0; i<=maxLevel(); i++)
        if (levelIndexSets_[i])
          levelIndexSets_[i]->update(*this, i);

      leafIndexSet_.update(*this);

    }

    //! \todo Please doc me !
    CollectiveCommunication<IdentityGrid> ccobj;

    //! Our set of level indices
    std::vector<IdentityGridLevelIndexSet<const IdentityGrid<HostGrid> >*> levelIndexSets_;

    //! \todo Please doc me !
    IdentityGridLeafIndexSet<const IdentityGrid<HostGrid> > leafIndexSet_;

    //! \todo Please doc me !
    IdentityGridGlobalIdSet<const IdentityGrid<HostGrid> > globalIdSet_;

    //! \todo Please doc me !
    IdentityGridLocalIdSet<const IdentityGrid<HostGrid> > localIdSet_;

    //! Stores the maximal difference of levels than elements containing a common vertex should have
    int maxLevelDiff_;

    //! Defines the possible steps in adaptation cycle
    enum AdaptationStep {nothingDone, preAdaptDone, adaptDone, postAdaptDone};

    //! Stores current adaptation step of subgrid
    AdaptationStep adaptationStep_;

    //! Stores current adaptation step of host grid
    AdaptationStep hostAdaptationStep_;

  }; // end Class IdentityGrid




  namespace Capabilities
  {
    //! \todo Please doc me !
    template<class HostGrid, int codim>
    struct hasEntity< IdentityGrid<HostGrid>, codim>
    {
      static const bool v = hasEntity<HostGrid,codim>::v;
    };


    //! \todo Please doc me !
    template<class HostGrid>
    struct isParallel< IdentityGrid<HostGrid> >
    {
      static const bool v = isParallel<HostGrid>::v;
    };


    //! \todo Please doc me !
    template<class HostGrid>
    struct hasHangingNodes< IdentityGrid<HostGrid> >
    {
      static const bool v = hasHangingNodes<HostGrid>::v;
    };

    //! \todo Please doc me !
    template<class HostGrid>
    struct isLevelwiseConforming< IdentityGrid<HostGrid> >
    {
      static const bool v = isLevelwiseConforming<HostGrid>::v;
    };

    //! \todo Please doc me !
    template<class HostGrid>
    struct isLeafwiseConforming< IdentityGrid<HostGrid> >
    {
      static const bool v = isLeafwiseConforming<HostGrid>::v;
    };
  }

} // namespace Dune

#endif
