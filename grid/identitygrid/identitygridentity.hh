// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_IDENTITYGRIDENTITY_HH
#define DUNE_IDENTITYGRIDENTITY_HH

/** \file
 * \brief The IdentityGridEntity class
 */

#include <dune/grid/common/referenceelements.hh>


namespace Dune {


  // Forward declarations

  template<int codim, int dim, class GridImp>
  class IdentityGridEntity;

  template<int codim, class GridImp>
  class IdentityGridEntityPointer;

  template<int codim, PartitionIteratorType pitype, class GridImp>
  class IdentityGridLevelIterator;

  template<class GridImp>
  class IdentityGridLevelIntersectionIterator;

  template<class GridImp>
  class IdentityGridLeafIntersectionIterator;

  template<class GridImp>
  class IdentityGridHierarchicIterator;




  template<int codim, int dim, class GridImp>
  class IdentityGridMakeableEntity :
    public GridImp::template Codim<codim>::Entity
  {
  public:

    // The codimension of this entitypointer wrt the host grid
    enum {CodimInHostGrid = GridImp::HostGridType::dimension - GridImp::dimension + codim};

    // EntityPointer to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::EntityPointer HostGridEntityPointer;


    //! \todo Please doc me !
    IdentityGridMakeableEntity(const GridImp* identityGrid, const HostGridEntityPointer& hostEntity) :
      GridImp::template Codim<codim>::Entity (IdentityGridEntity<codim, dim, const GridImp>(identityGrid,hostEntity)),
      identityGrid_(identityGrid)
    {}


    //! \todo Please doc me !
    void setToTarget(const HostGridEntityPointer& hostEntity) {
      this->realEntity.setToTarget(hostEntity);
    }


    //! \todo Please doc me !
    const HostGridEntityPointer& getTarget() {
      return this->realEntity.hostEntity_;
    }


  private:

    const GridImp* identityGrid_;
  };


  //**********************************************************************
  //
  // --IdentityGridEntity
  // --Entity
  //
  /** \brief The implementation of entities in a IdentityGrid
   *   \ingroup IdentityGrid
   *
   *  A Grid is a container of grid entities. An entity is parametrized by the codimension.
   *  An entity of codimension c in dimension d is a d-c dimensional object.
   *
   */
  template<int codim, int dim, class GridImp>
  class IdentityGridEntity :
    public EntityDefaultImplementation <codim,dim,GridImp,IdentityGridEntity>
  {
    friend class IdentityGridMakeableEntity<codim,dim,GridImp>;

    template <class GridImp_>
    friend class IdentityGridLevelIndexSet;

    template <class GridImp_>
    friend class IdentityGridLeafIndexSet;

    template <class GridImp_>
    friend class IdentityGridLocalIdSet;

    template <class GridImp_>
    friend class IdentityGridGlobalIdSet;

    friend class IdentityGridEntityPointer<codim,GridImp>;


  private:

    typedef typename GridImp::ctype ctype;

    // The codimension of this entitypointer wrt the host grid
    enum {CodimInHostGrid = GridImp::HostGridType::dimension - GridImp::dimension + codim};

    // EntityPointer to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::EntityPointer HostGridEntityPointer;


  public:

    typedef typename GridImp::template Codim<codim>::Geometry Geometry;


    //! Constructor for an entity in a given grid level
    IdentityGridEntity(const GridImp* identityGrid, const HostGridEntityPointer& hostEntity) :
      hostEntity_(hostEntity),
      identityGrid_(identityGrid),
      geo_(0),
      geoInFather_(0)
    {}


    //! \todo Please doc me !
    IdentityGridEntity(const IdentityGridEntity& original) :
      hostEntity_(original.hostEntity_),
      identityGrid_(original.identityGrid_),
      geo_(0),
      geoInFather_(0)
    {}


    //! Destructor
    ~IdentityGridEntity()
    {
      if (geo_!=0)
      {
        delete geo_;
        geo_ = 0;
      }
      if (geoInFather_!=0)
      {
        delete geoInFather_;
        geoInFather_ = 0;
      }
    }


    //! \todo Please doc me !
    IdentityGridEntity& operator=(const IdentityGridEntity& original)
    {
      if (this != &original)
      {
        if (geo_!=0)
        {
          delete geo_;
          geo_ = 0;
        }
        if (geoInFather_!=0)
        {
          delete geoInFather_;
          geoInFather_ = 0;
        }
        identityGrid_ = original.identityGrid_;
        hostEntity_ = original.hostEntity_;
      }
      return *this;
    }


    //! level of this element
    int level () const {
      return hostEntity_->level();
    }


    /** \brief The partition type for parallel computing
     */
    PartitionType partitionType () const {
      return hostEntity_->partitionType();
    }


    /** Intra-element access to entities of codimension cc > codim. Return number of entities
     * with codimension cc.
     */
    template<int cc> int count () const {
      return hostEntity_->template count<cc>();
    }


    //! geometry of this entity
    const Geometry& geometry () const
    {
      if (geo_==0)
        geo_ = new MakeableInterfaceObject<Geometry>(hostEntity_->geometry());
      return *geo_;
    }


    HostGridEntityPointer hostEntity_;


  private:

    //! \todo Please doc me !
    void setToTarget(const HostGridEntityPointer& target)
    {
      if(geo_!=0)
      {
        delete geo_;
        geo_ = 0;
      }
      if (geoInFather_!=0)
      {
        delete geoInFather_;
        geoInFather_ = 0;
      }
      hostEntity_ = target;
    }


    const GridImp* identityGrid_;

    //! the current geometry
    mutable MakeableInterfaceObject<Geometry> *geo_;
    mutable MakeableInterfaceObject<Geometry> *geoInFather_;
  };




  //***********************
  //
  //  --IdentityGridEntity
  //
  //***********************
  /** \brief Specialization for codim-0-entities.
   * \ingroup IdentityGrid
   *
   * This class embodies the topological parts of elements of the grid.
   * It has an extended interface compared to the general entity class.
   * For example, Entities of codimension 0  allow to visit all neighbors.
   */
  template<int dim, class GridImp>
  class IdentityGridEntity<0,dim,GridImp> :
    public EntityDefaultImplementation<0,dim,GridImp, IdentityGridEntity>
  {
  public:

    // The codimension of this entitypointer wrt the host grid
    enum {CodimInHostGrid = GridImp::HostGridType::dimension - GridImp::dimension};

    // EntityPointer to the equivalent entity in the host grid
    typedef typename GridImp::HostGridType::Traits::template Codim<CodimInHostGrid>::EntityPointer HostGridEntityPointer;

    typedef typename GridImp::template Codim<0>::Geometry Geometry;

    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    //! The Iterator over intersections on this level
    typedef IdentityGridLevelIntersectionIterator<GridImp> LevelIntersectionIterator;

    //! The Iterator over intersections on the leaf level
    typedef IdentityGridLeafIntersectionIterator<GridImp> LeafIntersectionIterator;

    //! Iterator over descendants of the entity
    typedef IdentityGridHierarchicIterator<GridImp> HierarchicIterator;


    //! Constructor for an entity in a given grid level
    IdentityGridEntity(const GridImp* identityGrid, const HostGridEntityPointer& hostEntity) :
      identityGrid_(identityGrid),
      geo_(0),
      geoInFather_(0),
      hostEntity_(hostEntity)
    {}


    //! \todo Please doc me !
    IdentityGridEntity(const IdentityGridEntity& original) :
      identityGrid_(original.identityGrid_),
      geo_(0),
      geoInFather_(0),
      hostEntity_(original.hostEntity_)
    {}


    //! Destructor
    ~IdentityGridEntity()
    {
      if (geo_!=0)
      {
        delete geo_;
        geo_ = 0;
      }
      if (geoInFather_!=0)
      {
        delete geoInFather_;
        geoInFather_ = 0;
      }
    }


    //! \todo Please doc me !
    IdentityGridEntity& operator=(const IdentityGridEntity& original)
    {
      if (this != &original)
      {
        if (geo_!=0)
        {
          delete geo_;
          geo_ = 0;
        }
        if (geoInFather_!=0)
        {
          delete geoInFather_;
          geoInFather_ = 0;
        }
        identityGrid_ = original.identityGrid_;
        hostEntity_ = original.hostEntity_;
      }
      return *this;
    }


    //! Level of this element
    int level () const
    {
      return hostEntity_->level();
    }


    /** \brief The partition type for parallel computing */
    PartitionType partitionType () const {
      return hostEntity_->partitionType();
    }


    //! Geometry of this entity
    const Geometry& geometry () const
    {
      if (geo_==0)
      {
        geo_ = new MakeableInterfaceObject<Geometry>(hostEntity_->geometry());
      }
      return *geo_;
    }


    /** \brief Return the number of subEntities of codimension cc.
     */
    template<int cc>
    int count () const
    {
      return hostEntity_->template count<cc>();
    }


    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... count<cc>()-1
     */
    template<int cc>
    typename GridImp::template Codim<cc>::EntityPointer subEntity (int i) const {
      return IdentityGridEntityPointer<cc,GridImp>(identityGrid_, hostEntity_->template subEntity<cc>(i));
    }


    //! First level intersection
    IdentityGridLevelIntersectionIterator<GridImp> ilevelbegin () const {
      return IdentityGridLevelIntersectionIterator<GridImp>(identityGrid_,
                                                            hostEntity_->ilevelbegin());
    }


    //! Reference to one past the last neighbor
    IdentityGridLevelIntersectionIterator<GridImp> ilevelend () const {
      return IdentityGridLevelIntersectionIterator<GridImp>(identityGrid_,
                                                            hostEntity_->ilevelend());
    }


    //! First leaf intersection
    IdentityGridLeafIntersectionIterator<GridImp> ileafbegin () const {
      return IdentityGridLeafIntersectionIterator<GridImp>(identityGrid_,
                                                           hostEntity_->ileafbegin());
    }


    //! Reference to one past the last leaf intersection
    IdentityGridLeafIntersectionIterator<GridImp> ileafend () const {
      return IdentityGridLeafIntersectionIterator<GridImp>(identityGrid_,
                                                           hostEntity_->ileafend());
    }


    //! returns true if Entity has NO children
    bool isLeaf() const {
      return hostEntity_->isLeaf();
    }


    //! Inter-level access to father element on coarser grid.
    //! Assumes that meshes are nested.
    IdentityGridEntityPointer<0,GridImp> father () const {
      return IdentityGridEntityPointer<0,GridImp>(identityGrid_, hostEntity_->father());
    }


    /** \brief Location of this element relative to the reference element element of the father.
     * This is sufficient to interpolate all dofs in conforming case.
     * Nonconforming may require access to neighbors of father and
     * computations with local coordinates.
     * On the fly case is somewhat inefficient since dofs  are visited several times.
     * If we store interpolation matrices, this is tolerable. We assume that on-the-fly
     * implementation of numerical algorithms is only done for simple discretizations.
     * Assumes that meshes are nested.
     */
    const LocalGeometry& geometryInFather () const {
      if (geoInFather_==0)
        geoInFather_ = new MakeableInterfaceObject<LocalGeometry>(hostEntity_->geometryInFather());
      return *geoInFather_;
    }


    /** \brief Inter-level access to son elements on higher levels<=maxlevel.
     * This is provided for sparsely stored nested unstructured meshes.
     * Returns iterator to first son.
     */
    IdentityGridHierarchicIterator<GridImp> hbegin (int maxLevel) const
    {
      return IdentityGridHierarchicIterator<const GridImp>(identityGrid_, *this, maxLevel);
    }


    //! Returns iterator to one past the last son
    IdentityGridHierarchicIterator<GridImp> hend (int maxLevel) const
    {
      return IdentityGridHierarchicIterator<const GridImp>(identityGrid_, *this, maxLevel, true);
    }


    //! \todo Please doc me !
    bool wasRefined () const
    {
      if (identityGrid_->adaptationStep!=GridImp::adaptDone)
        return false;

      int level = this->level();
      int index = identityGrid_->levelIndexSet(level).index(*this);
      return identityGrid_->refinementMark_[level][index];
    }


    //! \todo Please doc me !
    bool mightBeCoarsened () const
    {
      return true;
    }


    // /////////////////////////////////////////
    //   Internal stuff
    // /////////////////////////////////////////


    //! \todo Please doc me !
    void setToTarget(const HostGridEntityPointer& target)
    {
      if(geo_!=0)
      {
        delete geo_;
        geo_ = 0;
      }
      if (geoInFather_!=0)
      {
        delete geoInFather_;
        geoInFather_ = 0;
      }
      hostEntity_ = target;
    }


    const GridImp* identityGrid_;

    //! the current geometry
    mutable MakeableInterfaceObject<Geometry> *geo_;

    //! \todo Please doc me !
    mutable MakeableInterfaceObject<LocalGeometry> *geoInFather_;

    //! \todo Please doc me !
    HostGridEntityPointer hostEntity_;


  private:

    typedef typename GridImp::ctype ctype;

  }; // end of IdentityGridEntity codim = 0


} // namespace Dune


#endif
