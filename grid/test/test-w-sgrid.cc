// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/sgrid.hh>
#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkintersectionit.cc>


#include <dune/identitygrid/identitygrid.hh>



using namespace Dune;


// test subgrid for given dimension
template <int dim>
void testDim()
{
  typedef SGrid<dim,dim> GridType;
  int n[dim];
  double h[dim];

  for (int i=0; i<dim; ++i)
  {
    n[i] = 1;
    h[i] = 1.0;
  }

  GridType grid(n,h);

  grid.globalRefine(1);

  typedef typename GridType::template Codim<0>::LevelIterator HostElementIterator;
  typedef typename IdentityGrid<GridType>::template Codim<0>::LevelIterator ElementIterator;

  IdentityGrid<GridType> subGrid(grid);

  gridcheck(subGrid);
  checkIntersectionIterator(subGrid);

  // refine subgrid
  {
    ElementIterator eIt    = subGrid.template lbegin<0>(0);
    ElementIterator eEndIt = subGrid.template lend<0>(0);
    for (; eIt!=eEndIt; ++eIt) {
      if (eIt->geometry()[0][0] < 0.25)
        subGrid.mark(1, eIt);
    }

    subGrid.preAdapt();
    subGrid.adapt();
    subGrid.postAdapt();
  }

  // check locally refined subgrid
  std::cout << "Check locally refined subgrid with dim=" << dim << std::endl;

  gridcheck(subGrid);
  checkIntersectionIterator(subGrid);

}


int main (int argc, char *argv[]) try
{
  testDim<1>();
  testDim<2>();
  testDim<3>();
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Exception e) {

  std::cout << e << std::endl;
  return 1;
}
