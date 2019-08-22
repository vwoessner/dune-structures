#ifndef DUNE_STRUCTURES_VISUALIZATION_HH
#define DUNE_STRUCTURES_VISUALIZATION_HH

/** Visualization tools */

template<typename GV>
struct RankDummyContainer
{
  DummyContainer(Dune::MPIHelper& helper, GV gv) : rank(helper.rank()), size_(gv.size(0))
  {}

  double operator[](std::size_t i) const
  {
    return rank;
  }

  std::size_t size() const
  {
    return size_;
  }

  double rank;
  std::size_t size_;
};


/** Writes the MPI rank as a field into a VTK file */
template<typename GV>
void write_rankdata(Dune::VTKWriter<GV>& vtkwriter, Dune::MPIHelper& helper, GV gv)
{
  vtkwriter.addCellData(RankDummyContainer(helper, gv), "mpirank");
}

#endif
