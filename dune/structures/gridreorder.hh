#ifndef DUNE_STRUCTURES_GRIDREORDER_HH
#define DUNE_STRUCTURES_GRIDREORDER_HH

/** A tool to reorder the elements in a grid according to a sorting criterion.
 *
 * This was introduced as a debugging tool for the fibre-reinforced operator
 * to ensure that indices across the fibre curve are monotonuous and therefore
 * skeleton terms have no inside/outside confusion.
 */

#include<dune/common/parametertree.hh>
#include<dune/grid/common/gridfactory.hh>

#include<map>
#include<tuple>
#include<vector>


// A comparison object to define a total ordering on Dune::FieldVectors for the purpose
// of below algorithm to identify subsurface boundary faces with surface cells.
struct LexicographicFieldVectorComparator
{
  template<typename T, int dim>
  bool operator() (const Dune::FieldVector<T, dim>& a, const Dune::FieldVector<T, dim>& b) const
  {
    const T tolerance = 1e-6;
    using std::abs;
    for (int i=0; i<dim; ++i)
      if (abs(a[i] - b[i]) > tolerance)
        return a[i] < b[i];

    // If we ever make it here, a and b are fuzzy equal. In that case we do not really care about ordering.
    return true;
  }
};


template<typename Grid, typename ES, typename Physical, typename Compare>
auto reorder_impl(Grid&& grid, ES&& es, Physical&& physical, Compare&& compare)
{
  std::cout << "Starting a reordering strategy..." << std::endl;
  static constexpr int dim = std::decay_t<ES>::Traits::dimension;

  // Get a vector of EntitySeed objects
  std::cout << "Creating a vector of entity seeds..." << std::endl;
  std::vector<typename std::decay_t<ES>::Traits::Element::EntitySeed> seeds;
  for (auto e : elements(es))
    seeds.push_back(e.seed());

  // Sort the seeds by the given comparator
  std::cout << "Sort this vector..." << std::endl;
  std::sort(seeds.begin(), seeds.end(), compare);

  // Construct a new grid through a GridFactory
  std::cout << "Initialize a GridFactory with all vertices..." << std::endl;
  Dune::GridFactory<typename std::decay_t<ES>::Traits::Grid> factory;
  auto& is = es.indexSet();

  // Add all the vertices to the grid
  {
	 std::vector<typename std::decay_t<ES>::Traits::Element::Geometry::GlobalCoordinate> vertex_positions(es.size(dim));

     for (auto vertex : vertices(es))
	   vertex_positions[is.index(vertex)] = vertex.geometry().corner(0);

     for (auto pos: vertex_positions)
	   factory.insertVertex(pos);
  }

  std::cout << "Adding elements to the GridFactory..." << std::endl;
  for (auto seed : seeds)
  {
	std::vector<unsigned int> vertex_indices;
	auto entity = grid->entity(seed);
	for(std::size_t i=0; i<entity.subEntities(dim); ++i)
	  vertex_indices.push_back(is.index(entity.template subEntity<dim>(i)));

	factory.insertElement(entity.geometry().type(), vertex_indices);
  }

  std::cout << "Creating the reordered grid..." << std::endl;
  auto newgrid = factory.createGrid();

  return std::make_tuple(std::decay_t<Grid>(newgrid),
		                 std::decay_t<ES>(newgrid->leafGridView()),
						 physical);
}


template<typename Grid, typename ES, typename Physical>
auto reorder_grid(Grid&& grid, ES&& es, Physical&& physical, const Dune::ParameterTree& params)
{
  auto strategy = params.get<std::string>("reorder", "none");

  if (strategy == "lefttoright")
	return reorder_impl(grid, es, physical,
	    [grid](const auto& seed1, const auto& seed2)
	  	{
	  	  auto center1 = grid->entity(seed1).geometry().center();
	  	  auto center2 = grid->entity(seed2).geometry().center();
	  	  return LexicographicFieldVectorComparator()(center1, center2);
	    });

  if (strategy == "righttoleft")
    return reorder_impl(grid, es, physical,
        [grid](const auto& seed1, const auto& seed2)
        {
          auto center1 = grid->entity(seed1).geometry().center();
          auto center2 = grid->entity(seed2).geometry().center();
          return !LexicographicFieldVectorComparator()(center1, center2);
         });

  if (strategy != "none")
	DUNE_THROW(Dune::Exception, "Unknown Grid reordering strategy!");

  return std::make_tuple(grid, es, physical);
}

#endif
