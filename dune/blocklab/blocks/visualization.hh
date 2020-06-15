#ifndef DUNE_BLOCKLAB_BLOCKS_VISUALIZATION_HH
#define DUNE_BLOCKLAB_BLOCKS_VISUALIZATION_HH

#include<dune/blocklab/blocks/blockbase.hh>
#include<dune/common/parametertree.hh>
#include<dune/grid/io/file/vtk.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>

#include<filesystem>
#include<memory>
#include<variant>


namespace Dune::BlockLab {

  namespace impl {

    template<typename GV, bool instationary>
    struct VTKWriterChooser
    {
      using type = Dune::SubsamplingVTKWriter<GV>;
      using ptype = std::shared_ptr<type>;
    };


    template<typename GV>
    struct VTKWriterChooser<GV, true>
    {
      using type = Dune::VTKSequenceWriter<GV>;
      using ptype = std::shared_ptr<type>;
    };

  }

  template<typename P, typename V, std::size_t i>
  class VisualizationBlock
    : public ParentBlockBase<P, V>
  {
    public:
    using Traits = BlockTraits<P, V, i>;

    using VTKWriter = std::variant<typename impl::VTKWriterChooser<typename Traits::GridView, true>::ptype,
                                   typename impl::VTKWriterChooser<typename Traits::GridView, false>::ptype>;

    template<typename Context>
    VisualizationBlock(Context& ctx, const Dune::ParameterTree& config)
      : ParentBlockBase<P, V>(ctx, config)
      , instationary(config.get<bool>("instationary", false))
      , intervals(config.get<int>("intervals", 1))
      , time(0.0)
      , name(config.get<std::string>("filename"))
      , path(config.get<std::string>("path", ""))
    {}

    VisualizationBlock(bool instationary,
		       int intervals,
                       const std::string& name,
                       const std::string& path)
      : instationary(instationary)
      , intervals(intervals)
      , time(0.0)
      , name(name)
      , path(path)
    {}

    virtual ~VisualizationBlock() = default;

    virtual void update_parameter(std::string name, typename Traits::Parameter param) override
    {
      if (name == "time")
        time = std::get<double>(param);

      for (auto block : this->blocks)
        block->update_parameter(name, param);
    }

    virtual void setup() override
    {
      // Make sure that the output directory exists
      std::filesystem::create_directory(std::filesystem::current_path().append(path));

      // Instantiate the VTK writer instance
      auto vector = this->solver->template getVector<i>();
      auto gv = vector->gridFunctionSpace().gridView();

      if (instationary)
      {
         auto inner = std::make_shared<typename impl::VTKWriterChooser<typename Traits::GridView, false>::type>(gv, Dune::RefinementIntervals(intervals));
         vtkwriter = std::make_shared<typename impl::VTKWriterChooser<typename Traits::GridView, true>::type>(inner, name, path, "");
      }
      else
         vtkwriter = std::make_shared<typename impl::VTKWriterChooser<typename Traits::GridView, false>::type>(gv, Dune::RefinementIntervals(intervals));

      // Recrusively call setup on child block
      for(auto block : this->blocks)
	block->setup();
    }

    virtual void apply() override
    {
      for (auto block: this->blocks)
        block->apply();

      if (instationary)
        std::get<typename impl::VTKWriterChooser<typename Traits::GridView, true>::ptype>(vtkwriter)->write(time, Dune::VTK::appendedraw);
      else
        std::get<typename impl::VTKWriterChooser<typename Traits::GridView, false>::ptype>(vtkwriter)->write(name, Dune::VTK::ascii);
    }

    template<typename Container>
    void add_dataset(std::shared_ptr<Container> container)
    {
      if (instationary)
        Dune::PDELab::addSolutionToVTKWriter(
            *(std::get<typename impl::VTKWriterChooser<typename Traits::GridView, true>::ptype>(vtkwriter)),
            container->gridFunctionSpaceStorage(),
            container);
      else
        Dune::PDELab::addSolutionToVTKWriter(
            *(std::get<typename impl::VTKWriterChooser<typename Traits::GridView, false>::ptype>(vtkwriter)),
            container->gridFunctionSpaceStorage(),
            container);
    };

    template<typename Container>
    void add_celldata(std::shared_ptr<Container> container, std::string name)
    {
      if (instationary)
        std::get<typename impl::VTKWriterChooser<typename Traits::GridView, true>::ptype>(vtkwriter)->addCellData(*container, name);
      else
        std::get<typename impl::VTKWriterChooser<typename Traits::GridView, false>::ptype>(vtkwriter)->addCellData(*container, name);
    }

    private:
    bool instationary;
    int intervals;
    double time;
    std::string name;
    std::string path;
    VTKWriter vtkwriter;
  };


  template<typename P, typename V, std::size_t i>
  class VectorVisualizationBlock
    : public BlockBase<P, V, i>
  {
    public:
    using Traits = BlockTraits<P, V, i>;

    template<typename Context>
    VectorVisualizationBlock(Context&, const Dune::ParameterTree&)
    {}

    virtual ~VectorVisualizationBlock() = default;

    virtual void setup() override
    {
      auto vector = this->solver->template getVector<i>();
      std::dynamic_pointer_cast<VisualizationBlock<P, V, i>>(this->parent)->add_dataset(vector);
    }
  };

} // namespace Dune::BlockLab

#endif
